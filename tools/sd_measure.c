/* sd_measure.c — comprehensive surface-detection measurement harness.
 *
 * Loads a CT cube + label, optionally PREPROCESSES the CT (histogram-valley
 * air-cut and/or percentile contrast stretch — the fysics-style front-end,
 * self-calibrated on THIS cube since we have no metadata.json), computes the
 * structure-tensor sheetness ONCE, then sweeps the binarization threshold and
 * reports the FULL competition composite for each operating point:
 *   pred voxels, Dice, SurfaceDice@2, VOI + VOI_score, TopoScore (native tiled),
 *   b0/b1/b2 (pred) vs b0/b1/b2 (gt), and the blended Score.
 *
 * Purpose: zero in on WHERE the score is lost (it's topology) and whether better
 * front-end preprocessing moves it. Usage:
 *   sd_measure IMG.tif LAB.tif [pre=none|aircut|stretch|both] [aircut_frac=auto]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "postproc/morph.h"
#include "eval/metrics.h"
#include "eval/score.h"
#include "fysics.h"

static int find_valley(const u8 *img, size_t n) {
  long h[256] = {0};
  for (size_t i = 0; i < n; i++) h[img[i]]++;
  // smooth the histogram (box r=3) then find the lowest bin between the two
  // dominant low/high modes — classic air<->material valley.
  double s[256];
  for (int v = 0; v < 256; v++) {
    long acc = 0; int c = 0;
    for (int k = -3; k <= 3; k++) { int j = v + k; if (j >= 0 && j < 256) { acc += h[j]; c++; } }
    s[v] = (double)acc / c;
  }
  // dark mode = argmax in [0,96], then valley = argmin from dark..dark+96
  int dark = 0; for (int v = 1; v < 96; v++) if (s[v] > s[dark]) dark = v;
  int valley = dark; double best = s[dark];
  for (int v = dark; v < dark + 120 && v < 240; v++) if (s[v] < best) { best = s[v]; valley = v; }
  return valley;
}

static void percentiles(const u8 *img, size_t n, double plo, double phi, int *lo, int *hi) {
  long h[256] = {0}; for (size_t i = 0; i < n; i++) h[img[i]]++;
  long target_lo = (long)(plo * n), target_hi = (long)(phi * n), acc = 0;
  *lo = 0; *hi = 255;
  for (int v = 0; v < 256; v++) { acc += h[v]; if (acc >= target_lo) { *lo = v; break; } }
  acc = 0;
  for (int v = 0; v < 256; v++) { acc += h[v]; if (acc >= target_hi) { *hi = v; break; } }
}

int main(int argc, char **argv) {
  setvbuf(stdout, NULL, _IOLBF, 0);
  if (argc < 3) { fprintf(stderr, "usage: %s IMG.tif LAB.tif [pre=none|aircut|stretch|both]\n", argv[0]); return 2; }
  const char *pre = argc > 3 ? argv[3] : "none";
  const char *field = argc > 4 ? argv[4] : "st";   // st | gated
  int do_aircut = strstr(pre, "aircut") || strstr(pre, "both");
  int do_stretch = strstr(pre, "stretch") || strstr(pre, "both");
  int do_denoise = strstr(pre, "denoise") || strstr(pre, "full") != NULL;
  int do_deconv  = strstr(pre, "deconv")  || strstr(pre, "full") != NULL;
  int do_gate = strcmp(field, "gated") == 0;

  int nz, ny, nx, lz, ly, lx;
  u8 *img = tiff_load_u8(argv[1], &nz, &ny, &nx);
  u8 *lab = tiff_load_u8(argv[2], &lz, &ly, &lx);
  if (!img || !lab) { fprintf(stderr, "load failed\n"); return 1; }
  size_t n = (size_t)nz * ny * nx;
  printf("cube %dx%dx%d  pre=%s  field=%s\n", nz, ny, nx, pre, field);

  // --- preprocessing (self-calibrated on this cube) -------------------------
  int valley = find_valley(img, n), plo = 0, phi = 255;
  if (do_stretch) percentiles(img, n, 0.02, 0.998, &plo, &phi);
  printf("calib: valley=%d  stretch[%d,%d]  aircut=%d stretch=%d\n", valley, plo, phi, do_aircut, do_stretch);
  long zeroed = 0;
  f32 *volf = (f32 *)malloc(n * sizeof(f32));
  for (size_t i = 0; i < n; i++) {
    int v = img[i];
    if (do_aircut && v <= valley) { volf[i] = 0.f; zeroed++; continue; }
    if (do_stretch) {
      double t = (double)(v - plo) / (phi - plo > 1 ? phi - plo : 1);
      if (t < 0) t = 0; if (t > 1) t = 1; v = (int)(t * 255 + 0.5);
    }
    volf[i] = (f32)v;
  }
  if (do_aircut) printf("aircut: zeroed %ld voxels (%.1f%%)\n", zeroed, 100.0 * zeroed / n);

  // --- fysics front-end: DECONV (deblur the PSF that breaks thin sheets) then
  // guided DENOISE (kill speckle that punches pinholes). Both self-calibrate on
  // this cube: noise_ref from fy_estimate_noise -> guided eps; PSF default ~1vox
  // (the net nabu blur per fysics docs) since we lack metadata.json. ----------
  if (do_deconv) {
    fy_physics p = {0};
    p.psf_sigma_vox = 1.0;   // no metadata: use the documented ~1-voxel net blur
    f32 *tmp = (f32 *)malloc(n * sizeof(f32));
    int rc = fy_deconvolve_matched(volf, tmp, nz, ny, nx, &p, /*tikhonov auto*/ 0.0);
    if (rc == 0) { memcpy(volf, tmp, n * sizeof(f32)); printf("deconv: matched-Wiener psf=1.0 OK\n"); }
    else printf("deconv: FAILED (rc=%d)\n", rc);
    free(tmp);
  }
  if (do_denoise) {
    fy_noise_model nm = {0};
    fy_estimate_noise(volf, nz, ny, nx, 5, 10.0, 100.0 /*ref intensity ~u8*/, &nm);
    double eps = fy_guided_eps_for_noise(nm.noise_ref);
    f32 *tmp = (f32 *)malloc(n * sizeof(f32));
    int rc = fy_guided_denoise(volf, tmp, nz, ny, nx, /*radius*/ 2, eps);
    if (rc == 0) { memcpy(volf, tmp, n * sizeof(f32)); printf("denoise: guided r=2 noise_ref=%.2f eps=%.3f OK\n", nm.noise_ref, eps); }
    else printf("denoise: FAILED (rc=%d)\n", rc);
    free(tmp);
  }

  // --- structure-tensor sheetness (once) ------------------------------------
  fprintf(stderr, "[bc] sheet detect start\n");
  f32 *sheet = (f32 *)malloc(n * sizeof(f32));
  if (st_sheet_detect(volf, nz, ny, nx, NULL, sheet, NULL) != 0) { fprintf(stderr, "sheet detect failed\n"); return 1; }
  fprintf(stderr, "[bc] sheet detect done\n");

  // --- CONTRAST GATE (field=gated): pure (l0-l1)/(l0+eps) is anisotropy only,
  // so it saturates in flat/noisy regions. Gate it by local edge strength so
  // sheetness requires a REAL edge, not just a directional bias. Edge strength =
  // |grad| of the (preprocessed) volume, smoothed, normalized to its p90. ------
  if (do_gate) {
    f32 *gm = (f32 *)malloc(n * sizeof(f32));
    #define IDX(zz,yy,xx) (((size_t)(zz)*ny+(yy))*nx+(xx))
    #pragma omp parallel for
    for (int z = 0; z < nz; z++) for (int y = 0; y < ny; y++) for (int x = 0; x < nx; x++) {
      int zp = z>0?z-1:z, zn = z<nz-1?z+1:z, yp=y>0?y-1:y, yn=y<ny-1?y+1:y, xp=x>0?x-1:x, xn=x<nx-1?x+1:x;
      float gx = volf[IDX(z,y,xn)] - volf[IDX(z,y,xp)];
      float gy = volf[IDX(z,yn,x)] - volf[IDX(z,yp,x)];
      float gz = volf[IDX(zn,y,x)] - volf[IDX(zp,y,x)];
      gm[IDX(z,y,x)] = sqrtf(gx*gx + gy*gy + gz*gz);
    }
    // p90 of |grad| via a coarse histogram (max gradient for u8 ~ 2*255*sqrt3)
    long gh[512] = {0}; float gmax = 0;
    for (size_t i = 0; i < n; i++) if (gm[i] > gmax) gmax = gm[i];
    float gscale = gmax > 0 ? 511.0f / gmax : 0;
    for (size_t i = 0; i < n; i++) { int b = (int)(gm[i]*gscale); if(b>511)b=511; gh[b]++; }
    long acc = 0, tgt = (long)(0.90 * n); float p90 = gmax;
    for (int b = 0; b < 512; b++) { acc += gh[b]; if (acc >= tgt) { p90 = (b+0.5f)/gscale; break; } }
    if (p90 < 1) p90 = 1;
    printf("contrast gate: |grad| p90=%.1f (max=%.1f)\n", p90, gmax);
    for (size_t i = 0; i < n; i++) {
      float g = gm[i] / p90; if (g > 1) g = 1;   // soft edge-strength gate in [0,1]
      sheet[i] *= g;
    }
    free(gm);
    #undef IDX
  }
  free(volf); free(img);

  // sheetness histogram (where does the signal live?)
  long sh[10] = {0};
  for (size_t i = 0; i < n; i++) { int b = (int)(sheet[i] * 10); if (b > 9) b = 9; if (b < 0) b = 0; sh[b]++; }
  printf("sheetness deciles:");
  for (int b = 0; b < 10; b++) printf(" %.0f.%d=%.0f%%", b / 1.0, 0, 100.0 * sh[b] / n);
  printf("\n\n");

  // --- FAST sweep (eval_surface only) to locate the operating point ---------
  printf("%-6s %9s %7s %7s  %5s %6s %5s  (gt b0=?)\n", "thr", "pred", "Dice", "SurfD", "b0", "b1", "b2");
  float thrs[] = {0.40f, 0.50f, 0.60f, 0.70f, 0.80f, 0.90f};
  int nthr = (int)(sizeof(thrs) / sizeof(thrs[0]));
  u8 *pred = (u8 *)malloc(n);
  u8 *best_pred = (u8 *)malloc(n);
  float best_thr = 0; double best_surfd = -1;
  for (int t = 0; t < nthr; t++) {
    float thr = thrs[t];
    size_t on = 0;
    for (size_t i = 0; i < n; i++) { pred[i] = sheet[i] >= thr; on += pred[i]; }
    remove_small_components(pred, nz, ny, nx, TOPO_CONN26, 1000);
    majority_filter(pred, pred, nz, ny, nx, 7);
    fill_holes(pred, nz, ny, nx);
    surface_eval e = eval_surface(pred, lab, nz, ny, nx, 2, 1, 2);
    printf("%-6.2f %9zu %7.3f %7.3f  %5ld %6ld %5ld  (gt b0=%ld)\n",
           thr, on, e.dice, e.surface_dice, e.pred_b0, e.pred_b1, e.pred_b2, e.gt_b0);
    fflush(stdout);
    if (e.surface_dice > best_surfd) { best_surfd = e.surface_dice; best_thr = thr; memcpy(best_pred, pred, n); }
  }

  // --- full official composite at the best operating point ------------------
  printf("\n=== full composite @ thr=%.2f (best SurfaceDice) ===\n", best_thr);
  comp_score s = competition_score(best_pred, lab, nz, ny, nx, 2, 1, 2);
  printf("SurfaceDice : %.4f\n", s.surface_dice);
  printf("VOI         : %.3f  -> VOI_score %.4f\n", s.voi, s.voi_score);
  printf("TopoScore   : %.4f\n", s.topo_score);
  printf("pred betti  : b0=%ld b1=%ld b2=%ld\n", s.pred_b0, s.pred_b1, s.pred_b2);
  printf("gt   betti  : b0=%ld b1=%ld b2=%ld\n", s.gt_b0, s.gt_b1, s.gt_b2);
  printf("SCORE       : %.4f  (0.30*topo + 0.35*surfD + 0.35*voi)\n", s.score);
  fflush(stdout);
  free(pred); free(best_pred); free(sheet); free(lab);
  return 0;
}
