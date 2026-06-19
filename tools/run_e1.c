/* run_e1 — experiment E1 harness: sheet detection -> signed affinity -> Mutex
 * Watershed -> segmentation metrics, on a synthetic concentric-shell volume that
 * stands in for a scroll cross-section (cylindrical "wraps" with a controllable
 * inter-wrap gap). Validates the whole chain and, crucially, whether the signed
 * partition keeps adjacent wraps in separate segments (low voi_merge).
 *
 * Swap the synthetic generator for nrrd_read_f32(<cube>) + the instance labels to
 * run E1 on real Vesuvius data.
 *
 *   run_e1 [--gap G] [--thick T] [--noise N] [--kattract A] [--krepel R]
 *          [--out PREFIX]
 */
#include "io/nrrd.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/affinity.h"
#include "segmentation/partition.h"
#include "eval/metrics.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static u32 lcg = 12345u;
static double frand(void) { lcg = lcg * 1664525u + 1013904223u; return (lcg >> 8) * (1.0 / 16777216.0); }
static double gauss(void) {  // Box-Muller
  double u1 = frand() + 1e-9, u2 = frand();
  return sqrt(-2.0 * log(u1)) * cos(6.283185307 * u2);
}

int main(int argc, char **argv) {
  int N = 96;
  double gap = 2.0, thick = 2.0, noise = 12.0;
  affinity_params ap = affinity_default_params();
  const char *out = NULL;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--gap") && i + 1 < argc) gap = atof(argv[++i]);
    else if (!strcmp(argv[i], "--thick") && i + 1 < argc) thick = atof(argv[++i]);
    else if (!strcmp(argv[i], "--noise") && i + 1 < argc) noise = atof(argv[++i]);
    else if (!strcmp(argv[i], "--kattract") && i + 1 < argc) ap.k_attract = (f32)atof(argv[++i]);
    else if (!strcmp(argv[i], "--krepel") && i + 1 < argc) ap.k_repel = (f32)atof(argv[++i]);
    else if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
  }

  int nz = N, ny = N, nx = N;
  size_t n = (size_t)nz * ny * nx, nynx = (size_t)ny * nx;
  double c = N / 2.0, R0 = 12.0, spacing = thick + gap;
  int K = 0;
  double Rk[32];
  for (double R = R0; R < c - 4.0 && K < 32; R += spacing) Rk[K++] = R;

  f32 *vol = (f32 *)malloc(n * sizeof(f32));
  u32 *gt = (u32 *)malloc(n * sizeof(u32));
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = (size_t)z * nynx + (size_t)y * nx + x;
        double r = sqrt((y - c) * (y - c) + (x - c) * (x - c));
        int lab = 0;
        for (int k = 0; k < K; k++)
          if (fabs(r - Rk[k]) <= thick / 2.0) { lab = k + 1; break; }
        gt[i] = (u32)lab;
        double v = (lab ? 220.0 : 30.0) + noise * gauss();
        vol[i] = (f32)(v < 0 ? 0 : v > 255 ? 255 : v);
      }

  printf("synthetic: %d^3, %d shells, thick=%.1f gap=%.1f noise=%.1f\n", N, K, thick, gap, noise);

  // 1) sheet detection
  f32 *sheet = (f32 *)malloc(n * sizeof(f32));
  f32 *norm = (f32 *)malloc(3 * n * sizeof(f32));
  st_params sp = st_default_params();
  if (st_sheet_detect(vol, nz, ny, nx, &sp, sheet, norm)) { fprintf(stderr, "detect oom\n"); return 1; }

  // 2) foreground mask (intensity gate) + signed affinity graph
  u8 *mask = (u8 *)malloc(n);
  size_t fg = 0;
  for (size_t i = 0; i < n; i++) { mask[i] = vol[i] > 128.0f; fg += mask[i]; }

  sgraph g;
  s32 *node_of;
  u32 *voxel_of;
  if (affinity_build_voxel(mask, sheet, norm, nz, ny, nx, &ap, &g, &node_of, &voxel_of)) {
    fprintf(stderr, "affinity oom\n"); return 1;
  }
  printf("graph: %d nodes (fg=%zu), %d edges  [k_attract=%.2f k_repel=%.2f]\n",
         g.nnodes, fg, g.nedges, ap.k_attract, ap.k_repel);

  // 3) Mutex Watershed partition
  u32 *clab = (u32 *)malloc((size_t)(g.nnodes ? g.nnodes : 1) * sizeof(u32));
  int ncl = mws_partition(&g, clab);
  printf("clusters: %d (vs %d ground-truth wraps)\n", ncl, K);

  // 4) assemble seg volume (bg=0, foreground = cluster+1) and score
  u32 *seg = (u32 *)calloc(n, sizeof(u32));
  for (int v = 0; v < g.nnodes; v++) seg[voxel_of[v]] = clab[v] + 1;
  eval_result r = eval_seg(seg, gt, n);
  printf("VOI: split=%.4f merge=%.4f total=%.4f bits   ARE=%.4f (P=%.4f R=%.4f)\n",
         r.voi_split, r.voi_merge, r.voi, r.are, r.precision, r.recall);
  printf("  -> voi_merge is the wrap-collapse signal (want ~0).\n");

  if (out) {
    char path[512];
    snprintf(path, sizeof path, "%s_sheetness.nrrd", out);
    nrrd_write_f32(path, sheet, nz, ny, nx);
    // seg as u8 (mod 256 just for visualization)
    u8 *segu8 = (u8 *)malloc(n);
    for (size_t i = 0; i < n; i++) segu8[i] = (u8)(seg[i] % 256);
    snprintf(path, sizeof path, "%s_seg.nrrd", out);
    nrrd_write_u8(path, segu8, nz, ny, nx);
    free(segu8);
    printf("wrote %s_sheetness.nrrd and %s_seg.nrrd\n", out, out);
  }

  sgraph_free(&g);
  free(vol); free(gt); free(sheet); free(norm); free(mask);
  free(node_of); free(voxel_of); free(clab); free(seg);
  return 0;
}
