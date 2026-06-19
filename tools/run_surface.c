/* run_surface.c — fully classical (no-ML) surface-detection baseline on the
 * Kaggle Vesuvius surface-detection data, end to end:
 *   load CT TIFF -> structure-tensor sheetness -> threshold -> postproc
 *   (dust removal -> iterated median/majority filter -> fill holes)
 *   -> ignore-aware score vs the label TIFF.
 *
 * This is the experiment that tests "did the detector need to be ML?" — see
 * docs/surface-detection-kaggle.md. Usage:
 *   run_surface IMAGE.tif LABEL.tif [thresh=0.5] [median_iters=7]
 *               [min_voxels=1000] [tol=2] [PRED_OUT.tif]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "postproc/morph.h"
#include "eval/metrics.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s IMAGE.tif LABEL.tif [thresh=0.5] [median_iters=7]"
                    " [min_voxels=1000] [tol=2]\n", argv[0]);
    return 2;
  }
  const char *img_path = argv[1], *lab_path = argv[2];
  float thresh   = argc > 3 ? (float)atof(argv[3]) : 0.5f;
  int   iters    = argc > 4 ? atoi(argv[4]) : 7;
  size_t min_vox = argc > 5 ? (size_t)atoll(argv[5]) : 1000;
  int   tol      = argc > 6 ? atoi(argv[6]) : 2;

  int nz, ny, nx, lz, ly, lx;
  u8 *img = tiff_load_u8(img_path, &nz, &ny, &nx);
  if (!img) { fprintf(stderr, "failed to load image %s\n", img_path); return 1; }
  u8 *lab = tiff_load_u8(lab_path, &lz, &ly, &lx);
  if (!lab) { fprintf(stderr, "failed to load label %s\n", lab_path); free(img); return 1; }
  if (lz != nz || ly != ny || lx != nx) {
    fprintf(stderr, "dim mismatch: image %dx%dx%d vs label %dx%dx%d\n", nz,ny,nx, lz,ly,lx);
    free(img); free(lab); return 1;
  }
  size_t n = (size_t)nz * ny * nx;
  printf("loaded %dx%dx%d (%zu voxels)\n", nz, ny, nx, n);

  // --- structure-tensor sheetness on the CT ---------------------------------
  f32 *volf = (f32 *)malloc(n * sizeof(f32));
  for (size_t i = 0; i < n; i++) volf[i] = (f32)img[i];
  f32 *sheet = (f32 *)malloc(n * sizeof(f32));
  if (!volf || !sheet) { fprintf(stderr, "oom\n"); return 1; }
  if (st_sheet_detect(volf, nz, ny, nx, NULL, sheet, NULL) != 0) {
    fprintf(stderr, "sheet detect failed\n"); return 1;
  }
  free(volf); free(img);

  // --- threshold -> binary surface mask -------------------------------------
  u8 *pred = (u8 *)malloc(n);
  size_t raw_on = 0;
  for (size_t i = 0; i < n; i++) { pred[i] = sheet[i] >= thresh; raw_on += pred[i]; }
  free(sheet);
  printf("sheetness >= %.3f : %zu voxels (%.1f%%)\n", thresh, raw_on, 100.0*raw_on/n);

  // --- classical postproc (the competition "win condition") -----------------
  size_t removed = remove_small_components(pred, nz, ny, nx, TOPO_CONN26, min_vox);
  majority_filter(pred, pred, nz, ny, nx, iters);   // iterated 3x3x3 median
  fill_holes(pred, nz, ny, nx);
  printf("postproc: dust<%zu removed %zu comps, median x%d, fill_holes\n",
         min_vox, removed, iters);

  // --- ignore-aware evaluation (label: 0 bg, 1 surface, 2 ignore) -----------
  surface_eval e = eval_surface(pred, lab, nz, ny, nx, tol, /*surface*/1, /*ignore*/2);
  printf("\n=== surface eval (tol=%d) ===\n", tol);
  printf("valid voxels : %zu  (%.1f%% of volume)\n", e.valid_vox, 100.0*e.valid_vox/n);
  printf("pred / gt vox: %zu / %zu\n", e.pred_vox, e.gt_vox);
  printf("Dice         : %.4f\n", e.dice);
  printf("Surface Dice : %.4f  (@tol %d)\n", e.surface_dice, tol);
  printf("pred topo    : b0=%ld b1=%ld b2=%ld   (gt b0=%ld)\n",
         e.pred_b0, e.pred_b1, e.pred_b2, e.gt_b0);

  free(pred); free(lab);
  return 0;
}
