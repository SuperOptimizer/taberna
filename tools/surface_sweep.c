/* surface_sweep.c — tune the classical ridge detector on one labeled cube.
 *
 * The structure tensor (sheetness + normal) is the expensive step, so compute it
 * ONCE, then sweep the cheap ridge-NMS + light postproc + scoring over a grid of
 * (sheetness floor, intensity floor). Unlike the thick-mask competition pipeline,
 * the ridge output is ~1-2 voxels thick (matching the label), so we do NOT run the
 * eroding median filter here — just dust removal + hole fill.
 *
 * Usage: surface_sweep IMAGE.tif LABEL.tif [step=1.0] [min_voxels=100] [tol=2]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/metrics.h"

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s IMAGE.tif LABEL.tif [step=1.0] [min_voxels=100] [tol=2]\n", argv[0]);
    return 2;
  }
  float step    = argc > 3 ? (float)atof(argv[3]) : 1.0f;
  size_t min_vox= argc > 4 ? (size_t)atoll(argv[4]) : 100;
  int   tol     = argc > 5 ? atoi(argv[5]) : 2;

  int nz, ny, nx, lz, ly, lx;
  u8 *img = tiff_load_u8(argv[1], &nz, &ny, &nx);
  u8 *lab = tiff_load_u8(argv[2], &lz, &ly, &lx);
  if (!img || !lab || lz != nz || ly != ny || lx != nx) {
    fprintf(stderr, "load/dim error\n"); return 1;
  }
  size_t n = (size_t)nz * ny * nx;

  // gt surface fraction (the density the detector should aim for)
  size_t gt = 0, valid = 0;
  for (size_t i = 0; i < n; i++) { gt += lab[i]==1; valid += lab[i]!=2; }
  printf("%dx%dx%d  gt surface=%zu (%.1f%% of valid %zu)\n",
         nz,ny,nx, gt, 100.0*gt/valid, valid);

  // --- structure tensor ONCE ------------------------------------------------
  f32 *volf = (f32*)malloc(n*sizeof(f32));
  for (size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet = (f32*)malloc(n*sizeof(f32));
  f32 *normal= (f32*)malloc(3*n*sizeof(f32));
  if (!volf||!sheet||!normal){fprintf(stderr,"oom\n");return 1;}
  fprintf(stderr,"computing structure tensor...\n");
  if (st_sheet_detect(volf, nz,ny,nx, NULL, sheet, normal)!=0){fprintf(stderr,"ST failed\n");return 1;}

  // --- sweep ----------------------------------------------------------------
  float s_grid[] = {0.05f, 0.15f, 0.30f, 0.50f};
  float i_grid[] = {0.f, 80.f, 120.f, 160.f};
  u8 *pred = (u8*)malloc(n);
  printf("\n s_min  i_min   kept    predV    Dice   SurfDice  b0    b1\n");
  printf("-----------------------------------------------------------------\n");
  for (size_t a=0;a<sizeof(s_grid)/sizeof(*s_grid);a++)
    for (size_t b=0;b<sizeof(i_grid)/sizeof(*i_grid);b++) {
      float s_min=s_grid[a], i_min=i_grid[b];
      size_t kept = ridge_nms(volf, sheet, normal, nz,ny,nx, s_min, i_min, step, pred);
      remove_small_components(pred, nz,ny,nx, TOPO_CONN26, min_vox);
      fill_holes(pred, nz,ny,nx);
      surface_eval e = eval_surface(pred, lab, nz,ny,nx, tol, 1, 2);
      printf(" %.2f   %5.0f  %7zu  %7zu  %.4f  %.4f  %4ld  %5ld\n",
             s_min, i_min, kept, e.pred_vox, e.dice, e.surface_dice, e.pred_b0, e.pred_b1);
    }

  free(volf);free(sheet);free(normal);free(pred);free(img);free(lab);
  return 0;
}
