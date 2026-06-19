/* eval_seg — compare two label volumes (seg vs ground truth), both NRRD, and
 * print VOI split/merge, adapted-Rand, and (treating nonzero as foreground)
 * surface Dice.
 *
 *   eval_seg <seg.nrrd> <gt.nrrd> [--tol T]
 */
#include "io/nrrd.h"
#include "eval/metrics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s <seg.nrrd> <gt.nrrd> [--tol T]\n", argv[0]); return 2; }
  int tol = 2;
  for (int i = 3; i < argc; i++)
    if (!strcmp(argv[i], "--tol") && i + 1 < argc) tol = atoi(argv[++i]);

  int az, ay, ax, bz, by, bx;
  f32 *sf = nrrd_read_f32(argv[1], &az, &ay, &ax);
  f32 *gf = nrrd_read_f32(argv[2], &bz, &by, &bx);
  if (!sf || !gf) { fprintf(stderr, "read failed\n"); return 1; }
  if (az != bz || ay != by || ax != bx) { fprintf(stderr, "dim mismatch\n"); return 1; }
  size_t n = (size_t)az * ay * ax;

  u32 *seg = malloc(n * sizeof(u32));
  u32 *gt = malloc(n * sizeof(u32));
  u8 *sm = malloc(n), *gm = malloc(n);
  for (size_t i = 0; i < n; i++) {
    seg[i] = (u32)(sf[i] + 0.5f);
    gt[i] = (u32)(gf[i] + 0.5f);
    sm[i] = seg[i] != 0;
    gm[i] = gt[i] != 0;
  }

  eval_result r = eval_seg(seg, gt, n);
  double sd = surface_dice_at_tol(sm, gm, az, ay, ax, tol);
  printf("VOI: split=%.4f merge=%.4f total=%.4f bits   ARE=%.4f (P=%.4f R=%.4f)\n",
         r.voi_split, r.voi_merge, r.voi, r.are, r.precision, r.recall);
  printf("SurfaceDice@%d = %.4f\n", tol, sd);

  free(sf); free(gf); free(seg); free(gt); free(sm); free(gm);
  return 0;
}
