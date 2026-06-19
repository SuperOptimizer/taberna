/* run_pipeline — end-to-end unwrap on one cube: sheet detection -> signed affinity
 * -> Mutex Watershed -> winding field -> Archimedean fit -> deformation field,
 * reporting the Tier-0 metrics along the way. The block-level driver the scale-up
 * (stitching across bricks) wraps later.
 *
 *   run_pipeline <vol.nrrd> <umbilicus.txt> [--out PREFIX] [--dr DR]
 */
#include "io/nrrd.h"
#include "annotate/umbilicus.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/affinity.h"
#include "segmentation/partition.h"
#include "unwrap/winding_field.h"
#include "unwrap/spiral_fit.h"
#include "unwrap/deform.h"
#include "unwrap/wmetrics.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s <vol.nrrd> <umbilicus.txt> [--out PREFIX] [--dr DR]\n", argv[0]);
    return 2;
  }
  const char *out = NULL;
  float dr = 8.0f;
  for (int i = 3; i < argc; i++) {
    if (!strcmp(argv[i], "--out") && i + 1 < argc) out = argv[++i];
    else if (!strcmp(argv[i], "--dr") && i + 1 < argc) dr = (float)atof(argv[++i]);
  }

  int nz, ny, nx;
  f32 *vol = nrrd_read_f32(argv[1], &nz, &ny, &nx);
  if (!vol) { fprintf(stderr, "read failed: %s\n", argv[1]); return 1; }
  size_t n = (size_t)nz * ny * nx;
  umbilicus umb;
  if (umbilicus_load(argv[2], &umb)) { fprintf(stderr, "umbilicus load failed\n"); return 1; }
  printf("volume %dx%dx%d, umbilicus %d control points\n", nz, ny, nx, umb.n);

  // 1) sheet detection
  f32 *sheet = malloc(n * sizeof(f32)), *norm = malloc(3 * n * sizeof(f32));
  st_params sp = st_default_params();
  st_sheet_detect(vol, nz, ny, nx, &sp, sheet, norm);

  // 2) foreground + signed affinity + 3) Mutex Watershed
  u8 *mask = malloc(n);
  for (size_t i = 0; i < n; i++) mask[i] = vol[i] > 128.0f;
  sgraph g; s32 *node_of; u32 *voxel_of;
  affinity_params ap = affinity_default_params();
  affinity_build_voxel(mask, sheet, norm, nz, ny, nx, &ap, &g, &node_of, &voxel_of);
  u32 *clab = malloc((size_t)(g.nnodes ? g.nnodes : 1) * sizeof(u32));
  int ncl = mws_partition(&g, clab);
  printf("segmentation: %d nodes, %d edges -> %d clusters\n", g.nnodes, g.nedges, ncl);

  // 4) winding field
  f32 *winding = malloc(n * sizeof(f32));
  wfield_params wp = winding_default_params();
  wp.dr_per_winding = dr;
  winding_field_solve(mask, nz, ny, nx, &umb, &wp, NULL, NULL, winding);
  double mvf = winding_mvf(mask, winding, nz, ny, nx, &umb, 180, (nx < ny ? nx : ny) / 2);
  printf("winding field: MVF (monotonicity violations) = %.4f  (want ~0)\n", mvf);

  // 5) Archimedean spiral fit
  spiral_params spr = spiral_fit_from_field(mask, winding, nz, ny, nx, &umb, 8);
  printf("spiral fit: r = %.3f + %.5f*theta   rms=%.3f voxels  (n=%d samples)\n",
         spr.a, spr.b, spr.rms, spr.nsamples);

  // 6) deformation field + invertibility guard
  f32 *disp = malloc(3 * n * sizeof(f32));
  deform_build(mask, winding, nz, ny, nx, &umb, spr, disp);
  double folds = jacobian_fold_fraction(disp, nz, ny, nx);
  printf("deformation: det-J fold fraction = %.4f  (want 0 = invertible)\n", folds);

  if (out) {
    char path[512];
    snprintf(path, sizeof path, "%s_sheetness.nrrd", out); nrrd_write_f32(path, sheet, nz, ny, nx);
    snprintf(path, sizeof path, "%s_winding.nrrd", out);   nrrd_write_f32(path, winding, nz, ny, nx);
    u32 *seg = calloc(n, sizeof(u32));
    u8 *segu8 = malloc(n);
    for (int v = 0; v < g.nnodes; v++) seg[voxel_of[v]] = clab[v] + 1;
    for (size_t i = 0; i < n; i++) segu8[i] = (u8)(seg[i] % 256);
    snprintf(path, sizeof path, "%s_seg.nrrd", out); nrrd_write_u8(path, segu8, nz, ny, nx);
    free(seg); free(segu8);
    printf("wrote %s_{sheetness,winding,seg}.nrrd\n", out);
  }

  sgraph_free(&g); umbilicus_free(&umb);
  free(vol); free(sheet); free(norm); free(mask);
  free(node_of); free(voxel_of); free(clab); free(winding); free(disp);
  return 0;
}
