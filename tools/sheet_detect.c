/* sheet_detect — run structure-tensor sheet detection on a NRRD volume and write
 * the sheetness field (and optionally the normals) as NRRD, for inspection / as
 * the input to the signed-affinity experiment (E1).
 *
 * Usage:
 *   sheet_detect <in.nrrd> <out_sheetness.nrrd> [out_normals.nrrd]
 *                [--sigma-grad S] [--sigma-tensor R] [--eps E]
 */
#include "io/nrrd.h"
#include "segmentation/sheet_tensor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr,
            "usage: %s <in.nrrd> <out_sheetness.nrrd> [out_normals.nrrd]\n"
            "          [--sigma-grad S] [--sigma-tensor R] [--eps E]\n",
            argv[0]);
    return 2;
  }
  const char *in = argv[1];
  const char *out_sheet = argv[2];
  const char *out_norm = NULL;

  st_params p = st_default_params();
  for (int i = 3; i < argc; i++) {
    if (!strcmp(argv[i], "--sigma-grad") && i + 1 < argc) p.sigma_grad = (f32)atof(argv[++i]);
    else if (!strcmp(argv[i], "--sigma-tensor") && i + 1 < argc) p.sigma_tensor = (f32)atof(argv[++i]);
    else if (!strcmp(argv[i], "--eps") && i + 1 < argc) p.eps = (f32)atof(argv[++i]);
    else if (argv[i][0] != '-') out_norm = argv[i];
    else { fprintf(stderr, "unknown arg: %s\n", argv[i]); return 2; }
  }

  int nz, ny, nx;
  f32 *vol = nrrd_read_f32(in, &nz, &ny, &nx);
  if (!vol) { fprintf(stderr, "failed to read %s\n", in); return 1; }
  size_t n = (size_t)nz * ny * nx;
  printf("loaded %s : %d x %d x %d (z,y,x) = %zu voxels\n", in, nz, ny, nx, n);
  printf("params: sigma_grad=%.2f sigma_tensor=%.2f eps=%.3f\n",
         p.sigma_grad, p.sigma_tensor, p.eps);

  f32 *sheet = (f32 *)malloc(n * sizeof(f32));
  f32 *norm = out_norm ? (f32 *)malloc(3 * n * sizeof(f32)) : NULL;
  if (!sheet || (out_norm && !norm)) { fprintf(stderr, "oom\n"); return 1; }

  if (st_sheet_detect(vol, nz, ny, nx, &p, sheet, norm)) {
    fprintf(stderr, "sheet detect failed (oom)\n");
    return 1;
  }

  // stats
  double mn = 1e30, mx = -1e30, sum = 0.0;
  size_t hi = 0;
  for (size_t i = 0; i < n; i++) {
    double v = sheet[i];
    if (v < mn) mn = v;
    if (v > mx) mx = v;
    sum += v;
    if (v >= 0.5) hi++;
  }
  printf("sheetness: min=%.3f max=%.3f mean=%.3f  frac>=0.5 = %.3f%%\n",
         mn, mx, sum / (double)n, 100.0 * (double)hi / (double)n);

  if (nrrd_write_f32(out_sheet, sheet, nz, ny, nx))
    fprintf(stderr, "warning: failed to write %s\n", out_sheet);
  else
    printf("wrote %s\n", out_sheet);

  if (out_norm) {
    // store the 3-vector field as a 4D-ish raw: write as f32 with x size 3*nx is
    // misleading; instead write a (3*nz, ny, nx) stack is also confusing. For now
    // emit the per-voxel |normal| sanity field is unnecessary — write interleaved
    // as a flat (nz, ny, 3*nx) raw so it's loadable; documented for the caller.
    if (nrrd_write_f32(out_norm, norm, nz, ny, 3 * nx))
      fprintf(stderr, "warning: failed to write %s\n", out_norm);
    else
      printf("wrote %s (normals interleaved nx,ny,nz per voxel; axis0 = 3*nx)\n", out_norm);
  }

  free(vol); free(sheet); free(norm);
  return 0;
}
