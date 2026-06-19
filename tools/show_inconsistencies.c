/* show_inconsistencies — the "where to annotate next" report. Given a winding
 * field (NRRD) and a point-collection file, sample the field at each absolute
 * annotation and flag collections whose sampled winding disagrees with the
 * annotation (residual > tol) — these are the spots the human should correct
 * (docs/unwrapping-plan.md §3, the targeted loop). A first proxy for the
 * loop-closure holonomy check.
 *
 *   show_inconsistencies <winding.nrrd> <points.txt> [--tol T]
 */
#include "io/nrrd.h"
#include "annotate/point_collection.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: %s <winding.nrrd> <points.txt> [--tol T]\n", argv[0]); return 2; }
  double tol = 0.5;
  for (int i = 3; i < argc; i++)
    if (!strcmp(argv[i], "--tol") && i + 1 < argc) tol = atof(argv[++i]);

  int nz, ny, nx;
  f32 *w = nrrd_read_f32(argv[1], &nz, &ny, &nx);
  if (!w) { fprintf(stderr, "read failed: %s\n", argv[1]); return 1; }
  size_t nynx = (size_t)ny * nx;

  pc_set s;
  if (pc_load(argv[2], &s)) { fprintf(stderr, "points load failed: %s\n", argv[2]); return 1; }

  int flagged = 0;
  for (int c = 0; c < s.ncols; c++) {
    point_collection *col = &s.cols[c];
    double max_res = 0.0;
    int worst = -1;
    for (int k = 0; k < col->npts; k++) {
      anno_point *p = &col->pts[k];
      int x = (int)lround(p->x), y = (int)lround(p->y), z = (int)lround(p->z);
      if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) continue;
      double sampled = w[(size_t)z * nynx + (size_t)y * nx + x];
      double res = col->absolute ? fabs(sampled - p->wind) : 0.0;
      if (res > max_res) { max_res = res; worst = k; }
    }
    if (col->absolute && max_res > tol) {
      flagged++;
      printf("INCONSISTENT  collection '%s'  max winding residual=%.3f at point %d\n",
             col->id, max_res, worst);
    }
  }
  printf("%d/%d collections flagged for correction (tol=%.2f)\n", flagged, s.ncols, tol);

  pc_free(&s);
  free(w);
  return 0;
}
