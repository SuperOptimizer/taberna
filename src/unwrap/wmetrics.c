/* wmetrics.c — see wmetrics.h. */
#include "unwrap/wmetrics.h"

#include <math.h>
#include <stddef.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

double winding_mvf(const u8 *mask, const f32 *winding, int nz, int ny, int nx,
                   const umbilicus *umb, int nrays, int max_radius) {
  size_t nynx = (size_t)ny * nx;
  const double TWO_PI = 6.283185307179586;
  int z = nz / 2;
  f32 cy, cx;
  umbilicus_center(umb, (f32)z, &cy, &cx);
  if (nrays < 1) nrays = 180;

  size_t steps = 0, violations = 0;
  for (int a = 0; a < nrays; a++) {
    double ang = TWO_PI * a / nrays;
    double cs = cos(ang), sn = sin(ang);
    int have_prev = 0;
    double prev = 0.0;
    for (int r = 1; r <= max_radius; r++) {
      int x = (int)lround(cx + r * cs);
      int y = (int)lround(cy + r * sn);
      if (x < 0 || x >= nx || y < 0 || y >= ny) break;
      size_t i = IDX(z, y, x);
      if (!mask[i]) { have_prev = 0; continue; }
      double w = winding[i];
      if (have_prev) {
        steps++;
        if (w < prev - 1e-4) violations++;
      }
      prev = w;
      have_prev = 1;
    }
  }
  return steps ? (double)violations / (double)steps : 0.0;
}

double winding_satisfied_points(const pc_set *s, const f32 *winding,
                                int nz, int ny, int nx, f32 tol,
                                int *n_satisfied, int *n_total) {
  size_t nynx = (size_t)ny * nx;
  int sat = 0, tot = 0;
  for (int c = 0; c < s->ncols; c++) {
    if (!s->cols[c].absolute) continue;  // v0: only absolute collections
    for (int k = 0; k < s->cols[c].npts; k++) {
      anno_point *pt = &s->cols[c].pts[k];
      int x = (int)lround(pt->x), y = (int)lround(pt->y), z = (int)lround(pt->z);
      if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) continue;
      tot++;
      double w = winding[IDX(z, y, x)];
      if (fabs(w - pt->wind) <= tol) sat++;
    }
  }
  if (n_satisfied) *n_satisfied = sat;
  if (n_total) *n_total = tot;
  return tot ? (double)sat / (double)tot : 0.0;
}
