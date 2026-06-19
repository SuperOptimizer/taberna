/* spiral_fit.c — see spiral_fit.h. */
#include "unwrap/spiral_fit.h"

#include <math.h>
#include <stdlib.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

spiral_params spiral_fit_lsq(const double *theta, const double *radius, int n) {
  spiral_params sp = {0, 0, 0, n};
  if (n < 2) return sp;
  // ordinary least squares: radius ~ a + b*theta
  double st = 0, sr = 0, stt = 0, str = 0;
  for (int i = 0; i < n; i++) {
    st += theta[i]; sr += radius[i];
    stt += theta[i] * theta[i]; str += theta[i] * radius[i];
  }
  double dn = (double)n;
  double denom = dn * stt - st * st;
  if (fabs(denom) < 1e-12) return sp;
  sp.b = (dn * str - st * sr) / denom;
  sp.a = (sr - sp.b * st) / dn;

  double ss = 0;
  for (int i = 0; i < n; i++) {
    double e = radius[i] - (sp.a + sp.b * theta[i]);
    ss += e * e;
  }
  sp.rms = sqrt(ss / dn);
  return sp;
}

spiral_params spiral_fit_from_field(const u8 *mask, const f32 *winding,
                                    int nz, int ny, int nx,
                                    const umbilicus *umb, int stride) {
  size_t nynx = (size_t)ny * nx;
  const double TWO_PI = 6.283185307179586;
  if (stride < 1) stride = 1;

  // count then collect
  size_t cap = 0;
  for (size_t i = 0; i < (size_t)nz * nynx; i += stride) if (mask[i]) cap++;
  double *th = malloc((cap ? cap : 1) * sizeof(double));
  double *rd = malloc((cap ? cap : 1) * sizeof(double));
  int m = 0;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if ((i % (size_t)stride) != 0 || !mask[i]) continue;
        f32 theta, radius;
        umbilicus_polar(umb, (f32)z, (f32)y, (f32)x, &theta, &radius);
        th[m] = TWO_PI * (double)winding[i];
        rd[m] = radius;
        m++;
      }
  spiral_params sp = spiral_fit_lsq(th, rd, m);
  free(th); free(rd);
  return sp;
}
