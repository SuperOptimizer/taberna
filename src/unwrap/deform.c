/* deform.c — see deform.h. */
#include "unwrap/deform.h"

#include <math.h>
#include <stddef.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

int deform_build(const u8 *mask, const f32 *winding, int nz, int ny, int nx,
                 const umbilicus *umb, spiral_params sp, f32 *disp) {
  size_t nynx = (size_t)ny * nx;
  const double TWO_PI = 6.283185307179586;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        disp[3 * i + 0] = disp[3 * i + 1] = disp[3 * i + 2] = 0.0f;
        if (!mask[i]) continue;
        f32 cy, cx;
        umbilicus_center(umb, (f32)z, &cy, &cx);
        double dy = (double)y - cy, dx = (double)x - cx;
        double radius = sqrt(dy * dy + dx * dx);
        if (radius < 1e-6) continue;
        double theta_total = TWO_PI * (double)winding[i];
        double r_ideal = sp.a + sp.b * theta_total;
        double scale = r_ideal / radius;  // radial remap toward the ideal spiral
        double tx = cx + dx * scale, ty = cy + dy * scale;
        disp[3 * i + 0] = (f32)(tx - x);  // dx
        disp[3 * i + 1] = (f32)(ty - y);  // dy
        disp[3 * i + 2] = 0.0f;           // dz (radial remap only)
      }
  return 0;
}

double jacobian_fold_fraction(const f32 *disp, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  size_t folds = 0, total = 0;
  for (int z = 1; z < nz - 1; z++)
    for (int y = 1; y < ny - 1; y++)
      for (int x = 1; x < nx - 1; x++) {
        // central differences of each displacement component along each axis
        size_t xp = IDX(z, y, x + 1), xm = IDX(z, y, x - 1);
        size_t yp = IDX(z, y + 1, x), ym = IDX(z, y - 1, x);
        size_t zp = IDX(z + 1, y, x), zm = IDX(z - 1, y, x);
        double J[3][3];
        for (int c = 0; c < 3; c++) {
          J[c][0] = 0.5 * (disp[3 * xp + c] - disp[3 * xm + c]);  // d/dx
          J[c][1] = 0.5 * (disp[3 * yp + c] - disp[3 * ym + c]);  // d/dy
          J[c][2] = 0.5 * (disp[3 * zp + c] - disp[3 * zm + c]);  // d/dz
        }
        J[0][0] += 1.0; J[1][1] += 1.0; J[2][2] += 1.0;  // I + grad(disp)
        double det = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
                     J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
                     J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
        total++;
        if (det <= 0.0) folds++;
      }
  return total ? (double)folds / (double)total : 0.0;
}
