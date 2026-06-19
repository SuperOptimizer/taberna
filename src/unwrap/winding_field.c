/* winding_field.c — see winding_field.h. */
#include "unwrap/winding_field.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

wfield_params winding_default_params(void) {
  return (wfield_params){.dr_per_winding = 8.0f, .iters = 50, .omega = 0.3f, .warm_start = 0};
}

int winding_field_solve(const u8 *mask, int nz, int ny, int nx,
                        const umbilicus *umb, const wfield_params *p_in,
                        const f32 *seed_value, const u8 *seed_mask, f32 *winding) {
  wfield_params p = p_in ? *p_in : winding_default_params();
  size_t n = (size_t)nz * ny * nx, nynx = (size_t)ny * nx;
  const double TWO_PI = 6.283185307179586;

  // initialize: analytic polar estimate, OR (warm_start) keep the passed `winding`
  // as the initial guess for coarse-to-fine multigrid.
  #pragma omp parallel for schedule(static)
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (!mask[i]) { winding[i] = 0.0f; continue; }
        if (!p.warm_start) {
          f32 theta, radius;
          umbilicus_polar(umb, (f32)z, (f32)y, (f32)x, &theta, &radius);
          winding[i] = (f32)(radius / p.dr_per_winding + theta / TWO_PI);
        }
        if (seed_mask && seed_mask[i]) winding[i] = seed_value[i];
      }

  // Masked relaxation, IN PLACE via red-black Gauss-Seidel: no second n-sized
  // buffer (Jacobi needed an 8 GB ping-pong at LOD4), and GS converges ~2x faster
  // because it reuses neighbors already updated this sweep. Same-color voxels are
  // never 6-adjacent (a step flips (z+y+x) parity), so each color sweep is
  // parallel-safe. Honors Dirichlet seeds and the warm-start init unchanged.
  for (int it = 0; it < p.iters; it++) {
    for (int color = 0; color < 2; color++) {
      #pragma omp parallel for schedule(static)
      for (int z = 0; z < nz; z++)
        for (int y = 0; y < ny; y++)
          for (int x = 0; x < nx; x++) {
            if (((z + y + x) & 1) != color) continue;
            size_t i = IDX(z, y, x);
            if (!mask[i]) continue;                       // stays 0
            if (seed_mask && seed_mask[i]) continue;      // Dirichlet: pinned
            double acc = 0.0;
            int cnt = 0;
            int nb[6][3] = {{z, y, x - 1}, {z, y, x + 1}, {z, y - 1, x},
                            {z, y + 1, x}, {z - 1, y, x}, {z + 1, y, x}};
            for (int k = 0; k < 6; k++) {
              int zz = nb[k][0], yy = nb[k][1], xx = nb[k][2];
              if (xx < 0 || xx >= nx || yy < 0 || yy >= ny || zz < 0 || zz >= nz) continue;
              size_t j = IDX(zz, yy, xx);
              if (!mask[j]) continue;
              acc += winding[j];
              cnt++;
            }
            if (cnt) winding[i] = (f32)((1.0 - p.omega) * winding[i] + p.omega * (acc / cnt));
          }
    }
  }
  return 0;
}
