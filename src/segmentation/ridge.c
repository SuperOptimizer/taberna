/* ridge.c — see ridge.h. */
#include "segmentation/ridge.h"

#include <math.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

/* Trilinear sample of `v` at continuous (zf,yf,xf), clamped to the volume. */
static f32 trilin(const f32 *v, int nz, int ny, int nx, f32 zf, f32 yf, f32 xf) {
  size_t nynx = (size_t)ny * nx;
  if (zf < 0) zf = 0; if (zf > nz - 1) zf = (f32)(nz - 1);
  if (yf < 0) yf = 0; if (yf > ny - 1) yf = (f32)(ny - 1);
  if (xf < 0) xf = 0; if (xf > nx - 1) xf = (f32)(nx - 1);
  int z0 = (int)zf, y0 = (int)yf, x0 = (int)xf;
  int z1 = z0 < nz - 1 ? z0 + 1 : z0;
  int y1 = y0 < ny - 1 ? y0 + 1 : y0;
  int x1 = x0 < nx - 1 ? x0 + 1 : x0;
  f32 dz = zf - z0, dy = yf - y0, dx = xf - x0;
  f32 c00 = v[IDX(z0,y0,x0)]*(1-dx) + v[IDX(z0,y0,x1)]*dx;
  f32 c01 = v[IDX(z0,y1,x0)]*(1-dx) + v[IDX(z0,y1,x1)]*dx;
  f32 c10 = v[IDX(z1,y0,x0)]*(1-dx) + v[IDX(z1,y0,x1)]*dx;
  f32 c11 = v[IDX(z1,y1,x0)]*(1-dx) + v[IDX(z1,y1,x1)]*dx;
  f32 c0 = c00*(1-dy) + c01*dy;
  f32 c1 = c10*(1-dy) + c11*dy;
  return c0*(1-dz) + c1*dz;
}

size_t ridge_nms(const f32 *ct, const f32 *sheet, const f32 *normal,
                 int nz, int ny, int nx,
                 f32 s_min, f32 i_min, f32 step, u8 *out) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  size_t kept = 0;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (sheet[i] < s_min || ct[i] < i_min) { out[i] = 0; continue; }
        f32 nxc = normal[3*i+0], nyc = normal[3*i+1], nzc = normal[3*i+2];
        f32 c = ct[i];
        f32 fwd = trilin(ct, nz, ny, nx, z + step*nzc, y + step*nyc, x + step*nxc);
        f32 bwd = trilin(ct, nz, ny, nx, z - step*nzc, y - step*nyc, x - step*nxc);
        if (c < fwd || c < bwd) { out[i] = 0; continue; }   // not a local max
        // Sub-voxel peak of the parabola through (bwd, c, fwd) at offsets (-1,0,1)·step.
        // denom < 0 => concave (a real ridge); peak offset t* in units of step.
        f32 denom = bwd - 2.0f*c + fwd;
        if (denom < -1e-3f) {
          f32 tstar = 0.5f * (bwd - fwd) / denom;
          if (fabsf(tstar) > 0.5f) { out[i] = 0; continue; }  // peak sits in a neighbor
        }
        // (near-flat plateau max: denom ~ 0 -> keep this voxel)
        out[i] = 1; kept++;
      }
  (void)n;
  return kept;
}
