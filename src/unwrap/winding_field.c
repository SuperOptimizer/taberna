/* winding_field.c — see winding_field.h. */
#include "unwrap/winding_field.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

void winding_contour_warp(const u8 *mask, int nz, int ny, int nx,
                          const umbilicus *umb, double pitch, f32 *winding) {
  const int NB = 180;
  size_t nynx = (size_t)ny * nx;
  const double TWO_PI = 6.283185307179586;
  #pragma omp parallel for schedule(static)
  for (int z = 0; z < nz; z++) {
    f32 cyf, cxf;
    umbilicus_center(umb, (f32)z, &cyf, &cxf);
    double Rb[180];
    for (int i = 0; i < NB; i++) Rb[i] = 0;
    // outer-boundary radius per angle bin = the scroll envelope at this z
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = (size_t)z * nynx + (size_t)y * nx + x;
        if (!mask[i]) continue;
        double dy = y - cyf, dx = x - cxf, r = sqrt(dy * dy + dx * dx);
        int b = (int)((atan2(dy, dx) + M_PI) / TWO_PI * NB);
        if (b < 0) b = 0; if (b >= NB) b = NB - 1;
        if (r > Rb[b]) Rb[b] = r;
      }
    // fill empty bins + circularly smooth the envelope
    for (int pass = 0; pass < 10; pass++) {
      double t[180];
      for (int i = 0; i < NB; i++) {
        double a = Rb[(i - 1 + NB) % NB], b = Rb[i], c = Rb[(i + 1) % NB];
        int cn = 0; double s = 0;
        if (a > 0) { s += a; cn++; } if (b > 0) { s += b; cn++; } if (c > 0) { s += c; cn++; }
        t[i] = cn ? s / cn : 0;
      }
      for (int i = 0; i < NB; i++) Rb[i] = t[i];
    }
    double Rref = 0; int cc = 0;
    for (int i = 0; i < NB; i++) if (Rb[i] > 0) { Rref += Rb[i]; cc++; }
    Rref = cc ? Rref / cc : 1;
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = (size_t)z * nynx + (size_t)y * nx + x;
        if (!mask[i]) { winding[i] = 0; continue; }
        double dy = y - cyf, dx = x - cxf, r = sqrt(dy * dy + dx * dx), th = atan2(dy, dx);
        int b = (int)((th + M_PI) / TWO_PI * NB);
        if (b < 0) b = 0; if (b >= NB) b = NB - 1;
        double Rt = Rb[b] > 1 ? Rb[b] : Rref;
        winding[i] = (f32)(r * (Rref / Rt) / pitch + th / TWO_PI);
      }
  }
}

wfield_params winding_default_params(void) {
  return (wfield_params){.dr_per_winding = 8.0f, .iters = 50, .omega = 0.3f,
                         .warm_start = 0, .forcing = NULL, .anchor_lambda = 0.0f};
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
  // Optional Tikhonov anchor: a copy of the initial field the solution is softly
  // pulled toward (preserves the angular monodromy the forcing would otherwise wash
  // out). Only allocated when actually used.
  f32 *anchor = NULL;
  double lam = (p.anchor_lambda > 0.0f) ? (double)p.anchor_lambda : 0.0;
  if (lam > 0.0) {
    anchor = (f32 *)malloc(n * sizeof(f32));
    if (!anchor) return -1;
    memcpy(anchor, winding, n * sizeof(f32));
  }

  // Early-exit when the field stops moving: a converged relaxation is identical
  // whether you stop at convergence or keep sweeping, so this only skips redundant
  // iterations (output-preserving). eps is tiny vs the field's range (~radius/pitch).
  const double eps = 1e-3;
  for (int it = 0; it < p.iters; it++) {
    double maxd = 0.0;
    for (int color = 0; color < 2; color++) {
      #pragma omp parallel for schedule(static) reduction(max:maxd)
      for (int z = 0; z < nz; z++)
        for (int y = 0; y < ny; y++)
          for (int x = 0; x < nx; x++) {
            // unit-stride x (keeps cache lines + vectorization); skip wrong parity
            if (((z + y + x) & 1) != color) continue;
            size_t i = IDX(z, y, x);
            if (!mask[i]) continue;                       // stays 0
            if (seed_mask && seed_mask[i]) continue;      // Dirichlet: pinned
            double acc = 0.0;
            double cnt = 0.0;
            if (p.tensor6 || p.normal) {
              // FULL anisotropic operator div(D grad W) = Dxx d2_xx + ... + 2Dxy d2_xy + ...
              // axial part keeps the (positive) W_p coefficient = 2*(Dxx+Dyy+Dzz); the off-diagonal
              // cross terms couple DIAGONAL neighbors (only added when all 4 in-plane diagonals are
              // material) and are lagged to the current W (stable Gauss-Seidel deferred correction).
              double Dzz, Dyy, Dxx, Dyz, Dxz, Dxy;
              if (p.tensor6) {
                const f32 *T = p.tensor6;
                Dzz = T[6*i+0]; Dyy = T[6*i+1]; Dxx = T[6*i+2];
                Dyz = T[6*i+3]; Dxz = T[6*i+4]; Dxy = T[6*i+5];
              } else {
                // build D = I-(1-a)*s*n n^T on the fly from the stored normal (no 6*N array)
                double nx = p.normal[3*i+0], ny = p.normal[3*i+1], nz = p.normal[3*i+2];
                double s = p.sheetness ? (double)p.sheetness[i] : 1.0; if (s<0) s=0; if (s>1) s=1;
                double k = (1.0 - (double)p.tensor_alpha) * s, n2 = nx*nx+ny*ny+nz*nz;
                if (n2 < 1e-6) { Dxx=Dyy=Dzz=1.0; Dxy=Dxz=Dyz=0.0; }
                else { Dxx=1.0-k*nx*nx; Dyy=1.0-k*ny*ny; Dzz=1.0-k*nz*nz;
                       Dxy=-k*nx*ny; Dxz=-k*nx*nz; Dyz=-k*ny*nz; }
              }
              #define MAT(zz,yy,xx) ((xx)>=0&&(xx)<nx&&(yy)>=0&&(yy)<ny&&(zz)>=0&&(zz)<nz&&mask[IDX(zz,yy,xx)])
              #define WV(zz,yy,xx) (winding[IDX(zz,yy,xx)])
              if (MAT(z,y,x-1)) { acc += Dxx*WV(z,y,x-1); cnt += Dxx; }
              if (MAT(z,y,x+1)) { acc += Dxx*WV(z,y,x+1); cnt += Dxx; }
              if (MAT(z,y-1,x)) { acc += Dyy*WV(z,y-1,x); cnt += Dyy; }
              if (MAT(z,y+1,x)) { acc += Dyy*WV(z,y+1,x); cnt += Dyy; }
              if (MAT(z-1,y,x)) { acc += Dzz*WV(z-1,y,x); cnt += Dzz; }
              if (MAT(z+1,y,x)) { acc += Dzz*WV(z+1,y,x); cnt += Dzz; }
              // cross terms: 2*Dab*d2_ab, d2_ab=(W[++]+W[--]-W[+-]-W[-+])/4; scaled by cross_relax
              // (lagged -> can exceed diagonal slack at strong anisotropy; <1 keeps GS stable).
              double cr = (p.cross_relax > 0.0f) ? 0.5*(double)p.cross_relax : 0.0;
              if (cr>0.0) {
              if (Dxy!=0 && MAT(z,y+1,x+1)&&MAT(z,y-1,x-1)&&MAT(z,y-1,x+1)&&MAT(z,y+1,x-1))
                acc += cr*Dxy*(WV(z,y+1,x+1)+WV(z,y-1,x-1)-WV(z,y-1,x+1)-WV(z,y+1,x-1));
              if (Dxz!=0 && MAT(z+1,y,x+1)&&MAT(z-1,y,x-1)&&MAT(z-1,y,x+1)&&MAT(z+1,y,x-1))
                acc += cr*Dxz*(WV(z+1,y,x+1)+WV(z-1,y,x-1)-WV(z-1,y,x+1)-WV(z+1,y,x-1));
              if (Dyz!=0 && MAT(z+1,y+1,x)&&MAT(z-1,y-1,x)&&MAT(z-1,y+1,x)&&MAT(z+1,y-1,x))
                acc += cr*Dyz*(WV(z+1,y+1,x)+WV(z-1,y-1,x)-WV(z-1,y+1,x)-WV(z+1,y-1,x));
              }
              #undef MAT
              #undef WV
            } else {
            int nb[6][3] = {{z, y, x - 1}, {z, y, x + 1}, {z, y - 1, x},
                            {z, y + 1, x}, {z - 1, y, x}, {z + 1, y, x}};
            for (int k = 0; k < 6; k++) {
              int zz = nb[k][0], yy = nb[k][1], xx = nb[k][2];
              if (xx < 0 || xx >= nx || yy < 0 || yy >= ny || zz < 0 || zz >= nz) continue;
              size_t j = IDX(zz, yy, xx);
              if (!mask[j]) continue;
              double wgt = 1.0;
              if (p.aniso) {
                // edge weight = harmonic-ish mean of the two voxels' axis diffusivity.
                // aniso layout per voxel: [wz, wy, wx]; map axis x->2, y->1, z->0.
                int ax = (k < 2) ? 2 : (k < 4) ? 1 : 0;
                double wi = p.aniso[3 * i + ax], wj = p.aniso[3 * j + ax];
                wgt = 0.5 * (wi + wj);
              }
              acc += wgt * winding[j];
              cnt += wgt;
            }
            }
            if (cnt > 0.0) {
              // Poisson gradient-matching with optional anchor and anisotropic weights:
              //   (sum_w + lam)*W[i] = sum_w*neighbors - div(g) + lam*anchor[i]
              // forcing[i]=div of desired gradient (NULL->0->plain Laplacian); the
              // lam term softly holds W near the analytic init (keeps the monodromy).
              double num = acc - (p.forcing ? (double)p.forcing[i] : 0.0);
              double den = cnt;
              if (anchor) { num += lam * anchor[i]; den += lam; }
              double target = num / den;
              f32 nv = (f32)((1.0 - p.omega) * winding[i] + p.omega * target);
              double d = fabs((double)nv - (double)winding[i]);
              if (d > maxd) maxd = d;
              winding[i] = nv;
            }
          }
    }
    if (maxd < eps) break;   // converged
  }
  free(anchor);
  return 0;
}
