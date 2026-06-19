/* sheet_tensor.c — see sheet_tensor.h. */
#include "segmentation/sheet_tensor.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

st_params st_default_params(void) {
  return (st_params){.sigma_grad = 1.0f, .sigma_tensor = 2.0f, .eps = 1.0f};
}

// ---- separable Gaussian blur (replicate borders) ----------------------------

// Build a normalized 1D Gaussian kernel; returns radius, writes kernel[0..2r].
static int gauss_kernel(double sigma, double **kernel_out) {
  int r = (int)ceil(3.0 * sigma);
  if (r < 1) r = 1;
  double *k = (double *)malloc((size_t)(2 * r + 1) * sizeof(double));
  double sum = 0.0;
  for (int i = -r; i <= r; i++) {
    double v = exp(-(double)i * i / (2.0 * sigma * sigma));
    k[i + r] = v;
    sum += v;
  }
  for (int i = 0; i < 2 * r + 1; i++) k[i] /= sum;
  *kernel_out = k;
  return r;
}

#define CLAMP(v, lo, hi) ((v) < (lo) ? (lo) : ((v) > (hi) ? (hi) : (v)))

static void blur_axis(const f32 *src, f32 *dst, int nz, int ny, int nx,
                      const double *k, int r, int axis) {
  size_t nynx = (size_t)ny * nx;
  #pragma omp parallel for schedule(static)
  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        double acc = 0.0;
        for (int t = -r; t <= r; t++) {
          int zz = z, yy = y, xx = x;
          if (axis == 0) xx = CLAMP(x + t, 0, nx - 1);
          else if (axis == 1) yy = CLAMP(y + t, 0, ny - 1);
          else zz = CLAMP(z + t, 0, nz - 1);
          acc += k[t + r] * (double)src[IDX(zz, yy, xx)];
        }
        dst[IDX(z, y, x)] = (f32)acc;
      }
    }
  }
}

// In-place separable Gaussian using one scratch buffer.
static void gaussian_blur(f32 *vol, int nz, int ny, int nx, double sigma, f32 *scratch) {
  if (sigma <= 0.0) return;
  double *k = NULL;
  int r = gauss_kernel(sigma, &k);
  blur_axis(vol, scratch, nz, ny, nx, k, r, 0);  // x
  blur_axis(scratch, vol, nz, ny, nx, k, r, 1);  // y
  blur_axis(vol, scratch, nz, ny, nx, k, r, 2);  // z
  memcpy(vol, scratch, (size_t)nz * ny * nx * sizeof(f32));
  free(k);
}

// ---- 3x3 symmetric eigensolver (cyclic Jacobi) ------------------------------

static void mat3_mul(const double a[3][3], const double b[3][3], double out[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += a[i][k] * b[k][j];
      out[i][j] = s;
    }
}

// Eigendecompose symmetric A. On return w[] holds eigenvalues sorted DESCENDING
// and V's columns the corresponding (orthonormal) eigenvectors.
static void jacobi3(double A[3][3], double w[3], double V[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) V[i][j] = (i == j) ? 1.0 : 0.0;

  for (int sweep = 0; sweep < 50; sweep++) {
    // pick the largest off-diagonal magnitude among (0,1),(0,2),(1,2)
    int p = 0, q = 1;
    double mx = fabs(A[0][1]);
    if (fabs(A[0][2]) > mx) { mx = fabs(A[0][2]); p = 0; q = 2; }
    if (fabs(A[1][2]) > mx) { mx = fabs(A[1][2]); p = 1; q = 2; }
    if (mx < 1e-20) break;

    double app = A[p][p], aqq = A[q][q], apq = A[p][q];
    double theta = (aqq - app) / (2.0 * apq);
    double t = (theta >= 0.0) ? 1.0 / (theta + sqrt(theta * theta + 1.0))
                              : -1.0 / (-theta + sqrt(theta * theta + 1.0));
    double c = 1.0 / sqrt(t * t + 1.0);
    double s = t * c;

    // Jacobi rotation J in the (p,q) plane; A <- J^T A J, V <- V J.
    double J[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    J[p][p] = c; J[q][q] = c; J[p][q] = s; J[q][p] = -s;
    double Jt[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Jt[i][j] = J[j][i];
    double tmp[3][3], An[3][3], Vn[3][3];
    mat3_mul(Jt, A, tmp);
    mat3_mul(tmp, J, An);
    mat3_mul(V, J, Vn);
    memcpy(A, An, sizeof An);
    memcpy(V, Vn, sizeof Vn);
  }

  for (int i = 0; i < 3; i++) w[i] = A[i][i];

  // sort eigenvalues (and eigenvector columns) descending
  for (int i = 0; i < 2; i++) {
    int mi = i;
    for (int j = i + 1; j < 3; j++)
      if (w[j] > w[mi]) mi = j;
    if (mi != i) {
      double tw = w[i]; w[i] = w[mi]; w[mi] = tw;
      for (int r = 0; r < 3; r++) {
        double tv = V[r][i]; V[r][i] = V[r][mi]; V[r][mi] = tv;
      }
    }
  }
}

// ---- sheet detection --------------------------------------------------------

int st_sheet_detect(const f32 *vol, int nz, int ny, int nx,
                    const st_params *p_in, f32 *sheetness, f32 *normal) {
  st_params p = p_in ? *p_in : st_default_params();
  size_t n = (size_t)nz * ny * nx;
  size_t nynx = (size_t)ny * nx;

  f32 *work = (f32 *)malloc(n * sizeof(f32));
  f32 *scratch = (f32 *)malloc(n * sizeof(f32));
  f32 *jxx = (f32 *)malloc(n * sizeof(f32));
  f32 *jyy = (f32 *)malloc(n * sizeof(f32));
  f32 *jzz = (f32 *)malloc(n * sizeof(f32));
  f32 *jxy = (f32 *)malloc(n * sizeof(f32));
  f32 *jxz = (f32 *)malloc(n * sizeof(f32));
  f32 *jyz = (f32 *)malloc(n * sizeof(f32));
  if (!work || !scratch || !jxx || !jyy || !jzz || !jxy || !jxz || !jyz) {
    free(work); free(scratch); free(jxx); free(jyy); free(jzz);
    free(jxy); free(jxz); free(jyz);
    return -1;
  }

  memcpy(work, vol, n * sizeof(f32));
  gaussian_blur(work, nz, ny, nx, p.sigma_grad, scratch);

  // gradients (central differences, replicate borders) -> tensor outer products
  #pragma omp parallel for schedule(static)
  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        int xm = x > 0 ? x - 1 : 0,      xp = x < nx - 1 ? x + 1 : nx - 1;
        int ym = y > 0 ? y - 1 : 0,      yp = y < ny - 1 ? y + 1 : ny - 1;
        int zm = z > 0 ? z - 1 : 0,      zp = z < nz - 1 ? z + 1 : nz - 1;
        double gx = 0.5 * ((double)work[IDX(z, y, xp)] - (double)work[IDX(z, y, xm)]);
        double gy = 0.5 * ((double)work[IDX(z, yp, x)] - (double)work[IDX(z, ym, x)]);
        double gz = 0.5 * ((double)work[IDX(zp, y, x)] - (double)work[IDX(zm, y, x)]);
        size_t i = IDX(z, y, x);
        jxx[i] = (f32)(gx * gx); jyy[i] = (f32)(gy * gy); jzz[i] = (f32)(gz * gz);
        jxy[i] = (f32)(gx * gy); jxz[i] = (f32)(gx * gz); jyz[i] = (f32)(gy * gz);
      }
    }
  }

  // smooth each tensor component at the integration scale
  gaussian_blur(jxx, nz, ny, nx, p.sigma_tensor, scratch);
  gaussian_blur(jyy, nz, ny, nx, p.sigma_tensor, scratch);
  gaussian_blur(jzz, nz, ny, nx, p.sigma_tensor, scratch);
  gaussian_blur(jxy, nz, ny, nx, p.sigma_tensor, scratch);
  gaussian_blur(jxz, nz, ny, nx, p.sigma_tensor, scratch);
  gaussian_blur(jyz, nz, ny, nx, p.sigma_tensor, scratch);

  #pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; i++) {
    double A[3][3] = {
        {jxx[i], jxy[i], jxz[i]},
        {jxy[i], jyy[i], jyz[i]},
        {jxz[i], jyz[i], jzz[i]},
    };
    double w[3], V[3][3];
    jacobi3(A, w, V);
    double l0 = w[0], l1 = w[1];  // l0 >= l1 >= l2
    if (sheetness) sheetness[i] = (f32)((l0 - l1) / (l0 + p.eps));
    if (normal) {
      // eigenvector of the largest eigenvalue = across-sheet normal
      double nxx = V[0][0], nyy = V[1][0], nzz = V[2][0];
      double len = sqrt(nxx * nxx + nyy * nyy + nzz * nzz);
      if (len < 1e-12) len = 1.0;
      normal[3 * i + 0] = (f32)(nxx / len);
      normal[3 * i + 1] = (f32)(nyy / len);
      normal[3 * i + 2] = (f32)(nzz / len);
    }
  }

  free(work); free(scratch); free(jxx); free(jyy); free(jzz);
  free(jxy); free(jxz); free(jyz);
  return 0;
}
