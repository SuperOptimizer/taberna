/* hessian.h — multi-scale Hessian (Frangi) sheet/plate detection.
 *
 * Alternative front end to the structure tensor. The Hessian's eigenstructure is
 * centered on the sheet RIDGE (second-derivative peak = the bright papyrus
 * centerline), so it places the surface better than the gradient-based structure
 * tensor, and a multi-scale max over sigma catches sheets of varying thickness
 * (less porosity from scale mismatch).
 *
 * Bright plate (sheet): eigenvalues |l1|<=|l2|<=|l3| have l3 large and NEGATIVE
 * (curvature across the bright sheet) and l1,l2 ~ 0. Frangi plate measure:
 *   sheetness = [l3<0] * exp(-Ra^2/2a^2) * (1 - exp(-S^2/2c^2)),  Ra=|l2|/|l3|,
 *   S = ||lambda|| (structure strength). normal = eigenvector of l3.
 * Multi-scale: take the max sheetness over the scale set, keeping that scale's
 * normal. Volumes z-major, x-fastest.
 */
#ifndef TABERNA_HESSIAN_H
#define TABERNA_HESSIAN_H

#include "common/types.h"

typedef struct {
  const f32 *sigmas;   // scale set (Gaussian sigma per scale)
  int nsigma;
  f32 alpha;           // plate-vs-line sensitivity (Ra term), e.g. 0.5
  f32 beta_c;          // structure-strength cutoff c; <=0 => auto (half max ||H||)
} hess_params;

/* Multi-scale Frangi plate sheetness in [0,1] (`sheetness`, nz*ny*nx; may be NULL)
 * and across-sheet normal (`normal`, 3*nz*ny*nx interleaved x,y,z, unit; may be
 * NULL) at the best-responding scale per voxel. Returns 0 on success. */
int hessian_sheet_detect(const f32 *vol, int nz, int ny, int nx,
                         const hess_params *p, f32 *sheetness, f32 *normal);

#endif // TABERNA_HESSIAN_H
