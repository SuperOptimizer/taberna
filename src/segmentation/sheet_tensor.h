/* sheet_tensor.h — classical (no-ML) structure-tensor sheet detection.
 *
 * The first stage of taberna's unwrapping pipeline (see docs/unwrapping-plan.md).
 * For each voxel it produces two things the rest of the pipeline needs:
 *   - a SHEETNESS scalar in [0,1] — how plane-like the local neighborhood is, and
 *   - the SHEET NORMAL — the unit direction across the sheet.
 *
 * Method: optionally pre-smooth (scale sigma_grad), compute the gradient g, form
 * the structure tensor T = smooth(g g^T) at integration scale sigma_tensor (= the
 * scale ρ; keep it BELOW the inter-wrap gap, ~1-2 voxels, or it averages the
 * normals of two adjacent wraps together), then eigendecompose the 3x3 symmetric T
 * per voxel: eigenvalues l0 >= l1 >= l2. A planar sheet has one large eigenvalue
 * (across the sheet) and two small (in-plane), so:
 *     sheetness = (l0 - l1) / (l0 + eps)        in [0,1]
 *     normal    = eigenvector of l0             (the across-sheet direction)
 *
 * All volumes are (nz,ny,nx), z-major, x-fastest: voxel (z,y,x) at v[(z*ny+y)*nx+x].
 * NOTE: this runs on whatever scalar field the caller supplies — in taberna that is
 * the matter-compressor-decompressed u8 CT (~1% MAE), so the detector must tolerate
 * a few graylevels of compression noise (see experiment E2b in the plan).
 */
#ifndef TABERNA_SHEET_TENSOR_H
#define TABERNA_SHEET_TENSOR_H

#include "common/types.h"

typedef struct {
  f32 sigma_grad;    // pre-smoothing scale; <= 0 disables (typical 0.5-1.0)
  f32 sigma_tensor;  // integration scale ρ for the tensor (typical 1.0-2.0)
  f32 eps;           // sheetness denominator floor (typical 1e-3 * value range)
} st_params;

// Sensible defaults (sigma_grad=1.0, sigma_tensor=2.0, eps=1.0).
st_params st_default_params(void);

// Run sheet detection. `vol` is the (nz,ny,nx) input scalar field.
//   sheetness : OUT nz*ny*nx, may be NULL if not wanted.
//   normal    : OUT 3*nz*ny*nx, interleaved (nx,ny,nz) per voxel (unit length;
//               sign is arbitrary — it's an axis, not an oriented direction);
//               may be NULL.
// `p` NULL -> defaults. Returns 0 on success, non-zero on allocation failure.
int st_sheet_detect(const f32 *vol, int nz, int ny, int nx,
                    const st_params *p, f32 *sheetness, f32 *normal);

#endif // TABERNA_SHEET_TENSOR_H
