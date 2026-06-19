/* ced.h — coherence-enhancing (sheet) diffusion.
 *
 * Our ridge detector is porous because the sheetness/intensity field has gaps
 * along the sheets; thickened, those gaps become tunnels (b1) that the metric
 * punishes hard. CED fixes them at the source: anisotropic diffusion that smooths
 * strongly ALONG the local sheet (the tangent plane, perpendicular to the
 * structure-tensor normal n) and weakly ACROSS it, so weak/broken sheet signal
 * gets connected without blurring adjacent wraps together. (This is the 7th-place
 * classical front-end, in a normal-only form.)
 *
 * Diffusion tensor from the unit normal n:  D = c_perp*I - (c_perp-c_norm)*n n^T,
 * i.e. diffusivity c_perp in the sheet plane, c_norm across. We use c_perp tied to
 * sheetness (diffuse more where the sheet is confident). Evolved by an explicit
 * div(D grad u) scheme. Volumes z-major, x-fastest.
 */
#ifndef TABERNA_CED_H
#define TABERNA_CED_H

#include "common/types.h"
#include "segmentation/sheet_tensor.h"

typedef struct {
  st_params st;     // structure-tensor scales for the normal field
  f32 c_norm;       // across-sheet diffusivity (small, e.g. 0.05)
  f32 dt;           // explicit time step (e.g. 0.12; <= 1/6 for 3D stability)
  int iters;        // number of diffusion steps (e.g. 12)
} ced_params;

ced_params ced_default_params(void);

/* Diffuse `vol` in place (nz,ny,nx). The normal field is computed once from the
 * input. Returns 0 on success. */
int ced_diffuse(f32 *vol, int nz, int ny, int nx, const ced_params *p);

#endif // TABERNA_CED_H
