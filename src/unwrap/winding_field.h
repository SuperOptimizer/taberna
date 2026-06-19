/* winding_field.h — Eulerian winding-number field over a volume.
 *
 * The dense backbone of the unwrap (docs/unwrapping-plan.md §2 stage 4/5): a
 * scalar field `winding` whose value is the (continuous) wrap index at each
 * foreground voxel. Level sets are sheets; the field value is the spiral
 * coordinate the Archimedean fit consumes.
 *
 * v0 implementation (to be improved): initialize from the umbilicus polar
 * coordinates (radius/dr + theta/2pi gives a continuous, monotone-in-radius
 * estimate), then a light masked Jacobi relaxation denoises it while honoring
 * absolute-winding seeds as Dirichlet conditions. The proper structure-tensor
 * orientation-integration Poisson solve replaces the relaxation later; the
 * interface stays the same.
 */
#ifndef TABERNA_WINDING_FIELD_H
#define TABERNA_WINDING_FIELD_H

#include "common/types.h"
#include "annotate/umbilicus.h"

typedef struct {
  f32 dr_per_winding;  // radial pitch (voxels per wrap); default 8
  int iters;           // Jacobi sweeps; default 50
  f32 omega;           // relaxation blend in [0,1]; default 0.3
  int warm_start;      // if nonzero, use the passed `winding` as the initial guess
                       // (coarse-to-fine multigrid) instead of the analytic init
  const f32 *forcing;  // optional (size nz*ny*nx), default NULL: the divergence of
                       // the target winding gradient (div of oriented sheet-normal /
                       // pitch). Turns the Laplacian smoother into a gradient-matching
                       // Poisson solve so level sets follow the ACTUAL sheets, not a
                       // radial assumption. See winding_normal_forcing().
  f32 anchor_lambda;   // default 0. When >0, softly pull the solution toward the
                       // initial field. Two uses: (1) with `forcing`, preserves the
                       // spiral monodromy a normal-divergence solve would collapse;
                       // (2) per-tile NATIVE re-solve warm-started from a global coarse
                       // winding — the anchor keeps each tile locked to the global
                       // coordinate (no seams) while local LOD0 detail refines it.
} wfield_params;

wfield_params winding_default_params(void);

// Contour-warp winding init: the real scroll is NOT a circular spiral. Warp the radial
// coordinate by the per-z outer-boundary shape R(theta) so iso-winding contours follow
// the actual (egg-shaped) scroll envelope: W = (r*Rref/R(theta))/pitch + theta/2pi.
// Robust (outer envelope only, no per-voxel normals). Writes `winding` (nz*ny*nx); feed
// it back as a warm_start init to winding_field_solve to refine. `pitch` = voxels/wrap.
void winding_contour_warp(const u8 *mask, int nz, int ny, int nx,
                          const umbilicus *umb, double pitch, f32 *winding);

// Solve the winding field over the foreground (mask != 0). `seed_value`/`seed_mask`
// (both size nz*ny*nx, may be NULL) pin absolute-winding Dirichlet seeds. Writes
// `winding` (nz*ny*nx). Returns 0 on success.
int winding_field_solve(const u8 *mask, int nz, int ny, int nx,
                        const umbilicus *umb, const wfield_params *p,
                        const f32 *seed_value, const u8 *seed_mask, f32 *winding);

#endif // TABERNA_WINDING_FIELD_H
