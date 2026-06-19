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
} wfield_params;

wfield_params winding_default_params(void);

// Solve the winding field over the foreground (mask != 0). `seed_value`/`seed_mask`
// (both size nz*ny*nx, may be NULL) pin absolute-winding Dirichlet seeds. Writes
// `winding` (nz*ny*nx). Returns 0 on success.
int winding_field_solve(const u8 *mask, int nz, int ny, int nx,
                        const umbilicus *umb, const wfield_params *p,
                        const f32 *seed_value, const u8 *seed_mask, f32 *winding);

#endif // TABERNA_WINDING_FIELD_H
