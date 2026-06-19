/* deform.h — build the recorded deformation field that maps the scanned scroll to
 * its ideal Archimedean spiral, and the invertibility guard (det-J > 0).
 *
 * deform_build computes, per foreground voxel, a displacement that moves it
 * radially so its radius matches the ideal spiral radius for its winding value
 * (r_ideal = a + b * 2*pi*winding), keeping theta and z. The displacement field
 * is interleaved (dx,dy,dz) per voxel — the same convention as the normal field.
 * jacobian_fold_fraction reports the fraction of voxels where det(I + grad(disp))
 * <= 0 (folds / non-invertibility) — the strongest Tier-0 wrap-merge guard.
 */
#ifndef TABERNA_DEFORM_H
#define TABERNA_DEFORM_H

#include "common/types.h"
#include "annotate/umbilicus.h"
#include "unwrap/spiral_fit.h"

// Build displacement (3*nz*ny*nx, interleaved dx,dy,dz). Returns 0 on success.
int deform_build(const u8 *mask, const f32 *winding, int nz, int ny, int nx,
                 const umbilicus *umb, spiral_params sp, f32 *disp);

// Fraction of interior voxels with det(Jacobian) <= 0 (folds). disp interleaved.
double jacobian_fold_fraction(const f32 *disp, int nz, int ny, int nx);

#endif // TABERNA_DEFORM_H
