/* spiral_fit.h — fit an ideal Archimedean spiral r = a + b*theta_total to the
 * recovered winding field (docs/unwrapping-plan.md §2 stage 6). theta_total is the
 * accumulated angle = 2*pi*winding, so b is the radial pitch per radian
 * (= dr_per_winding / 2pi). Closed-form least squares.
 */
#ifndef TABERNA_SPIRAL_FIT_H
#define TABERNA_SPIRAL_FIT_H

#include "common/types.h"
#include "annotate/umbilicus.h"

typedef struct {
  double a, b;       // r = a + b * theta_total
  double rms;        // RMS residual of the fit (voxels)
  int    nsamples;
} spiral_params;

// Least-squares fit from paired samples (theta_total[i], radius[i]).
spiral_params spiral_fit_lsq(const double *theta_total, const double *radius, int n);

// Convenience: collect (theta_total, radius) samples over the foreground from a
// winding field + umbilicus and fit. `stride` subsamples voxels (>=1).
spiral_params spiral_fit_from_field(const u8 *mask, const f32 *winding,
                                    int nz, int ny, int nx,
                                    const umbilicus *umb, int stride);

#endif // TABERNA_SPIRAL_FIT_H
