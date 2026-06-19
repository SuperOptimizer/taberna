/* wmetrics.h — winding-field correctness metrics (Tier-0, mostly GT-free), see
 * docs/unwrapping-plan.md §6.
 *   - winding_mvf: Monotonicity Violation Fraction. Cast rays from the umbilicus;
 *     the winding field should increase monotonically outward. Returns the
 *     fraction of radial steps where it does not (ideal 0). GT-free.
 *   - winding_satisfied_points: Paul-style "satisfied" check against absolute
 *     winding annotations — fraction of annotated points whose sampled field
 *     winding is within `tol` of their annotation. Needs annotations only.
 */
#ifndef TABERNA_WMETRICS_H
#define TABERNA_WMETRICS_H

#include "common/types.h"
#include "annotate/umbilicus.h"
#include "annotate/point_collection.h"

double winding_mvf(const u8 *mask, const f32 *winding, int nz, int ny, int nx,
                   const umbilicus *umb, int nrays, int max_radius);

// Returns the satisfied fraction; also writes counts via out params (may be NULL).
double winding_satisfied_points(const pc_set *s, const f32 *winding,
                                int nz, int ny, int nx, f32 tol,
                                int *n_satisfied, int *n_total);

#endif // TABERNA_WMETRICS_H
