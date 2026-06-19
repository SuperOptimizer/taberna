/* umbilicus.h — the scroll's center axis: a polyline of (z, y, x) control points,
 * the single highest-leverage human annotation (see docs/unwrapping-plan.md §3).
 * Provides the cylindrical coordinate system (theta, radius about the center at
 * each z) that the winding field and the Archimedean fit are built on.
 *
 * File format (our own, trivial): one control point per line, "z y x" (floats);
 * blank lines and lines beginning with '#' are ignored. Order does not matter;
 * points are sorted by z on load. Compatible to write from any viewer.
 */
#ifndef TABERNA_UMBILICUS_H
#define TABERNA_UMBILICUS_H

#include "common/types.h"

typedef struct {
  int  n;
  f32 *z, *y, *x;  // control points, sorted ascending by z
} umbilicus;

int  umbilicus_load(const char *path, umbilicus *u);
void umbilicus_free(umbilicus *u);
int  umbilicus_save(const char *path, const umbilicus *u);

// Estimate the scroll axis automatically from a material mask (no annotation).
// Splits z into `nctrl` bands; each control point is the band centroid refined to
// the center of radial symmetry (the point minimizing angular variance of the
// material's radial extent), then smoothed across z. Robust to off-center / tilted
// scrolls. Returns 0 and fills `*out` (caller frees with umbilicus_free), else -1.
int  umbilicus_estimate(const u8 *mask, int nz, int ny, int nx, int nctrl, umbilicus *out);

// Center (cy,cx) at height z (linear interpolation; clamped past the ends).
void umbilicus_center(const umbilicus *u, f32 z, f32 *cy, f32 *cx);

// Cylindrical coordinates of voxel (z,y,x): theta in (-pi,pi], radius >= 0.
void umbilicus_polar(const umbilicus *u, f32 z, f32 y, f32 x, f32 *theta, f32 *radius);

#endif // TABERNA_UMBILICUS_H
