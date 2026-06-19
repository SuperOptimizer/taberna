/* point_collection.h — winding-annotated point collections + link constraints,
 * the targeted human corrections that drive the fit (docs/unwrapping-plan.md §3).
 *
 * A collection groups points that should agree on absolute winding (up to their
 * relative offsets). An ABSOLUTE collection's `wind` values are absolute winding
 * numbers (calibrate the field's integration constant); a RELATIVE collection's
 * `wind` values are only meaningful as differences (link wraps across breaks).
 * Link constraints are must-link / cannot-link pairs given as two 3D positions,
 * consumed by the signed-graph partition (cannot-link => hard mutex).
 *
 * File format (our own, trivial line-based):
 *   collection <id> absolute|relative
 *   point <z> <y> <x> <wind>
 *   ...
 *   cannotlink <z y x> <z y x>
 *   mustlink   <z y x> <z y x>
 *   # comments and blank lines ignored
 */
#ifndef TABERNA_POINT_COLLECTION_H
#define TABERNA_POINT_COLLECTION_H

#include "common/types.h"

typedef struct { f32 z, y, x, wind; } anno_point;

typedef struct {
  char        id[64];
  int         absolute;  // 1 = absolute winding, 0 = relative
  anno_point *pts;
  int         npts, cap;
} point_collection;

typedef struct {
  f32 az, ay, ax, bz, by, bx;
  int cannot;            // 1 = cannot-link, 0 = must-link
} link_constraint;

typedef struct {
  point_collection *cols;
  int               ncols, cols_cap;
  link_constraint  *links;
  int               nlinks, links_cap;
} pc_set;

int  pc_load(const char *path, pc_set *s);
void pc_free(pc_set *s);

#endif // TABERNA_POINT_COLLECTION_H
