/* cubical.h — native cubical persistent homology for 3D scalar fields.
 *
 * taberna's own persistence engine. Built because clean topology of the
 * *continuous* sheetness/intensity field is core to getting perfect surfaces and
 * to unrolling (persistence-based simplification kills spurious tunnels/holes at
 * the source, not as binary post-hoc cleanup) — and because Betti matching is the
 * objective the surface metric rewards. See docs/surface-detection-kaggle.md.
 *
 * Construction (matches the Betti-Matching-3D oracle we validate against):
 *   - SUBLEVEL filtration (lower values enter first).
 *   - T-construction: voxels are the 0-cells (vertices); a k-cell enters at the
 *     MAX value of its constituent voxels. Foreground is therefore 6-connected
 *     (face-adjacency), background 26-connected — verified against the oracle.
 *   - Z/2 coefficients.
 *
 * v0 is correctness-first: a generic boundary-matrix reduction. Memory is O(cells)
 * ~ 4*nz*ny*nx, so this is for small / cropped / octant-sized volumes; the fast
 * union-find (dim 0) and cubical-specific dim 1/2 paths come next.
 *
 * Volumes are z-major, x-fastest: f[(z*ny + y)*nx + x].
 */
#ifndef TABERNA_CUBICAL_H
#define TABERNA_CUBICAL_H

#include <stddef.h>
#include "common/types.h"

#define TOPO_INF (3.4e38f)   /* death value for essential (never-dying) classes */

typedef struct {
  int dim;       /* homology dimension: 0, 1, 2 */
  f32 birth;     /* filtration value where the class is born */
  f32 death;     /* filtration value where it dies (TOPO_INF if essential) */
} pers_pair;

/* Compute the persistence diagram of `field` under the sublevel T-construction.
 * Returns a malloc'd array of pairs (caller frees) and writes the count to
 * *npairs. Essential classes have death == TOPO_INF. */
pers_pair *cubical_persistence(const f32 *field, int nz, int ny, int nx,
                               int *npairs);

#endif // TABERNA_CUBICAL_H
