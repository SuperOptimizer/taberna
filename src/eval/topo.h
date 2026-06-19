/* topo.h — computational-topology measures for binary surface masks.
 *
 * The Kaggle Surface-Detection metric is topology-dominated: a single tunnel or
 * cavity tanks the TopoScore (see docs/surface-detection-kaggle.md). The winning
 * solutions all leaned on cheap topology *detection* (Euler number) to find holes
 * before repairing them. This module provides the exact discrete invariants.
 *
 * Connectivity convention (the only self-consistent pairing used here):
 *   foreground = 26-connected   (matches the cubical complex: voxels sharing even
 *                                a corner are one cell-complex-connected blob)
 *   background =  6-connected   (its complement)
 * With this pairing the three Betti numbers and the Euler characteristic satisfy
 *   chi = b0 - b1 + b2     (b0 components, b1 tunnels/handles, b2 enclosed cavities)
 * and we compute chi exactly (cubical V-E+F-C), b0/b2 by labeling, b1 by the
 * identity. Verified on synthetic shapes in tools/test_topo.c (ball: chi 1; hollow
 * ball: chi 2, b2 1; solid torus: chi 0, b1 1).
 *
 * All volumes are z-major, x-fastest: idx = (z*ny + y)*nx + x. Masks are u8, any
 * nonzero = foreground.
 */
#ifndef TABERNA_TOPO_H
#define TABERNA_TOPO_H

#include <stddef.h>
#include "common/types.h"

typedef enum { TOPO_CONN6 = 6, TOPO_CONN26 = 26 } topo_conn;

/* Label connected components of the foreground (nonzero) voxels of `mask`.
 * Writes 1-based component ids into `labels` (0 = background) and returns the
 * component count. `labels` must hold nz*ny*nx u32s. */
u32 cc_label(const u8 *mask, int nz, int ny, int nx, topo_conn conn, u32 *labels);

/* Exact Euler characteristic of the foreground cubical complex (V - E + F - C). */
long euler_characteristic(const u8 *mask, int nz, int ny, int nx);

typedef struct {
  long b0;    // connected components (26-conn foreground)
  long b1;    // independent tunnels / handles (the Betti number the metric hates)
  long b2;    // enclosed cavities (6-conn background components not touching border)
  long chi;   // Euler characteristic = b0 - b1 + b2
} topo_betti;

topo_betti betti_numbers(const u8 *mask, int nz, int ny, int nx);

/* Betti numbers under the (6-conn foreground, 26-conn background) pairing — the
 * convention used by the Betti-Matching-3D cubical complex (corner-touching
 * voxels are NOT connected). b0 = 6-conn components; b2 = enclosed 26-conn
 * background cavities; chi via the voxel-as-vertex cubical complex (V−E+F−C with
 * 6-adjacency edges, unit-square 2-cells, unit-cube 3-cells); b1 = b0+b2−chi.
 * This is what the official TopoScore needs (see src/eval/score.c). */
topo_betti betti_numbers_6(const u8 *mask, int nz, int ny, int nx);

/* Tunnel screen: does the sub-box have nontrivial 1-cycles (b1 > 0)?  This is the
 * cheap per-window check the top solutions used to localize repairs. Returns b1 of
 * the cropped region. (z0..z0+dz) etc. clamped to the volume. */
long region_b1(const u8 *mask, int nz, int ny, int nx,
               int z0, int y0, int x0, int dz, int dy, int dx);

#endif // TABERNA_TOPO_H
