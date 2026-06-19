/* affinity.h — build a SIGNED region-adjacency graph from sheet detection, for
 * signed-graph partitioning into per-wrap segments (the core of experiment E1).
 *
 * Signed edge convention: w > 0 attractive ("same wrap"), w < 0 repulsive
 * ("different wraps, must not merge"). Between two adjacent sheet voxels we form
 *     across = |n · d|   (d = the +x/+y/+z edge direction, n = sheet normal)
 *     w = k_attract · (1 - across)  -  k_repel · (across · sheetness_boundary)
 * i.e. an in-plane contact (edge perpendicular to the normal) attracts, while a
 * contact running ALONG the normal through high-sheetness material repels (the
 * signature of crossing from one wrap into the touching next one).
 *
 * This first cut works at VOXEL resolution over the foreground (mask != 0) — the
 * separation between non-touching wraps comes for free from the gap voxels being
 * outside the mask; the repulsion term is what must catch TOUCHING wraps, and
 * tuning/strengthening it (OOF, phase-symmetry, human cannot-links) is exactly
 * what E1 evaluates. Supervoxel coarsening (snic.h) is the path to TB scale and
 * slots in later behind the same sgraph.
 *
 * Volumes are (nz,ny,nx), z-major x-fastest. normal is 3*n interleaved (nx,ny,nz).
 */
#ifndef TABERNA_AFFINITY_H
#define TABERNA_AFFINITY_H

#include "common/types.h"

typedef struct { u32 a, b; f32 w; } sg_edge;  // signed weight

typedef struct {
  sg_edge *edges;
  int      nedges;
  int      nnodes;
} sgraph;

typedef struct {
  f32 k_attract;  // in-plane attraction gain (default 1.0)
  f32 k_repel;    // across-sheet repulsion gain (default 1.0)
} affinity_params;

affinity_params affinity_default_params(void);

// Build a voxel-level signed RAG over the foreground (mask != 0). Allocates the
// edge list inside `g`, and (if non-NULL) a `node_of` map (size nz*ny*nx; node id
// per voxel or -1) and `voxel_of` map (size g->nnodes; voxel linear index per
// node). Returns 0 on success. Free with sgraph_free + free(node_of)/free(voxel_of).
int  affinity_build_voxel(const u8 *mask, const f32 *sheet, const f32 *normal,
                          int nz, int ny, int nx, const affinity_params *p,
                          sgraph *g, s32 **node_of, u32 **voxel_of);
void sgraph_free(sgraph *g);

#endif // TABERNA_AFFINITY_H
