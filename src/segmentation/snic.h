/* snic.h — 3D SNIC supervoxels for turning a dense scroll volume into a sparse
 * region-adjacency graph.
 *
 * SNIC (Simple Non-Iterative Clustering, Achanta & Susstrunk 2017) grows
 * compact superpixels from a seed grid using a single priority-queue sweep —
 * no iteration, deterministic, O(N). This is the 3D variant: it partitions a
 * dense (nz,ny,nx) scalar volume into ~N compact supervoxels and returns, for
 * each, a centroid, a mean value, voxel counts, and a list of adjacent
 * supervoxels. That adjacency list IS the graph: the first step in collapsing
 * the dense voxel grid into something sparse enough to fit + reason about
 * sheets, wraps, and (eventually) the scroll's spiral on.
 *
 * --- Pipeline role / structure-tensor coupling --------------------------------
 * `img` is just a scalar field. Feed it whatever scalar best exposes the
 * structure you want supervoxels to respect:
 *   - the raw CT volume (intensity-compact blobs), or
 *   - mc_segment's structure-tensor *sheetness* field `(l0-l1)/(l0+eps)`
 *     (supervoxels that hug papyrus sheets).
 * Per-supervoxel sheet *orientation* (the dominant tensor eigenvector), which
 * the spiral-undeformation work will want on each graph node, is a planned
 * extension — see the note by `snic()` in snic.c. The distance metric here is
 * the classic intensity+spatial one, unchanged from the reference.
 *
 * --- Memory layout (IMPORTANT) ------------------------------------------------
 * Re-indexed from the stabia reference (which was y-fastest) to taberna's
 * convention, matching matter-compressor's mc_segment: volumes are z-major,
 * x-fastest, voxel (z,y,x) at  v[(z*ny + y)*nx + x].
 *
 * Ported from https://github.com/spelufo/stabia (MIT, (c) 2023 Santiago
 * Pelufo), itself based on the reference SNIC implementation
 * (https://github.com/achanta/SNIC). See LICENSE note in snic.c.
 */
#ifndef TABERNA_SNIC_H
#define TABERNA_SNIC_H

#include <stdint.h>

typedef uint8_t  u8;
typedef int8_t   s8;
typedef uint16_t u16;
typedef int16_t  s16;
typedef uint32_t u32;
typedef int32_t  s32;
typedef uint64_t u64;
typedef int64_t  s64;
typedef float    f32;
typedef double   f64;

// There is no theoretical maximum for SNIC neighbors; a cube has 26, so with
// reasonable compactness we shouldn't exceed that by much. 56*2 keeps
// sizeof(Superpixel) a multiple of a cacheline.
#define SUPERPIXEL_MAX_NEIGHS (56 * 2)

// One graph node: a supervoxel. After snic() returns, (x,y,z) is the centroid,
// `c` the mean `img` value, `n` the voxel count, nlow/nmid/nhig the counts of
// voxels falling in the [<=lowmid], (lowmid,midhig], (>midhig] value buckets
// (handy for air / papyrus / dense-material splits), and `neighs` the ids of
// adjacent supervoxels (0-terminated, 0 is unused — labels are 1-based).
typedef struct Superpixel {
  f32 x, y, z, c;
  u32 n, nlow, nmid, nhig;
  u32 neighs[SUPERPIXEL_MAX_NEIGHS];
} Superpixel;

// Capacity of the neighbor list, for callers walking the graph.
int snic_superpixel_max_neighs(void);

// Number of supervoxels snic() will produce for a (nz,ny,nx) volume seeded with
// step `d_seed`. Use it to size the `superpixels` array (+1 for the unused 0
// slot).
int snic_superpixel_count(int nz, int ny, int nx, int d_seed);

// Append k2 to k1's neighbor list (idempotent). Returns 1 if the list was full
// (neighbor dropped), else 0. Exposed mainly for testing.
int superpixel_add_neighbors(Superpixel *superpixels, u32 k1, u32 k2);

// Partition `img` into compact supervoxels.
//   img          dense (nz,ny,nx) scalar field, v[(z*ny+y)*nx+x]
//   d_seed       seed-grid step (≈ target supervoxel edge length in voxels)
//   compactness  spatial-vs-value weight; higher => rounder, more uniform sizes
//   lowmid,midhig value-bucket thresholds for nlow/nmid/nhig accounting
//   labels       OUT, must be nz*ny*nx and zero-initialized; gets the 1-based
//                supervoxel id of each voxel
//   superpixels  OUT, must hold snic_superpixel_count()+1 entries, zeroed;
//                index by label (slot 0 unused)
// Returns the number of adjacencies that overflowed SUPERPIXEL_MAX_NEIGHS (0 is
// the healthy case).
int snic(f32 *img, int nz, int ny, int nx, int d_seed, f32 compactness,
         f32 lowmid, f32 midhig, u32 *labels, Superpixel *superpixels);

#endif // TABERNA_SNIC_H
