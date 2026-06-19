/* sheet_repair.h — make each sheet watertight by fitting it as a height-map
 * surface over its principal plane (the 1st/2nd/4th-place classical repair).
 *
 * The metric demands topologically clean sheets (GT has b1=0), but our ridge
 * detector leaves hundreds of tunnels that no amount of median/dilate removes.
 * This rebuilds each connected component as a single-valued height surface
 * w(u,v) over its PCA tangent plane:
 *   PCA(component) -> (tangent1=u, tangent2=v, normal=w);
 *   splat heights onto a (u,v) grid; flood-remove the outside-the-outline cells;
 *   Laplace-inpaint the interior gaps; rasterize the surface watertight.
 * A height surface over a simply-connected (u,v) domain is a disk (b1=0, b2=0),
 * so holes and tunnels are filled by construction. Same SVD->(u,v,w) primitive as
 * the Archimedean unwrap.
 *
 * Volumes z-major, x-fastest. Caveat: assumes each component is a height function
 * over one plane (true for local patches; strongly-folded components degrade).
 */
#ifndef TABERNA_SHEET_REPAIR_H
#define TABERNA_SHEET_REPAIR_H

#include "common/types.h"

/* Rebuild every >= min_voxels component of `in` as a watertight height surface,
 * writing the union into `out` (zeroed first; `in` and `out` must differ). */
void sheet_repair(const u8 *in, u8 *out, int nz, int ny, int nx, int min_voxels);

/* Windowed repair: run sheet_repair on overlapping `win`-sized blocks (stride
 * `win-overlap`) and OR the results, so each block's sheet patch is locally planar
 * (curved sheets break the global single-plane fit). out zeroed first. */
void sheet_repair_windowed(const u8 *in, u8 *out, int nz, int ny, int nx,
                           int win, int overlap, int min_voxels);

#endif // TABERNA_SHEET_REPAIR_H
