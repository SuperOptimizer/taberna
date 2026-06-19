/* morph.h — classical binary post-processing for surface masks.
 *
 * These are the operations that won the Kaggle Surface-Detection competition (the
 * segmentation model was commoditized; the score lived here — see
 * docs/surface-detection-kaggle.md). All operate on u8 masks (nonzero = fg),
 * z-major / x-fastest: idx = (z*ny + y)*nx + x.
 *
 * The headline primitive is the iterated 3x3x3 binary MEDIAN filter, which is
 * exactly a 27-neighborhood MAJORITY vote (a voxel turns on iff >= 14 of its 27
 * neighbors are on). Iterated, it is a discrete surface-tension / curvature flow:
 * it shaves thin spurs and bridges (the source of Betti-1 tunnels and false sheet
 * merges) and closes pinholes, while leaving flat sheet interior untouched — so it
 * raises the topology score without costing surface Dice. In the competition this
 * single step (run ~6-10x) beat binary closing, hole-patching and cavity fill (the
 * 10th/18th-place "few lines of code" boost, ~+0.01 private LB).
 */
#ifndef TABERNA_MORPH_H
#define TABERNA_MORPH_H

#include <stddef.h>
#include "common/types.h"
#include "eval/topo.h"   // topo_conn

/* 3x3x3 binary median == 27-neighborhood majority (threshold 14), applied `iters`
 * times. out-of-bounds neighbors count as background. `in` and `out` may alias the
 * same buffer; internal double-buffering handles it. out must hold nz*ny*nx. */
void majority_filter(const u8 *in, u8 *out, int nz, int ny, int nx, int iters);

/* General majority with an explicit on-threshold over the 27-neighborhood (14 =
 * median; lower = more dilating, higher = more eroding). */
void majority_filter_thresh(const u8 *in, u8 *out, int nz, int ny, int nx,
                            int iters, int on_thresh);

/* Zero foreground components smaller than `min_voxels` ("dust removal" — the
 * single biggest post-proc gain in the 1st-place ablation). In-place. Returns the
 * number of components removed. */
size_t remove_small_components(u8 *mask, int nz, int ny, int nx,
                               topo_conn conn, size_t min_voxels);

/* Fill enclosed background cavities (6-conn background not reachable from the
 * volume border becomes foreground) — i.e. drive b2 to 0. In-place. */
void fill_holes(u8 *mask, int nz, int ny, int nx);

/* Plug 1-voxel pinholes: a background voxel whose 6 face-neighbors are all
 * foreground becomes foreground. Cheap complement to fill_holes for the common
 * single-voxel case. In-place. Returns voxels plugged. */
size_t plug_pinholes(u8 *mask, int nz, int ny, int nx);

/* Morphology with a Euclidean-ball structuring element of radius `r` (voxels with
 * dz^2+dy^2+dx^2 <= r^2). `in`/`out` may alias. */
void dilate_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r);
void erode_ball (const u8 *in, u8 *out, int nz, int ny, int nx, int r);
void closing_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r); // dilate∘erode
void opening_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r); // erode∘dilate

#endif // TABERNA_MORPH_H
