/* topo_surgery.h — PH-guided tunnel filling (dim-1 topological surgery).
 *
 * Uses the native cubical PH (cubical_features) to LOCALIZE each foreground tunnel
 * (its representative H1 loop) and fills it, instead of the blunt global dilation
 * that destroys SurfaceDice. This is the #2 (dim-1 PH) work applied to make the
 * surface tracer's output watertight (b1 -> 0): trace gives connected-but-porous
 * sheets, this caps the residual tunnels precisely.
 *
 * The generic reduction is O(cells), so we tile into small windows and run PH only
 * in windows a fast Betti screen flags as having b1>0; the tunnel loop's voxels are
 * then dilate-filled (radius `radius`) into the mask, bridging the small hole
 * without thickening the flat sheet. Repeated `iters` times. Volumes z-major.
 */
#ifndef TABERNA_TOPO_SURGERY_H
#define TABERNA_TOPO_SURGERY_H

#include "common/types.h"

/* Fill foreground tunnels in `mask` in place. `win`/`overlap` set the PH window
 * tiling; `radius` the local fill ball; `iters` repeat count. Returns the total
 * number of tunnel windows processed. */
int fill_tunnels(u8 *mask, int nz, int ny, int nx,
                 int win, int overlap, int radius, int iters);

#endif // TABERNA_TOPO_SURGERY_H
