/* ridge.h — thin-sheet centerline extraction by non-maximum suppression along the
 * structure-tensor normal.
 *
 * Raw structure-tensor sheetness is a THICK response over the whole body of a
 * papyrus sheet, but the surface-detection label is the sheet's thin smoothed
 * centerline (the bright medial ridge). This recovers that ridge classically: a
 * voxel is kept iff (a) it is sheet-like enough and bright enough, and (b) its CT
 * intensity is a local maximum along the across-sheet normal n (sampling CT at
 * x ± step·n by trilinear interpolation). That non-max suppression collapses the
 * thick sheetness band to a ~1-voxel ridge that matches the label convention.
 *
 * Volumes are (nz,ny,nx), z-major, x-fastest. `normal` is the interleaved output
 * of st_sheet_detect: 3 components per voxel ordered (x,y,z); its sign is
 * arbitrary (an axis), which is fine since NMS samples both ±n.
 */
#ifndef TABERNA_RIDGE_H
#define TABERNA_RIDGE_H

#include <stddef.h>
#include "common/types.h"

/* Write a binary ridge mask into `out` (nz*ny*nx). Keep voxel i iff
 * sheet[i] >= s_min, ct[i] >= i_min, and ct[i] is >= the trilinear CT samples at
 * x ± step·normal[i]. Returns the number of voxels kept. */
size_t ridge_nms(const f32 *ct, const f32 *sheet, const f32 *normal,
                 int nz, int ny, int nx,
                 f32 s_min, f32 i_min, f32 step, u8 *out);

#endif // TABERNA_RIDGE_H
