/* trace.h — sheet surface tracing by advancing front (the watertight detector).
 *
 * Post-processing a porous ridge can't reach b1=0 (gaps in the data become
 * tunnels). This grows each sheet as a connected surface instead: from a seed,
 * an advancing front marches in the local tangent plane (perpendicular to the
 * structure-tensor normal) and at each step SNAPS back onto the sheet ridge
 * (intensity max along the normal). Marching *over* a gap and re-acquiring the
 * sheet on the far side BRIDGES porosity, and rasterizing the path keeps the
 * surface connected — so the traced sheet is watertight (b1=0) where the porous
 * ridge had holes. The normal-consistency gate stops the front jumping to an
 * adjacent wrap. This is taberna's surface tracer and the front end of unrolling.
 *
 * Volumes z-major, x-fastest.
 */
#ifndef TABERNA_TRACE_H
#define TABERNA_TRACE_H

#include "common/types.h"
#include "segmentation/sheet_tensor.h"

typedef struct {
  st_params st;          // structure-tensor scales (sheetness + normal)
  f32 seed_thresh;       // min sheetness to START a sheet (strict, e.g. 0.4)
  f32 cont_thresh;       // min sheetness to CONTINUE/bridge (permissive, e.g. 0.08)
  f32 step;              // in-plane march step in voxels (e.g. 1.0)
  f32 snap_radius;       // search range along normal to re-acquire the ridge (e.g. 2.0; < inter-wrap gap)
  f32 i_min;             // min CT intensity for a ridge point (e.g. 80)
  f32 normal_cos;        // min |n . n'| to accept a step (e.g. 0.6 ~ 53deg)
  int min_size;          // discard sheets smaller than this many voxels
} trace_params;

trace_params trace_default_params(void);

/* Trace watertight sheets in `vol` (CT, f32). Writes the union of traced sheets
 * (u8 0/1) into `out` (nz*ny*nx, zeroed first). Returns the number of sheets, or
 * <0 on allocation failure. */
int sheet_trace(const f32 *vol, int nz, int ny, int nx, const trace_params *p, u8 *out);

#endif // TABERNA_TRACE_H
