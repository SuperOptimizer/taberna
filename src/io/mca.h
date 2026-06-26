/* mca.h — read regions of a matter-compressor archive (.mca) into taberna's
 * u8 z-major volumes. Consumes matter-compressor's read API only (we never modify
 * mc). For running the pipeline on REAL scroll data (e.g. the full PHerc0332.mca).
 */
#ifndef TABERNA_MCA_H
#define TABERNA_MCA_H

#include <stddef.h>
#include "common/types.h"

/* Archive geometry (LOD0 voxel dims + build quality + #LODs). Returns 0 on success. */
int mca_dims(const char *path, int *nz, int *ny, int *nx, float *quality, int *nlods);

/* Persistent archive handle: open the .mca ONCE (one mmap + one node-table parse)
 * and read many regions through it, instead of re-opening per read. Use this for
 * multi-region jobs (tiling the whole scroll) — the per-read functions below each
 * re-open the archive. */
typedef struct mca_handle mca_handle;
mca_handle *mca_open(const char *path);
void        mca_close(mca_handle *h);
void        mca_handle_dims(const mca_handle *h, int *nz, int *ny, int *nx, float *quality, int *nlods);
/* Read a (dz,dy,dx) region at (z0,y0,x0) of `lod` through an open handle into a
 * malloc'd u8 buffer (z-major, x-fastest). Returns NULL on failure. */
u8 *mca_read(const mca_handle *h, int lod, int z0, int y0, int x0, int dz, int dy, int dx);

/* Read a (dz,dy,dx) region at (z0,y0,x0) of `lod` into a malloc'd u8 buffer
 * (z-major, x-fastest). Returns NULL on failure. */
u8 *mca_load_region(const char *path, int lod, int z0, int y0, int x0,
                    int dz, int dy, int dx);

/* Borrow the archive's provenance metadata blob (the JSON stamped into the .mca
 * carveout by mca_export). Returns a pointer owned by the handle (do NOT free) and
 * its length in *out_len, or NULL if the archive carries no metadata. */
const char *mca_metadata(const mca_handle *h, size_t *out_len);

/* Read the ROI trim origin from the archive metadata ("roi":{"origin":[z,y,x]}):
 * the world/source-volume voxel that the archive's (0,0,0) maps to. Coords from a
 * full-volume frame (e.g. VC3D segments) become archive coords by subtracting it.
 * Writes *oz,*oy,*ox (at LOD0). Returns 0 on success, <0 if absent/unparseable
 * (in which case origin is left at 0,0,0 = assume no trim). */
int mca_roi_origin(const mca_handle *h, int *oz, int *oy, int *ox);

/* Find one material-rich region of size (d,d,d) at `lod`: deterministic seeded
 * sampling keeping boxes with mean material fraction >= min_frac. Writes the
 * chosen origin to *z0,*y0,*x0. Returns 0 on success, <0 if none found. */
int mca_find_region(const char *path, int lod, int d, float min_frac, uint64_t seed,
                    int *z0, int *y0, int *x0);

#endif // TABERNA_MCA_H
