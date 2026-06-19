/* mca.h — read regions of a matter-compressor archive (.mca) into taberna's
 * u8 z-major volumes. Consumes matter-compressor's read API only (we never modify
 * mc). For running the pipeline on REAL scroll data (e.g. the full PHerc0332.mca).
 */
#ifndef TABERNA_MCA_H
#define TABERNA_MCA_H

#include "common/types.h"

/* Archive geometry (LOD0 voxel dims + build quality + #LODs). Returns 0 on success. */
int mca_dims(const char *path, int *nz, int *ny, int *nx, float *quality, int *nlods);

/* Read a (dz,dy,dx) region at (z0,y0,x0) of `lod` into a malloc'd u8 buffer
 * (z-major, x-fastest). Returns NULL on failure. */
u8 *mca_load_region(const char *path, int lod, int z0, int y0, int x0,
                    int dz, int dy, int dx);

/* Find one material-rich region of size (d,d,d) at `lod`: deterministic seeded
 * sampling keeping boxes with mean material fraction >= min_frac. Writes the
 * chosen origin to *z0,*y0,*x0. Returns 0 on success, <0 if none found. */
int mca_find_region(const char *path, int lod, int d, float min_frac, uint64_t seed,
                    int *z0, int *y0, int *x0);

#endif // TABERNA_MCA_H
