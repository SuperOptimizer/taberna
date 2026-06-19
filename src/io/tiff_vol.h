/* tiff_vol.h — load a multi-page TIFF volume into taberna's u8 z-major layout.
 *
 * Thin wrapper over the vendored tiff submodule's tiff_read_volume (LZW + raw,
 * multi-IFD). The Vesuvius surface-detection data ships as 320^3 u8 TIFF stacks
 * (CT image 0..255; label 0=bg, 1=surface, 2=ignore). Returned buffer is
 * z-major / x-fastest: v[(z*ny + y)*nx + x]. Caller frees with free().
 */
#ifndef TABERNA_TIFF_VOL_H
#define TABERNA_TIFF_VOL_H

#include "common/types.h"

/* Load a single-channel u8 TIFF volume. On success returns a malloc'd buffer of
 * nz*ny*nx bytes and sets the dims; on failure returns NULL. Rejects multi-channel
 * or non-u8 volumes (the surface-detection data is neither). */
u8 *tiff_load_u8(const char *path, int *nz, int *ny, int *nx);

#endif // TABERNA_TIFF_VOL_H
