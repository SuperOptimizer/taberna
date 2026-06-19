/* tiff_vol.c — see tiff_vol.h. */
#include "io/tiff_vol.h"

#include <stdlib.h>
#include "tiff.h"   // vendored third-party/tiff

u8 *tiff_load_u8(const char *path, int *nz, int *ny, int *nx) {
  tiff_volume v;
  if (tiff_read_volume(path, &v) != TIFF_OK) return NULL;
  if (v.channels != 1 || v.type != TIFF_U8) { tiff_volume_free(&v); return NULL; }
  if (nz) *nz = (int)v.depth;
  if (ny) *ny = (int)v.height;
  if (nx) *nx = (int)v.width;
  u8 *data = (u8 *)v.data;   // transfer ownership; already z-major u8
  v.data = NULL;
  tiff_volume_free(&v);
  return data;
}
