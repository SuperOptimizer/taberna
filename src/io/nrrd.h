/* nrrd.h — minimal NRRD (Nearly Raw Raster Data) reader/writer for taberna's
 * experiment phase. NRRD is a trivial ASCII-header + binary-blob format; the
 * friction-free Vesuvius feasibility data (instance-labels-harmonized, the GP
 * cubes) ships as .nrrd. We support what that data needs: 3D volumes, encodings
 * `raw` and `gzip`, sample types uint8 / uint16 / float (little-endian).
 *
 * Layout convention: NRRD stores the first listed axis fastest. We treat a 3D
 * volume's axes as (x, y, z) fastest-first, which matches taberna's z-major,
 * x-fastest convention used elsewhere (voxel (z,y,x) at v[(z*ny+y)*nx+x]).
 */
#ifndef TABERNA_NRRD_H
#define TABERNA_NRRD_H

#include <stddef.h>
#include "common/types.h"

typedef enum { NRRD_U8, NRRD_U16, NRRD_F32 } nrrd_type;

typedef struct {
  int       ndim;        // number of dimensions (we use/expect 3)
  int       sizes[4];    // axis sizes, fastest-first: sizes[0]=nx, [1]=ny, [2]=nz
  nrrd_type type;        // sample type
  int       type_size;   // bytes per sample (1/2/4)
  size_t    count;       // total number of samples
  void     *data;        // raw little-endian samples (owned)
} nrrd;

// Read a NRRD file. Returns 0 on success (fills *out, caller frees with
// nrrd_free), non-zero on error.
int  nrrd_read(const char *path, nrrd *out);
void nrrd_free(nrrd *n);

// Convenience: load a 3D NRRD and return it as a freshly-allocated f32 volume in
// taberna layout (v[(z*ny+y)*nx+x]). Sets *nz,*ny,*nx. uint8/uint16/float inputs
// are cast to f32 (no rescaling). Returns NULL on error. Caller free()s the result.
f32 *nrrd_read_f32(const char *path, int *nz, int *ny, int *nx);

// Write a 3D volume (taberna layout) as a `raw`-encoded NRRD.
// Returns 0 on success.
int  nrrd_write_f32(const char *path, const f32 *vol, int nz, int ny, int nx);
int  nrrd_write_u8 (const char *path, const u8  *vol, int nz, int ny, int nx);

#endif // TABERNA_NRRD_H
