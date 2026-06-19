/* extract_cube — crop a (z0,y0,x0)+(dz,dy,dx) subcube out of a NRRD volume and
 * write it as a NRRD. Use it to carve small test cubes out of the big public
 * volumes (e.g. instance-labels-harmonized) for the experiment phase.
 *
 *   extract_cube <in.nrrd> <out.nrrd> <z0> <y0> <x0> <dz> <dy> <dx>
 */
#include "io/nrrd.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  if (argc != 9) {
    fprintf(stderr, "usage: %s <in.nrrd> <out.nrrd> z0 y0 x0 dz dy dx\n", argv[0]);
    return 2;
  }
  int z0 = atoi(argv[3]), y0 = atoi(argv[4]), x0 = atoi(argv[5]);
  int dz = atoi(argv[6]), dy = atoi(argv[7]), dx = atoi(argv[8]);

  int nz, ny, nx;
  f32 *vol = nrrd_read_f32(argv[1], &nz, &ny, &nx);
  if (!vol) { fprintf(stderr, "read failed: %s\n", argv[1]); return 1; }
  if (z0 < 0 || y0 < 0 || x0 < 0 || z0 + dz > nz || y0 + dy > ny || x0 + dx > nx) {
    fprintf(stderr, "subcube out of bounds (volume is %dx%dx%d z,y,x)\n", nz, ny, nx);
    return 1;
  }

  f32 *out = malloc((size_t)dz * dy * dx * sizeof(f32));
  for (int z = 0; z < dz; z++)
    for (int y = 0; y < dy; y++)
      for (int x = 0; x < dx; x++)
        out[((size_t)z * dy + y) * dx + x] =
            vol[((size_t)(z0 + z) * ny + (y0 + y)) * nx + (x0 + x)];

  int rc = nrrd_write_f32(argv[2], out, dz, dy, dx);
  if (rc) fprintf(stderr, "write failed: %s\n", argv[2]);
  else printf("wrote %s : %d x %d x %d (z,y,x)\n", argv[2], dz, dy, dx);
  free(vol); free(out);
  return rc ? 1 : 0;
}
