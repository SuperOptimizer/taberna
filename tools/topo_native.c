/* topo_native.c — print taberna's native TopoScore for two binary mask TIFFs,
 * for validation against the official topometrics Betti-Matching TopoScore.
 * Usage: topo_native PRED.tif GT.tif */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "eval/score.h"

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: %s PRED.tif GT.tif\n", argv[0]); return 2; }
  int nz,ny,nx,gz,gy,gx;
  u8 *p = tiff_load_u8(argv[1], &nz,&ny,&nx);
  u8 *g = tiff_load_u8(argv[2], &gz,&gy,&gx);
  if (!p || !g || gz!=nz || gy!=ny || gx!=nx) { fprintf(stderr,"load/dim error\n"); return 1; }
  printf("%.6f\n", toposcore_native(p, g, nz, ny, nx));
  free(p); free(g);
  return 0;
}
