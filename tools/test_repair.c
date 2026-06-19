/* test_repair.c — sheet_repair must turn a holed/tunneled sheet watertight (b1=0)
 * while keeping it a single sheet. */
#include <stdio.h>
#include <stdlib.h>
#include "postproc/sheet_repair.h"
#include "eval/topo.h"

int main(void){
  int nz=20,ny=40,nx=40; size_t n=(size_t)nz*ny*nx;
  u8*m=calloc(n,1),*o=calloc(n,1);
  // a flat-ish sheet at z=10 spanning y,x in [4,35], thickness 2, with a hole
  for(int y=4;y<36;y++)for(int x=4;x<36;x++){ m[((size_t)10*ny+y)*nx+x]=1; m[((size_t)11*ny+y)*nx+x]=1; }
  // punch a hole (tunnel through the sheet) at center
  for(int z=0;z<nz;z++)for(int y=18;y<22;y++)for(int x=18;x<22;x++) m[((size_t)z*ny+y)*nx+x]=0;
  topo_betti before=betti_numbers_6(m,nz,ny,nx);
  printf("before repair: b0=%ld b1=%ld b2=%ld\n",before.b0,before.b1,before.b2);
  sheet_repair(m,o,nz,ny,nx,50);
  topo_betti after=betti_numbers_6(o,nz,ny,nx);
  printf("after  repair: b0=%ld b1=%ld b2=%ld  (want b0=1 b1=0)\n",after.b0,after.b1,after.b2);
  int ok = (after.b0==1 && after.b1==0);
  printf(ok?"PASS\n":"FAIL\n");
  free(m);free(o);
  return ok?0:1;
}
