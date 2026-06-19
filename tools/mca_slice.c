/* mca_slice — dump one z-slice of an .mca at a given lod/region as a TIFF (diagnostic). */
#include <stdio.h>
#include <stdlib.h>
#include "io/mca.h"
#include "io/tiff_vol.h"
int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: %s ARC OUT.tif lod z0 y0 x0 d\n",argv[0]);return 2;}
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),d=atoi(argv[7]);
  u8*img=mca_load_region(path,lod,z0,y0,x0,1,d,d);
  if(!img){fprintf(stderr,"read fail\n");return 1;}
  tiff_save_u8(outp,img,1,d,d); printf("wrote %dx%d lod%d @ (%d,%d,%d)\n",d,d,lod,z0,y0,x0);
  return 0;
}
