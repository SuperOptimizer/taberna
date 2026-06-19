#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "eval/nsd.h"
int main(int argc,char**argv){
  if(argc!=4){fprintf(stderr,"usage: %s GT.tif PRED.tif tol\n",argv[0]);return 2;}
  int nz,ny,nx,pz,py,px;
  u8*g=tiff_load_u8(argv[1],&nz,&ny,&nx),*p=tiff_load_u8(argv[2],&pz,&py,&px);
  if(!g||!p||pz!=nz||py!=ny||px!=nx){fprintf(stderr,"load/dim err\n");return 1;}
  printf("%.6f\n",surface_dice_nsd(g,p,nz,ny,nx,atof(argv[3])));
  return 0;
}
