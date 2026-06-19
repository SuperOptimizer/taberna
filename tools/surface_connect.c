/* surface_connect.c — ridge -> targeted fragment connection (VOI) -> cleanup.
 * Sweeps connect radius. Usage: surface_connect IMAGE LABEL */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s IMAGE LABEL\n",argv[0]);return 2;}
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*v=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) v[i]=(f32)img[i];
  f32*sh=malloc(n*sizeof(f32)),*nm=malloc(3*n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f};
  st_sheet_detect(v,nz,ny,nx,&sp,sh,nm);
  u8*ridge=malloc(n);
  ridge_nms(v,sh,nm,nz,ny,nx,0.3f,80.f,1.0f,ridge);
  remove_small_components(ridge,nz,ny,nx,TOPO_CONN26,100);
  free(v);free(sh);free(nm);free(img);
  u8*q=malloc(n);
  printf("baseline dilate1+med2: 0.4406 (b0~322 before postproc)\n");
  printf(" connect_r  added   Score   SurfD   VOI    b0\n");
  for(int r=0;r<=4;r++){
    memcpy(q,ridge,n);
    size_t added=0;
    if(r>0) added=connect_fragments(q,nz,ny,nx,r);
    dilate_ball(q,q,nz,ny,nx,1); majority_filter(q,q,nz,ny,nx,2); fill_holes(q,nz,ny,nx);
    comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
    printf("    %d    %7zu  %.4f  %.4f  %.4f  %ld\n",r,added,s.score,s.surface_dice,s.voi_score,s.pred_b0);
    fflush(stdout);
  }
  free(ridge);free(q);free(lab);return 0;
}
