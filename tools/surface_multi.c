/* surface_multi.c — multi-cue detector: structure-tensor sheetness * Hessian
 * plate-ness. Both fire on real sheets, disagree on noise -> the product is a
 * cleaner response (fewer false fragments -> better VOI). Reports VOI explicitly.
 * Usage: surface_multi IMAGE LABEL [s_min=0.25] [i_min=80] [dil=2] [med=4] */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/hessian.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s IMAGE LABEL [s_min][i_min][dil][med]\n",argv[0]);return 2;}
  float s_min=argc>3?atof(argv[3]):0.25f,i_min=argc>4?atof(argv[4]):80.f;
  int dil=argc>5?atoi(argv[5]):2,med=argc>6?atoi(argv[6]):4;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32*sst=malloc(n*sizeof(f32)),*nst=malloc(3*n*sizeof(f32)),*sh=malloc(n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f};
  fprintf(stderr,"structure tensor...\n"); st_sheet_detect(volf,nz,ny,nx,&sp,sst,nst);
  f32 sig[3]={1,2,3}; hess_params hp={sig,3,0.5f,15.f};
  fprintf(stderr,"hessian...\n"); hessian_sheet_detect(volf,nz,ny,nx,&hp,sh,NULL);
  for(size_t i=0;i<n;i++) sst[i]=sqrtf(sst[i]*sh[i]);   // geometric-mean combine
  u8*q=malloc(n);
  ridge_nms(volf,sst,nst,nz,ny,nx,s_min,i_min,1.0f,q);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0)dilate_ball(q,q,nz,ny,nx,dil);
  if(med>0)majority_filter(q,q,nz,ny,nx,med);
  fill_holes(q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("MULTI-CUE: Score %.4f  Topo %.4f  SurfD %.4f  VOI %.4f(=%.3f bits)  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.voi,s.pred_b0,s.pred_b1);
  free(volf);free(sst);free(nst);free(sh);free(q);free(img);free(lab);return 0;
}
