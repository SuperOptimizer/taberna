/* surface_hess.c — Hessian (Frangi) multi-scale plate detector -> ridge -> score.
 * Usage: surface_hess IMAGE LABEL [s_min=0.3] [i_min=80] [dilate=2] [median=4] [beta_c=15] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/hessian.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s IMAGE LABEL [s_min][i_min][dil][med][beta_c]\n",argv[0]);return 2;}
  float s_min=argc>3?atof(argv[3]):0.3f,i_min=argc>4?atof(argv[4]):80.f;
  int dil=argc>5?atoi(argv[5]):2,med=argc>6?atoi(argv[6]):4; float bc=argc>7?atof(argv[7]):15.f;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32*sheet=malloc(n*sizeof(f32)),*normal=malloc(3*n*sizeof(f32));
  f32 sig[3]={1.0f,2.0f,3.0f}; hess_params hp={sig,3,0.5f,bc};
  fprintf(stderr,"hessian multi-scale...\n");
  hessian_sheet_detect(volf,nz,ny,nx,&hp,sheet,normal);
  // sheetness stats
  double mx=0; for(size_t i=0;i<n;i++) if(sheet[i]>mx)mx=sheet[i];
  fprintf(stderr,"max sheetness=%.3f\n",mx);
  u8*q=malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,1.0f,q);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0)dilate_ball(q,q,nz,ny,nx,dil);
  if(med>0)majority_filter(q,q,nz,ny,nx,med);
  fill_holes(q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("HESSIAN: Score %.4f  Topo %.4f  SurfD %.4f  VOI %.4f  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(volf);free(sheet);free(normal);free(q);free(img);free(lab);return 0;
}
