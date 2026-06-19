/* surface_optimize.c — thorough EXACT-metric optimization. Detect sheetness once,
 * then sweep single-threshold postproc grids AND multi-threshold consensus (the
 * 2nd-place method: vote across thresholds). Reports the best exact score. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
static comp_score best; static char bestcfg[128];
static void try_cfg(u8*q,u8*lab,int d,const char*cfg){
  comp_score s=competition_score(q,lab,d,d,d,2,1,2);
  printf("  %-28s Score %.4f  Topo %.4f SurfD %.4f VOI %.4f b0=%ld b1=%ld\n",
         cfg,s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1); fflush(stdout);
  if(s.score>best.score){best=s; strncpy(bestcfg,cfg,127);}
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage: %s IMAGE LABEL\n",argv[0]);return 2;}
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"err\n");return 1;}
  int d=nz; size_t n=(size_t)d*d*d;
  f32*v=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) v[i]=(f32)img[i];
  f32*sh=malloc(n*sizeof(f32)),*nm=malloc(3*n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f}; fprintf(stderr,"detect...\n"); st_sheet_detect(v,d,d,d,&sp,sh,nm);
  u8*q=malloc(n),*r=malloc(n); best.score=-1;
  float ths[]={0.2f,0.3f,0.4f,0.5f};
  // single-threshold grid
  for(size_t t=0;t<4;t++){
    ridge_nms(v,sh,nm,d,d,d,ths[t],80.f,1.0f,r);
    remove_small_components(r,d,d,d,TOPO_CONN26,100);
    int dils[]={1,2}, meds[]={2,4,6};
    for(int a=0;a<2;a++)for(int b=0;b<3;b++){
      memcpy(q,r,n); dilate_ball(q,q,d,d,d,dils[a]); majority_filter(q,q,d,d,d,meds[b]); fill_holes(q,d,d,d);
      char c[128]; snprintf(c,128,"t%.1f dil%d med%d",ths[t],dils[a],meds[b]); try_cfg(q,lab,d,c);
    }
  }
  // multi-threshold consensus: vote (>=2 of 3 thresholds) -> cleaner, then postproc
  { u8*acc=calloc(n,1);
    for(size_t t=0;t<3;t++){ ridge_nms(v,sh,nm,d,d,d,ths[t],80.f,1.0f,r);
      for(size_t i=0;i<n;i++) acc[i]+=r[i]!=0; }
    for(size_t i=0;i<n;i++) q[i]= acc[i]>=2;
    remove_small_components(q,d,d,d,TOPO_CONN26,100);
    dilate_ball(q,q,d,d,d,2); majority_filter(q,q,d,d,d,4); fill_holes(q,d,d,d);
    try_cfg(q,lab,d,"consensus>=2 dil2 med4"); free(acc); }
  printf("\nBEST: %.4f  [%s]\n",best.score,bestcfg);
  return 0;
}
