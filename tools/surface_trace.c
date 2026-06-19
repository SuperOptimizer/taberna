/* surface_trace.c — watertight sheet tracing -> score. */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/trace.h"
#include "postproc/morph.h"
#include "eval/score.h"
#include "eval/topo.h"
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s IMAGE LABEL OUT.tif [seed=0.3][snap=2][min=200]\n",argv[0]);return 2;}
  float seed=argc>4?atof(argv[4]):0.3f, snap=argc>5?atof(argv[5]):2.0f; int minsz=argc>6?atoi(argv[6]):200;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"load/dim err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i]; free(img);
  trace_params tp=trace_default_params(); tp.seed_thresh=seed; tp.snap_radius=snap; tp.min_size=minsz;
  fprintf(stderr,"tracing...\n");
  u8*q=malloc(n);
  int ns=sheet_trace(volf,nz,ny,nx,&tp,q);
  fprintf(stderr,"%d fronts\n",ns);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,minsz);
  fill_holes(q,nz,ny,nx);
  tiff_save_u8(argv[3],q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("Score %.4f  Topo %.4f  SurfD %.4f  VOI %.4f  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(volf);free(q);free(lab);return 0;
}
