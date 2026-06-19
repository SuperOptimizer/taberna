/* surface_clean.c — detect once, sweep CLOSING-based cleanup vs the accurate
 * metric. Closing fills holes (drops b1) while the erosion restores thickness, so
 * it should cut tunnels with less SurfaceDice loss than dilation. */
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
  if(!img||!lab||lz!=nz){fprintf(stderr,"load/dim err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32*sheet=malloc(n*sizeof(f32)),*normal=malloc(3*n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f};
  fprintf(stderr,"structure tensor...\n");
  st_sheet_detect(volf,nz,ny,nx,&sp,sheet,normal);
  u8*ridge=malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,0.3f,80.f,1.0f,ridge);
  remove_small_components(ridge,nz,ny,nx,TOPO_CONN26,100);
  free(volf);free(sheet);free(normal);free(img);

  u8*q=malloc(n);
  printf(" close med   Score   Topo   SurfD   VOI    b1\n");
  printf("-----------------------------------------------------\n");
  int rs[]={2,3,4}, ms[]={0,2,4};
  double best=0; int br=0,bm=0;
  for(size_t a=0;a<sizeof(rs)/sizeof(*rs);a++)
   for(size_t b=0;b<sizeof(ms)/sizeof(*ms);b++){
     memcpy(q,ridge,n);
     closing_ball(q,q,nz,ny,nx,rs[a]);
     if(ms[b]) majority_filter(q,q,nz,ny,nx,ms[b]);
     fill_holes(q,nz,ny,nx);
     comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
     printf("   %d   %d   %.4f  %.4f  %.4f  %.4f  %5ld\n",
            rs[a],ms[b],s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b1);
     fflush(stdout);
     if(s.score>best){best=s.score;br=rs[a];bm=ms[b];}
   }
  printf("\nbest %.4f at close=%d median=%d\n",best,br,bm);
  free(ridge);free(q);free(lab);return 0;
}
