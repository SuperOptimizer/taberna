/* surface_opt.c — optimize the prediction against the native composite score.
 * Detect once (structure tensor + ridge), then sweep topology-cleanup configs,
 * scoring each with the native metric (VOI + NSD exact, TopoScore proxy). The
 * official score showed the entire headroom is topology, so this sweeps the knobs
 * that crush b1 (dilation thickness x median iterations) and reports the best.
 *
 * Usage: surface_opt IMAGE.tif LABEL.tif [s_min=0.3] [i_min=80]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
#include "eval/topo.h"

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr,"usage: %s IMAGE LABEL [s_min] [i_min]\n",argv[0]); return 2; }
  float s_min=argc>3?(float)atof(argv[3]):0.3f, i_min=argc>4?(float)atof(argv[4]):80.f;
  int nz,ny,nx,lz,ly,lx;
  u8 *img=tiff_load_u8(argv[1],&nz,&ny,&nx), *lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim error\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  // GT topology for reference
  u8 *gtm=malloc(n); for(size_t i=0;i<n;i++) gtm[i]= lab[i]==1;
  topo_betti gb=betti_numbers_6(gtm,nz,ny,nx);
  printf("GT surface: b0=%ld b1=%ld b2=%ld\n", gb.b0,gb.b1,gb.b2); free(gtm);

  f32 *volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet=malloc(n*sizeof(f32)), *normal=malloc(3*n*sizeof(f32));
  fprintf(stderr,"structure tensor...\n");
  st_params p={0.5f,1.0f,1.0f};
  st_sheet_detect(volf,nz,ny,nx,&p,sheet,normal);
  u8 *ridge=malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,1.0f,ridge);
  remove_small_components(ridge,nz,ny,nx,TOPO_CONN26,100);
  free(volf);free(sheet);free(normal);free(img);

  int dils[]={1,2,3}, meds[]={2,4,6,8,12};
  u8 *q=malloc(n);
  printf("\n dil med  predV    Score   Topo   SurfD   VOI    b1\n");
  printf("------------------------------------------------------------\n");
  double best=-1; int bd=0,bm=0;
  for(size_t a=0;a<sizeof(dils)/sizeof(*dils);a++)
   for(size_t b=0;b<sizeof(meds)/sizeof(*meds);b++){
     memcpy(q,ridge,n);
     dilate_ball(q,q,nz,ny,nx,dils[a]);
     majority_filter(q,q,nz,ny,nx,meds[b]);
     fill_holes(q,nz,ny,nx);
     comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
     printf("  %d  %2d  b0=%-4ld  %.4f  %.4f  %.4f  %.4f  %5ld\n",
            dils[a],meds[b], s.pred_b0, s.score, s.topo_score, s.surface_dice, s.voi_score, s.pred_b1);
     fflush(stdout);
     if(s.score>best){ best=s.score; bd=dils[a]; bm=meds[b]; }
   }
  printf("\nbest proxy score %.4f at dilate=%d median=%d\n", best, bd, bm);
  free(ridge);free(q);free(lab);
  return 0;
}
