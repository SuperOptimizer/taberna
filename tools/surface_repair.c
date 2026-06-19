/* surface_repair.c — full classical pipeline with watertight sheet repair:
 * CT -> ridge -> consolidate (dust/dilate/median) -> PCA height-map repair ->
 * write pred TIFF + native composite. The repair makes each sheet a height
 * surface (b1=0), attacking the topology wall the metric punishes.
 *
 * Usage: surface_repair IMAGE LABEL OUT.tif [s_min=0.3] [i_min=80]
 *                       [dilate=2] [median=4] [min_comp=500]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "postproc/sheet_repair.h"
#include "eval/score.h"

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s IMAGE LABEL OUT.tif [s_min][i_min][dilate][median][min_comp]\n",argv[0]);return 2;}
  float s_min=argc>4?(float)atof(argv[4]):0.3f, i_min=argc>5?(float)atof(argv[5]):80.f;
  int dil=argc>6?atoi(argv[6]):2, med=argc>7?atoi(argv[7]):4, minc=argc>8?atoi(argv[8]):500;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32*sheet=malloc(n*sizeof(f32)),*normal=malloc(3*n*sizeof(f32));
  fprintf(stderr,"structure tensor...\n");
  st_params p={0.5f,1.0f,1.0f};
  st_sheet_detect(volf,nz,ny,nx,&p,sheet,normal);
  u8*q=malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,1.0f,q);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0) dilate_ball(q,q,nz,ny,nx,dil);
  if(med>0) majority_filter(q,q,nz,ny,nx,med);
  free(volf);free(sheet);free(normal);free(img);

  fprintf(stderr,"sheet repair...\n");
  u8*rep=malloc(n);
  sheet_repair(q,rep,nz,ny,nx,minc);

  if(tiff_save_u8(argv[3],rep,nz,ny,nx)!=0){fprintf(stderr,"write fail\n");return 1;}
  comp_score s=competition_score(rep,lab,nz,ny,nx,2,1,2);
  printf("=== repaired (native; official_score.py for exact topo) ===\n");
  printf("Score %.4f  Topo %.4f  SurfD %.4f(exact)  VOI %.4f(exact)  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(q);free(rep);free(lab);
  return 0;
}
