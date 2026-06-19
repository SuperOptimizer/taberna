/* surface_ced.c — CED front-end: connect the sheet field before detection.
 * CT -> coherence-enhancing diffusion -> ridge -> consolidate -> score + write.
 * Tests whether reducing porosity at the source lowers b1 / lifts the score.
 *
 * Usage: surface_ced IMAGE LABEL OUT.tif [iters=12] [s_min=0.3] [i_min=80]
 *                    [dilate=2] [median=4]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ced.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s IMAGE LABEL OUT.tif [iters][s_min][i_min][dilate][median]\n",argv[0]);return 2;}
  int iters=argc>4?atoi(argv[4]):12;
  float s_min=argc>5?(float)atof(argv[5]):0.3f, i_min=argc>6?(float)atof(argv[6]):80.f;
  int dil=argc>7?atoi(argv[7]):2, med=argc>8?atoi(argv[8]):4;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i]; free(img);

  fprintf(stderr,"CED (%d iters)...\n",iters);
  ced_params cp=ced_default_params(); cp.iters=iters;
  ced_diffuse(volf,nz,ny,nx,&cp);

  fprintf(stderr,"ridge on diffused...\n");
  f32*sheet=malloc(n*sizeof(f32)),*normal=malloc(3*n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f};
  st_sheet_detect(volf,nz,ny,nx,&sp,sheet,normal);
  u8*q=malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,1.0f,q);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0) dilate_ball(q,q,nz,ny,nx,dil);
  if(med>0) majority_filter(q,q,nz,ny,nx,med);
  fill_holes(q,nz,ny,nx);
  free(volf);free(sheet);free(normal);

  tiff_save_u8(argv[3],q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("=== CED pipeline (native; official_score.py for exact topo) ===\n");
  printf("Score %.4f  Topo %.4f  SurfD %.4f  VOI %.4f  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(q);free(lab);
  return 0;
}
