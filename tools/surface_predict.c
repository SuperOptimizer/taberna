/* surface_predict.c — produce a surface-prediction TIFF and print the proxy
 * composite score. The TIFF can be fed to the OFFICIAL metric via
 * scripts/official_score.py for the authoritative number.
 *
 * Pipeline: CT -> structure-tensor sheetness + normal -> ridge NMS -> dust ->
 * dilate(r) -> median(iters) -> fill_holes. Writes pred as a multi-page u8 TIFF
 * (0/1) and scores it with the in-process proxy of the Kaggle composite.
 *
 * Usage: surface_predict IMAGE.tif LABEL.tif OUT_PRED.tif
 *        [s_min=0.3] [i_min=80] [dilate=1] [median=2] [tol=2]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"

int main(int argc, char **argv) {
  if (argc < 4) {
    fprintf(stderr,"usage: %s IMAGE.tif LABEL.tif OUT_PRED.tif "
                   "[s_min=0.3] [i_min=80] [dilate=1] [median=2] [tol=2]\n", argv[0]);
    return 2;
  }
  float s_min=argc>4?(float)atof(argv[4]):0.3f, i_min=argc>5?(float)atof(argv[5]):80.f;
  int dil=argc>6?atoi(argv[6]):1, med=argc>7?atoi(argv[7]):2, tol=argc>8?atoi(argv[8]):2;

  int nz,ny,nx,lz,ly,lx;
  u8 *img=tiff_load_u8(argv[1],&nz,&ny,&nx), *lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim error\n");return 1;}
  size_t n=(size_t)nz*ny*nx;

  f32 *volf=(f32*)malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet=(f32*)malloc(n*sizeof(f32)), *normal=(f32*)malloc(3*n*sizeof(f32));
  fprintf(stderr,"structure tensor...\n");
  st_params p={0.5f,1.0f,1.0f};
  st_sheet_detect(volf,nz,ny,nx,&p,sheet,normal);
  u8 *pred=(u8*)malloc(n);
  ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,1.0f,pred);
  remove_small_components(pred,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0) dilate_ball(pred,pred,nz,ny,nx,dil);
  if(med>0) majority_filter(pred,pred,nz,ny,nx,med);
  fill_holes(pred,nz,ny,nx);
  free(volf);free(sheet);free(normal);free(img);

  if(tiff_save_u8(argv[3],pred,nz,ny,nx)!=0){fprintf(stderr,"write %s failed\n",argv[3]);return 1;}
  fprintf(stderr,"wrote %s\n", argv[3]);

  comp_score s=competition_score(pred,lab,nz,ny,nx,tol,1,2);
  printf("=== composite (TopoScore is a proxy; official_score.py for exact) ===\n");
  printf("Score        : %.4f\n", s.score);
  printf("  TopoScore  : %.4f  (PROXY: clamped incl-excl; real = Betti-matching)\n", s.topo_score);
  printf("  SurfaceDice: %.4f  (EXACT: native Google NSD port)\n", s.surface_dice);
  printf("  VOI_score  : %.4f  (EXACT: 1/(1+0.3*VOI), VOI=%.4f)\n", s.voi_score, s.voi);
  printf("  pred betti : b0=%ld b1=%ld b2=%ld\n", s.pred_b0, s.pred_b1, s.pred_b2);

  free(pred);free(lab);
  return 0;
}
