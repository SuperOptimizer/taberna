/* st_sweep.c — sweep structure-tensor scales to sharpen the ridge normal and
 * improve placement precision. The ridge lands a systematic 1-2 voxels off the
 * label (Surface Dice jumps with tolerance), which points at the normal estimate;
 * the integration scale sigma_tensor must stay below the inter-wrap gap or it
 * averages adjacent-wrap normals together. Recomputes the tensor per (sigma_grad,
 * sigma_tensor) and scores ridge+dust+fill at tol 0/1/2.
 *
 * Usage: st_sweep IMAGE.tif LABEL.tif [s_min=0.3] [i_min=80] [step=1]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/metrics.h"

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr,"usage: %s IMAGE LABEL [s_min] [i_min] [step]\n",argv[0]); return 2; }
  float s_min=argc>3?(float)atof(argv[3]):0.3f, i_min=argc>4?(float)atof(argv[4]):80.f;
  float step =argc>5?(float)atof(argv[5]):1.0f;

  int nz,ny,nx,lz,ly,lx;
  u8 *img=tiff_load_u8(argv[1],&nz,&ny,&nx), *lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim error\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32 *volf=(f32*)malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet=(f32*)malloc(n*sizeof(f32)), *normal=(f32*)malloc(3*n*sizeof(f32));
  u8 *pred=(u8*)malloc(n);
  free(img);

  float sg_grid[]={0.5f,1.0f};
  float st_grid[]={0.5f,0.75f,1.0f,1.5f,2.0f};
  printf("sigma_grad sigma_tens  predV   SurfD@0  SurfD@1  SurfD@2   b1\n");
  printf("---------------------------------------------------------------\n");
  for(size_t a=0;a<sizeof(sg_grid)/sizeof(*sg_grid);a++)
    for(size_t b=0;b<sizeof(st_grid)/sizeof(*st_grid);b++){
      st_params p={sg_grid[a], st_grid[b], 1.0f};
      st_sheet_detect(volf,nz,ny,nx,&p,sheet,normal);
      ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,step,pred);
      remove_small_components(pred,nz,ny,nx,TOPO_CONN26,100);
      fill_holes(pred,nz,ny,nx);
      surface_eval e0=eval_surface(pred,lab,nz,ny,nx,0,1,2);
      surface_eval e1=eval_surface(pred,lab,nz,ny,nx,1,1,2);
      surface_eval e2=eval_surface(pred,lab,nz,ny,nx,2,1,2);
      printf("  %.2f      %.2f    %7zu  %.4f   %.4f   %.4f  %5ld\n",
             sg_grid[a], st_grid[b], e2.pred_vox,
             e0.surface_dice, e1.surface_dice, e2.surface_dice, e2.pred_b1);
      fflush(stdout);
    }
  free(volf);free(sheet);free(normal);free(pred);free(lab);
  return 0;
}
