/* ridge_clean.c — given the best ridge-detector params, sweep the *cleanup* knobs.
 * The NMS ridge is ~1 voxel thick and porous (huge b1). Giving it thickness lets
 * the iterated median filter close the porosity without eroding the sheet away.
 * Computes structure tensor + ridge ONCE, then sweeps (dilate radius, median iters).
 *
 * Usage: ridge_clean IMAGE.tif LABEL.tif [s_min=0.3] [i_min=80] [step=1] [tol=2]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/metrics.h"

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr,"usage: %s IMAGE LABEL [s_min] [i_min] [step] [tol]\n",argv[0]); return 2; }
  float s_min = argc>3?(float)atof(argv[3]):0.3f;
  float i_min = argc>4?(float)atof(argv[4]):80.f;
  float step  = argc>5?(float)atof(argv[5]):1.0f;
  int   tol   = argc>6?atoi(argv[6]):2;

  int nz,ny,nx,lz,ly,lx;
  u8 *img=tiff_load_u8(argv[1],&nz,&ny,&nx), *lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim error\n");return 1;}
  size_t n=(size_t)nz*ny*nx;

  f32 *volf=(f32*)malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet=(f32*)malloc(n*sizeof(f32)), *normal=(f32*)malloc(3*n*sizeof(f32));
  fprintf(stderr,"computing structure tensor...\n");
  st_sheet_detect(volf,nz,ny,nx,NULL,sheet,normal);
  u8 *ridge=(u8*)malloc(n);
  size_t kept = ridge_nms(volf,sheet,normal,nz,ny,nx,s_min,i_min,step,ridge);
  remove_small_components(ridge,nz,ny,nx,TOPO_CONN26,100);
  printf("ridge s_min=%.2f i_min=%.0f -> %zu voxels (after dust)\n\n", s_min,i_min,kept);
  free(sheet);free(normal);free(volf);free(img);

  u8 *p=(u8*)malloc(n);
  printf(" dilate  median   predV    Dice   SurfDice   b1\n");
  printf("------------------------------------------------------\n");
  int rs[]={0,1}, its[]={0,1,2,3};
  for(size_t a=0;a<sizeof(rs)/sizeof(*rs);a++)
    for(size_t b=0;b<sizeof(its)/sizeof(*its);b++){
      memcpy(p,ridge,n);
      if(rs[a]>0) dilate_ball(p,p,nz,ny,nx,rs[a]);
      if(its[b]>0) majority_filter(p,p,nz,ny,nx,its[b]);
      fill_holes(p,nz,ny,nx);
      surface_eval e=eval_surface(p,lab,nz,ny,nx,tol,1,2);
      printf("   %d      %d     %7zu  %.4f  %.4f   %5ld\n",
             rs[a],its[b],e.pred_vox,e.dice,e.surface_dice,e.pred_b1);
    }
  free(ridge);free(p);free(lab);
  return 0;
}
