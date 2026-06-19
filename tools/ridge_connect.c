/* ridge_connect.c — reconnect the fragmented ridge with normal-aware (in-plane)
 * closing and compare against the blunt isotropic dilate+median. Computes the
 * structure tensor + ridge ONCE (keeping the normal field), then evaluates several
 * connection strategies. Goal: cut b0/b1 (fragmentation/porosity) without the
 * across-sheet thickening that caps Dice.
 *
 * Usage: ridge_connect IMAGE.tif LABEL.tif [s_min=0.3] [i_min=80] [step=1] [tol=2]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/metrics.h"

static f32 *g_normal;
static int G_NZ, G_NY, G_NX;

static void score(const char *name, u8 *p, const u8 *lab, int tol) {
  surface_eval e = eval_surface(p, lab, G_NZ, G_NY, G_NX, tol, 1, 2);
  printf("%-26s predV=%7zu  Dice=%.4f  SurfD=%.4f  b0=%4ld  b1=%5ld\n",
         name, e.pred_vox, e.dice, e.surface_dice, e.pred_b0, e.pred_b1);
}

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr,"usage: %s IMAGE LABEL [s_min] [i_min] [step] [tol]\n",argv[0]); return 2; }
  float s_min=argc>3?(float)atof(argv[3]):0.3f, i_min=argc>4?(float)atof(argv[4]):80.f;
  float step =argc>5?(float)atof(argv[5]):1.0f; int tol=argc>6?atoi(argv[6]):2;

  int nz,ny,nx,lz,ly,lx;
  u8 *img=tiff_load_u8(argv[1],&nz,&ny,&nx), *lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz||ly!=ny||lx!=nx){fprintf(stderr,"load/dim error\n");return 1;}
  G_NZ=nz;G_NY=ny;G_NX=nx; size_t n=(size_t)nz*ny*nx;
  f32 *volf=(f32*)malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i];
  f32 *sheet=(f32*)malloc(n*sizeof(f32)); g_normal=(f32*)malloc(3*n*sizeof(f32));
  fprintf(stderr,"computing structure tensor...\n");
  st_params p={0.5f, 1.0f, 1.0f};
  st_sheet_detect(volf,nz,ny,nx,&p,sheet,g_normal);
  u8 *ridge=(u8*)malloc(n);
  ridge_nms(volf,sheet,g_normal,nz,ny,nx,s_min,i_min,step,ridge);
  remove_small_components(ridge,nz,ny,nx,TOPO_CONN26,100);
  free(volf);free(sheet);free(img);

  u8 *q=(u8*)malloc(n);
  printf("(gt surface b0=5)\n");
  memcpy(q,ridge,n); fill_holes(q,nz,ny,nx); score("baseline (ridge+fill)", q, lab, tol);

  for (int it=1; it<=3; it++) {
    char nm[64];
    memcpy(q,ridge,n); inplane_close(q,g_normal,nz,ny,nx,0.5f,it); fill_holes(q,nz,ny,nx);
    snprintf(nm,sizeof nm,"inplane_close(0.5, x%d)", it); score(nm, q, lab, tol);
  }
  memcpy(q,ridge,n); inplane_close(q,g_normal,nz,ny,nx,0.7f,2); fill_holes(q,nz,ny,nx);
  score("inplane_close(0.7, x2)", q, lab, tol);

  memcpy(q,ridge,n); dilate_ball(q,q,nz,ny,nx,1); majority_filter(q,q,nz,ny,nx,2); fill_holes(q,nz,ny,nx);
  score("isotropic dilate1+med2", q, lab, tol);

  free(ridge);free(q);free(g_normal);free(lab);
  return 0;
}
