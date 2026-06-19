/* surface_prep.c — preprocess the CT (fysics) before detection.
 * mode: none | denoise | musica | both. Then structure-tensor ridge + score.
 * Tests whether cleaning/enhancing the data (faint sheets are the porosity cause)
 * lifts the score. Usage: surface_prep IMAGE LABEL MODE [dil=2] [med=4] */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/tiff_vol.h"
#include "fysics.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/score.h"
int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s IMAGE LABEL MODE[none|denoise|musica|both] [dil][med]\n",argv[0]);return 2;}
  const char*mode=argv[3]; int dil=argc>4?atoi(argv[4]):2,med=argc>5?atoi(argv[5]):4;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*v=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) v[i]=(f32)img[i]; free(img);

  int den = (strstr(mode,"denoise")||!strcmp(mode,"both"));
  int mus = (strstr(mode,"musica")||!strcmp(mode,"both"));
  if(den){
    fy_noise_model nm; fy_estimate_noise(v,nz,ny,nx,5,10.0,100.0,&nm);
    double eps=fy_guided_eps_for_noise(nm.noise_ref);
    fprintf(stderr,"noise_ref=%.3f eps=%.4f\n",nm.noise_ref,eps);
    f32*o=malloc(n*sizeof(f32)); fy_guided_denoise(v,o,nz,ny,nx,2,eps); memcpy(v,o,n*sizeof(f32)); free(o);
  }
  if(mus){
    fprintf(stderr,"musica per-slice...\n");
    f32*o=malloc((size_t)ny*nx*sizeof(f32));
    for(int z=0;z<nz;z++){ fy_musica2d(v+(size_t)z*ny*nx,o,ny,nx,4,0.6f,0.1f);
      memcpy(v+(size_t)z*ny*nx,o,(size_t)ny*nx*sizeof(f32)); }
    free(o);
  }
  f32*sheet=malloc(n*sizeof(f32)),*normal=malloc(3*n*sizeof(f32));
  st_params sp={0.5f,1.0f,1.0f};
  fprintf(stderr,"detect...\n"); st_sheet_detect(v,nz,ny,nx,&sp,sheet,normal);
  u8*q=malloc(n);
  ridge_nms(v,sheet,normal,nz,ny,nx,0.3f,80.f,1.0f,q);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,100);
  if(dil>0)dilate_ball(q,q,nz,ny,nx,dil);
  if(med>0)majority_filter(q,q,nz,ny,nx,med);
  fill_holes(q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("PREP[%s]: Score %.4f  SurfD %.4f  VOI %.4f  b0=%ld b1=%ld\n",
         mode,s.score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(v);free(sheet);free(normal);free(q);free(lab);return 0;
}
