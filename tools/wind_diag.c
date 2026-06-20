/* wind_diag.c — diagnose the winding field on a CROSS-SECTION slice. Solves the
 * coarse winding (same as unroll_full) and dumps a mid-z slice three ways so we can
 * SEE whether the iso-winding bands follow the actual papyrus sheets:
 *   OUT_mat.tif  — material intensity (the sheets)
 *   OUT_frac.tif — frac(winding)*255 (iso-W ramps; band edges should sit BETWEEN
 *                  sheets if the coordinate is correct)
 *   OUT_over.tif — overlay: material dimmed, with a white line drawn at each integer
 *                  winding crossing (where the wrap index ticks over)
 *
 * Usage: wind_diag ARCHIVE.mca OUTBASE clod zc [citers] [pitch]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "io/tiff_vol.h"
#include "postproc/morph.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"
#include "unwrap/winding_field.h"
#include "segmentation/sheet_tensor.h"

/* Data-driven pitch: median radial gap between consecutive sheet starts (raw
 * intensity) along rays from the umbilicus. = dr_per_winding (one sheet per wrap). */
static double measure_pitch(const u8*v,int nz,int ny,int nx,const umbilicus*umb,double fb){
  int z=nz/2; f32 cyf,cxf; umbilicus_center(umb,(f32)z,&cyf,&cxf);
  size_t nynx=(size_t)ny*nx; int maxr=(int)(0.9*(ny<nx?ny:nx)*0.5); if(maxr<12) return fb;
  int *hist=calloc(maxr+2,sizeof(int)); long total=0;
  double *p=malloc(maxr*sizeof(double)),*sm=malloc(maxr*sizeof(double));
  for(int a=0;a<360;a++){ double ca=cos(a*M_PI/180),sa=sin(a*M_PI/180); int rmax=0;
    for(int r=0;r<maxr;r++){ int yy=(int)(cyf+r*sa),xx=(int)(cxf+r*ca);
      if(yy<0||yy>=ny||xx<0||xx>=nx)break; p[r]=v[(size_t)z*nynx+(size_t)yy*nx+xx]; rmax=r; }
    if(rmax<12) continue;
    for(int r=0;r<=rmax;r++){ double s=0;int c=0; for(int d=-2;d<=2;d++){int rr=r+d; if(rr>=0&&rr<=rmax){s+=p[rr];c++;}} sm[r]=s/c; }
    int last=-1;
    for(int r=3;r<=rmax-3;r++){
      if(sm[r]>sm[r-1]&&sm[r]>=sm[r+1]&&sm[r]>sm[r-3]+2&&sm[r]>sm[r+3]+2){   // prominent peak
        if(last>0){ int g=r-last; if(g>=2&&g<maxr){hist[g]++;total++;} } last=r; } }
  }
  free(p);free(sm);
  if(total<50){ free(hist); return fb; }
  long acc=0; int med=(int)fb; for(int g=2;g<maxr;g++){ acc+=hist[g]; if(acc*2>=total){med=g;break;} }
  free(hist); return med>=2?med:fb;
}

/* Build the Poisson forcing div(g), g = (sheetness * oriented sheet-normal)/pitch,
 * so the winding gradient follows the ACTUAL sheets. Normals are oriented outward
 * (away from the umbilicus axis = increasing-winding direction). Returns malloc'd
 * forcing[nz*ny*nx], or NULL on failure. */
static f32 *normal_forcing(const u8 *cv, const u8 *mask, int nz,int ny,int nx,
                           const umbilicus *umb, double pitch){
  size_t n=(size_t)nz*ny*nx, nynx=(size_t)ny*nx;
  f32 *vol=malloc(n*sizeof(f32)); f32 *nrm=malloc(3*n*sizeof(f32)); f32 *sh=malloc(n*sizeof(f32));
  if(!vol||!nrm||!sh){free(vol);free(nrm);free(sh);return NULL;}
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<n;i++) vol[i]=cv[i];
  st_params sp=st_default_params(); sp.sigma_tensor=1.5f; // keep below inter-wrap gap
  st_sheet_detect(vol,nz,ny,nx,&sp,sh,nrm); free(vol);
  // orient outward and pre-scale by sheetness: g = s * n_oriented / pitch (stored in nrm)
  #pragma omp parallel for schedule(static)
  for(int z=0;z<nz;z++){ f32 cyf,cxf; umbilicus_center(umb,(f32)z,&cyf,&cxf);
    for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
      size_t i=(size_t)z*nynx+(size_t)y*nx+x;
      f32 nxv=nrm[3*i+0],nyv=nrm[3*i+1],nzv=nrm[3*i+2];
      double rx=x-cxf, ry=y-cyf; double dot=nxv*rx+nyv*ry;     // axis ~ vertical, ignore z
      double sgn=(dot<0)?-1.0:1.0; double w=(double)sh[i]/pitch*sgn;
      nrm[3*i+0]=(f32)(nxv*w); nrm[3*i+1]=(f32)(nyv*w); nrm[3*i+2]=(f32)(nzv*w); }}
  free(sh);
  // forcing[i] = div(g) via masked central differences (interleaved nx,ny,nz)
  f32 *forcing=calloc(n,sizeof(f32));
  #pragma omp parallel for schedule(static)
  for(int z=0;z<nz;z++)for(int y=0;y<ny;y++)for(int x=0;x<nx;x++){
    size_t i=(size_t)z*nynx+(size_t)y*nx+x; if(!mask[i]) continue;
    double div=0;
    if(x>0&&x+1<nx&&mask[i-1]&&mask[i+1])      div+=(nrm[3*(i+1)+0]-nrm[3*(i-1)+0])*0.5;
    if(y>0&&y+1<ny&&mask[i-nx]&&mask[i+nx])    div+=(nrm[3*(i+nx)+1]-nrm[3*(i-nx)+1])*0.5;
    if(z>0&&z+1<nz&&mask[i-nynx]&&mask[i+nynx])div+=(nrm[3*(i+nynx)+2]-nrm[3*(i-nynx)+2])*0.5;
    forcing[i]=(f32)div; }
  free(nrm);
  return forcing;
}

int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: %s ARCHIVE OUTBASE clod zc [citers] [pitch] [usenormal]\n",argv[0]);return 2;}
  const char*path=argv[1],*base=argv[2]; int clod=atoi(argv[3]); long zc=atol(argv[4]);
  int citers=argc>5?atoi(argv[5]):-1; double pitch=argc>6?atof(argv[6]):-1;
  int usenormal=argc>7?atoi(argv[7]):0; double lambda=argc>8?atof(argv[8]):0.05;
  int warp=argc>9?atoi(argv[9]):0;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  double s=(double)(1<<clod); int cz=(int)(fz/s),cy=(int)(fy/s),cx=(int)(fx/s);
  fprintf(stderr,"coarse LOD%d %dx%dx%d\n",clod,cz,cy,cx);
  u8*cv=mca_read(arc,clod,0,0,0,cz,cy,cx); if(!cv){fprintf(stderr,"read fail\n");return 1;}
  size_t cn=(size_t)cz*cy*cx; u8*cm=malloc(cn);
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<cn;i++) cm[i]=cv[i]!=0;
  
  umbilicus umb; if(umbilicus_estimate(cm,cz,cy,cx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;}
  f32*cw=malloc(cn*sizeof(f32)); wfield_params wp=winding_default_params();
  wp.dr_per_winding=(f32)(pitch>0?pitch:measure_pitch(cv,cz,cy,cx,&umb,96.0/s));
  fprintf(stderr,"pitch=%.2f vox\n",wp.dr_per_winding); if(citers>=0) wp.iters=citers;
  f32 *forcing=NULL;
  if(usenormal){ fprintf(stderr,"computing sheet-normal forcing (lambda=%.3f)...\n",lambda);
    forcing=normal_forcing(cv,cm,cz,cy,cx,&umb,wp.dr_per_winding); wp.forcing=forcing;
    wp.anchor_lambda=(f32)lambda; }
  if(warp){ fprintf(stderr,"contour-warp init (deform to scroll shape)...\n");
    winding_contour_warp(cm,cz,cy,cx,&umb,wp.dr_per_winding,cw); wp.warm_start=1; }
  winding_field_solve(cm,cz,cy,cx,&umb,&wp,NULL,NULL,cw);
  free(forcing);

  int z=(int)(zc/s); if(z<0)z=0; if(z>=cz)z=cz-1;
  size_t off=(size_t)z*cy*cx;
  u8*mat=malloc((size_t)cy*cx),*frac=malloc((size_t)cy*cx),*over=malloc((size_t)cy*cx);
  for(int y=0;y<cy;y++)for(int x=0;x<cx;x++){
    size_t i=(size_t)y*cx+x; mat[i]=cv[off+i];
    double w=cw[off+i]; double f=w-floor(w);
    frac[i]= cm[off+i]? (u8)(f*255):0;
    // overlay: dim material; draw a white tick where an integer winding crosses
    // between this pixel and its +x neighbor (the wrap index changes there).
    int line=0;
    if(x+1<cx && cm[off+i] && cm[off+i+1]){ double w2=cw[off+i+1]; if(floor(w)!=floor(w2)) line=1; }
    over[i]= line?255:(u8)(mat[i]/3);
  }
  char p[512];
  snprintf(p,sizeof p,"%s_mat.tif",base);  tiff_save_u8(p,mat,1,cy,cx);
  snprintf(p,sizeof p,"%s_frac.tif",base); tiff_save_u8(p,frac,1,cy,cx);
  snprintf(p,sizeof p,"%s_over.tif",base); tiff_save_u8(p,over,1,cy,cx);
  float wmin=1e30f,wmax=-1e30f; for(size_t i=0;i<cn;i++) if(cm[i]){if(cw[i]<wmin)wmin=cw[i];if(cw[i]>wmax)wmax=cw[i];}
  // Objective bunching metric: along rays from the umbilicus, the winding should
  // increase ~linearly, so consecutive integer-winding crossings are evenly spaced.
  // Report the coefficient of variation (std/mean) of inter-crossing gaps: lower =
  // more uniform wraps = less bunching/crossing in the unroll.
  f32 ccy,ccx; umbilicus_center(&umb,(f32)z,&ccy,&ccx);
  double sg=0,sg2=0; long ng=0; int maxr=(int)(0.95*(cy<cx?cy:cx)*0.5);
  for(int a=0;a<720;a++){ double an=a*M_PI/360.0,ca=cos(an),sa=sin(an);
    double lastcross=-1, lastw=-1; int have=0;
    for(int r=2;r<maxr;r++){ int yy=(int)(ccy+r*sa),xx=(int)(ccx+r*ca);
      if(yy<0||yy>=cy||xx<0||xx>=cx) break; size_t i=off+(size_t)yy*cx+xx; if(!cm[i])continue;
      double w=cw[i]; if(have && floor(w)!=floor(lastw)){ if(lastcross>0){double g=r-lastcross; if(g>0.3&&g<3*wp.dr_per_winding){sg+=g;sg2+=g*g;ng++;}} lastcross=r; }
      else if(!have) lastcross=r; lastw=w; have=1; } }
  double mean=ng?sg/ng:0, var=ng?sg2/ng-mean*mean:0, cvar=mean>0?sqrt(var>0?var:0)/mean:0;
  printf("slice z=%d winding %.1f..%.1f (%.1f wraps) pitch=%.2f iters=%d normal=%d lambda=%.3f | gap mean=%.2f CV=%.3f (n=%ld)\n",
         z,wmin,wmax,wmax-wmin,wp.dr_per_winding,wp.iters,usenormal,usenormal?lambda:0.0,mean,cvar,ng);
  return 0;
}
