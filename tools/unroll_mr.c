/* unroll_mr.c — MULTI-RESOLUTION scroll unwrapping at native (LOD0) resolution.
 *
 * The full scroll is too big for LOD0 in RAM, so split geometry from content:
 *   (1) COARSE pass: read the whole scroll at `clod` (fits in RAM), build the
 *       global winding field W_coarse (every point's continuous wrap coordinate).
 *   (2) NATIVE pass: read a LOD0 region; place each material voxel's FULL-RES
 *       intensity into the unrolled image at its global winding coord, sampled
 *       (trilinear) from W_coarse. The coarse field gives a coordinate system
 *       shared by all LOD0 tiles, so they stitch into one consistent unroll while
 *       the detail stays native.
 *
 * Usage: unroll_mr ARCHIVE.mca OUT.tif clod z0 y0 x0 d [ppw=32]
 *        (z0,y0,x0,d are LOD0 voxel coords of the native region to unroll)
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

static f32 trilin(const f32 *w, int nz,int ny,int nx, double z,double y,double x){
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0; size_t nynx=(size_t)ny*nx;
  #define G(a,b,c) w[((size_t)(a)*ny+(b))*nx+(c)]
  double c00=G(z0,y0,x0)*(1-dx)+G(z0,y0,x1)*dx, c01=G(z0,y1,x0)*(1-dx)+G(z0,y1,x1)*dx;
  double c10=G(z1,y0,x0)*(1-dx)+G(z1,y0,x1)*dx, c11=G(z1,y1,x0)*(1-dx)+G(z1,y1,x1)*dx;
  #undef G
  return (f32)((c00*(1-dy)+c01*dy)*(1-dz)+(c10*(1-dy)+c11*dy)*dz);
}

int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: %s ARCHIVE.mca OUT.tif clod z0 y0 x0 d [ppw=32]\n",argv[0]);return 2;}
  const char*path=argv[1],*outp=argv[2];
  int clod=atoi(argv[3]); long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]); int d=atoi(argv[7]);
  int ppw=argc>8?atoi(argv[8]):32;
  int fz,fy,fx,nl; float qual;
  if(mca_dims(path,&fz,&fy,&fx,&qual,&nl)){fprintf(stderr,"open failed\n");return 1;}
  double s=(double)(1<<clod);
  int cz=(int)(fz/s),cy=(int)(fy/s),cx=(int)(fx/s);
  printf("scroll %dx%dx%d; coarse LOD%d = %dx%dx%d\n",fz,fy,fx,clod,cz,cy,cx);

  // ---- (1) coarse global winding ----
  fprintf(stderr,"reading whole scroll at LOD%d...\n",clod);
  u8 *cv=mca_load_region(path,clod,0,0,0,cz,cy,cx);
  if(!cv){fprintf(stderr,"coarse read failed\n");return 1;}
  size_t cn=(size_t)cz*cy*cx;
  u8 *cm=malloc(cn); for(size_t i=0;i<cn;i++) cm[i]=cv[i]>=80; majority_filter(cm,cm,cz,cy,cx,1);
  remove_small_components(cm,cz,cy,cx,TOPO_CONN26,50);
  umbilicus umb; umb.n=2; umb.z=malloc(8);umb.y=malloc(8);umb.x=malloc(8);
  umb.z[0]=0;umb.z[1]=(f32)(cz-1); umb.y[0]=umb.y[1]=(f32)(cy*0.5); umb.x[0]=umb.x[1]=(f32)(cx*0.5);
  // pitch estimate (coarse voxels) along rays at mid-z
  double gaps[100000]; int ng=0; int zmid=cz/2,maxr=(int)(0.95*cy*0.5);
  for(int a=0;a<256&&ng<90000;a++){double an=a*2*M_PI/256,ca=cos(an),sa=sin(an);int ins=0;double last=-1;
    for(double r=2;r<maxr;r+=0.5){int yy=(int)(cy*0.5+r*sa),xx=(int)(cx*0.5+r*ca); if(yy<0||yy>=cy||xx<0||xx>=cx)break;
      int on=cm[((size_t)zmid*cy+yy)*cx+xx]!=0; if(on&&!ins){if(last>0&&r-last>1.5)gaps[ng++]=r-last;last=r;} ins=on;}}
  double pitch=8; if(ng>20){for(int i=0;i<ng;i++)for(int j=i+1;j<ng;j++)if(gaps[j]<gaps[i]){double t=gaps[i];gaps[i]=gaps[j];gaps[j]=t;}pitch=gaps[ng/2];}
  printf("coarse pitch %.1f vox (%d crossings)\n",pitch,ng);
  fprintf(stderr,"solving coarse winding...\n");
  f32 *cw=malloc(cn*sizeof(f32)); wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)pitch;
  if(winding_field_solve(cm,cz,cy,cx,&umb,&wp,NULL,NULL,cw)){fprintf(stderr,"winding failed\n");return 1;}
  free(cv);free(cm);

  // ---- (2) native LOD0 region, placed at global winding coord ----
  fprintf(stderr,"reading native LOD0 region %d^3 @ (%ld,%ld,%ld)...\n",d,z0,y0,x0);
  u8 *img=mca_load_region(path,0,(int)z0,(int)y0,(int)x0,d,d,d);
  if(!img){fprintf(stderr,"native read failed\n");return 1;}
  // winding range over the region (sample coarse field at region corners->voxels)
  float wmin=1e30f,wmax=-1e30f;
  for(int z=0;z<d;z+=8)for(int y=0;y<d;y+=8)for(int x=0;x<d;x+=8){
    if(img[((size_t)z*d+y)*d+x]<80) continue;
    f32 w=trilin(cw,cz,cy,cx,(z0+z)/s,(y0+y)/s,(x0+x)/s);
    if(w<wmin)wmin=w; if(w>wmax)wmax=w; }
  printf("region global winding %.2f .. %.2f (~%.1f wraps)\n",wmin,wmax,wmax-wmin);
  int W=(int)((wmax-wmin)*ppw)+1; if(W<1)W=1; if(W>60000)W=60000; int H=d;
  double *acc=calloc((size_t)W*H,sizeof(double)); int *cnt=calloc((size_t)W*H,sizeof(int));
  #pragma omp parallel for schedule(static)
  for(int z=0;z<d;z++)for(int y=0;y<d;y++)for(int x=0;x<d;x++){
    size_t i=((size_t)z*d+y)*d+x; if(img[i]<80) continue;
    f32 w=trilin(cw,cz,cy,cx,(z0+z)/s,(y0+y)/s,(x0+x)/s);
    int u=(int)((w-wmin)*ppw); if(u<0||u>=W) continue;
    size_t p=(size_t)z*W+u;
    #pragma omp atomic
    acc[p]+=img[i];
    #pragma omp atomic
    cnt[p]++;
  }
  u8 *flat=calloc((size_t)W*H,1);
  for(size_t i=0;i<(size_t)W*H;i++) if(cnt[i]) flat[i]=(u8)(acc[i]/cnt[i]);
  tiff_save_u8(outp,flat,1,H,W);
  size_t fill=0; for(size_t i=0;i<(size_t)W*H;i++) fill+=flat[i]>0;
  printf("wrote NATIVE-res unroll %s (%dx%d), %.0f%% filled\n",outp,W,H,100.0*fill/((size_t)W*H));
  free(cw);free(img);free(acc);free(cnt);free(flat);free(umb.z);free(umb.y);free(umb.x);
  return 0;
}
