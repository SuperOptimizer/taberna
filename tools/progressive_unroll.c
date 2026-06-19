/* progressive_unroll.c — coarse-to-fine multigrid scroll unwrapping.
 *
 * Solve the winding field at the coarsest LOD (whole scroll, fast), then warm-start
 * each finer LOD from the upsampled coarser solution and refine. The global winding
 * geometry is carried down the pyramid as far as RAM allows (the finest GLOBAL
 * level), giving a consistent coordinate system; the native (LOD0) content is then
 * placed against it for a region. Reports the wrap count + quality (MVF) at each
 * level — the count should be stable across levels (it's physical), validating the
 * pitch estimate.
 *
 * Usage: progressive_unroll ARCHIVE.mca OUT.tif start_lod fine_lod
 *                           [z0 y0 x0 d ppw]   (LOD0 region to unroll natively)
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
#include "unwrap/wmetrics.h"

/* robust radial pitch: dominant period of the angle-averaged radial intensity
 * profile (autocorrelation), far less noisy than median crossing-gaps. */
static double robust_pitch(const u8 *v, int nz,int ny,int nx, double cy,double cx){
  int maxr=(int)(0.9*(ny<nx?ny:nx)*0.5); if(maxr<8) return 8;
  double *prof=calloc(maxr,sizeof(double)); int *pc=calloc(maxr,sizeof(int));
  int zmid=nz/2; size_t nynx=(size_t)ny*nx;
  for(int a=0;a<360;a++){ double an=a*M_PI/180.0,ca=cos(an),sa=sin(an);
    for(int r=0;r<maxr;r++){ int yy=(int)(cy+r*sa),xx=(int)(cx+r*ca);
      if(yy<0||yy>=ny||xx<0||xx>=nx)continue; prof[r]+=v[(size_t)zmid*nynx+(size_t)yy*nx+xx]; pc[r]++; }}
  for(int r=0;r<maxr;r++) if(pc[r]) prof[r]/=pc[r];
  double mean=0; for(int r=0;r<maxr;r++) mean+=prof[r]; mean/=maxr;
  for(int r=0;r<maxr;r++) prof[r]-=mean;            // de-mean for autocorr
  double bestc=-1e30; int bestL=8;
  for(int L=3;L<maxr/2;L++){ double c=0; for(int r=0;r+L<maxr;r++) c+=prof[r]*prof[r+L];
    if(c>bestc){bestc=c;bestL=L;} }
  free(prof);free(pc);
  return bestL;
}

static f32 trilin(const f32*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0;if(z>nz-1)z=nz-1;if(y<0)y=0;if(y>ny-1)y=ny-1;if(x<0)x=0;if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x,z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0; size_t nynx=(size_t)ny*nx;
  #define G(a,b,c) w[((size_t)(a)*ny+(b))*nx+(c)]
  double c00=G(z0,y0,x0)*(1-dx)+G(z0,y0,x1)*dx,c01=G(z0,y1,x0)*(1-dx)+G(z0,y1,x1)*dx;
  double c10=G(z1,y0,x0)*(1-dx)+G(z1,y0,x1)*dx,c11=G(z1,y1,x0)*(1-dx)+G(z1,y1,x1)*dx;
  #undef G
  return (f32)((c00*(1-dy)+c01*dy)*(1-dz)+(c10*(1-dy)+c11*dy)*dz);
}

int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: %s ARCHIVE.mca OUT.tif start_lod fine_lod [z0 y0 x0 d ppw]\n",argv[0]);return 2;}
  const char*path=argv[1],*outp=argv[2]; int slod=atoi(argv[3]),flod=atoi(argv[4]);
  int fz,fy,fx,nl; float ql; if(mca_dims(path,&fz,&fy,&fx,&ql,&nl)){fprintf(stderr,"open failed\n");return 1;}

  f32 *prevw=NULL; int pz=0,py=0,px=0;
  for(int lod=slod; lod>=flod; lod--){
    double s=(double)(1<<lod); int dz=(int)(fz/s),dy=(int)(fy/s),dx=(int)(fx/s);
    size_t nn=(size_t)dz*dy*dx;
    fprintf(stderr,"LOD%d: %dx%dx%d (%.2f GB vox)\n",lod,dz,dy,dx,nn/1e9);
    u8 *v=mca_load_region(path,lod,0,0,0,dz,dy,dx);
    if(!v){fprintf(stderr,"read failed\n");return 1;}
    u8 *m=malloc(nn); for(size_t i=0;i<nn;i++) m[i]=v[i]>=80; majority_filter(m,m,dz,dy,dx,1);
    remove_small_components(m,dz,dy,dx,TOPO_CONN26,50);
    double cyc=dy*0.5,cxc=dx*0.5;
    double pitch=robust_pitch(v,dz,dy,dx,cyc,cxc);
    umbilicus umb; umb.n=2; umb.z=malloc(8);umb.y=malloc(8);umb.x=malloc(8);
    umb.z[0]=0;umb.z[1]=(f32)(dz-1);umb.y[0]=umb.y[1]=(f32)cyc;umb.x[0]=umb.x[1]=(f32)cxc;
    f32 *w=malloc(nn*sizeof(f32));
    wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)pitch;
    if(prevw){ // warm-start: upsample prev winding into this grid
      #pragma omp parallel for schedule(static)
      for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){
        size_t i=((size_t)z*dy+y)*dx+x; w[i]= m[i]? trilin(prevw,pz,py,px,z*0.5,y*0.5,x*0.5):0; }
      wp.warm_start=1; wp.iters=20;
    }
    if(winding_field_solve(m,dz,dy,dx,&umb,&wp,NULL,NULL,w)){fprintf(stderr,"solve failed\n");return 1;}
    float wmin=1e30f,wmax=-1e30f; for(size_t i=0;i<nn;i++) if(m[i]){if(w[i]<wmin)wmin=w[i];if(w[i]>wmax)wmax=w[i];}
    double mvf=winding_mvf(m,w,dz,dy,dx,&umb,64,dz/2);
    printf("LOD%d  pitch=%.1f vox  wraps=%.1f (%.2f..%.2f)  MVF=%.4f\n",lod,pitch,wmax-wmin,wmin,wmax,mvf);
    fflush(stdout);
    free(v);free(m);free(umb.z);free(umb.y);free(umb.x);
    free(prevw); prevw=w; pz=dz;py=dy;px=dx;
  }
  printf("\nfinest global winding at LOD%d (%dx%dx%d)\n",flod,pz,py,px);

  // ---- native LOD0 region unroll, placed against the global winding ----
  if (argc>=9){
    long z0=atol(argv[5]),y0=atol(argv[6]),x0=atol(argv[7]); int d=atoi(argv[8]);
    int ppw=argc>9?atoi(argv[9]):32; double s=(double)(1<<flod);
    fprintf(stderr,"native LOD0 region %d^3 @ (%ld,%ld,%ld)...\n",d,z0,y0,x0);
    u8 *img=mca_load_region(path,0,(int)z0,(int)y0,(int)x0,d,d,d);
    if(img){
      size_t dn=(size_t)d*d*d;
      // native mask + LOCAL LOD0 winding: warm-start from upsampled coarse global
      // winding, then re-solve at native resolution (fixes inner-wrap distortion).
      u8 *lm=malloc(dn); for(size_t i=0;i<dn;i++) lm[i]=img[i]>=80; majority_filter(lm,lm,d,d,d,1);
      remove_small_components(lm,d,d,d,TOPO_CONN26,100);
      f32 *lw=malloc(dn*sizeof(f32));
      #pragma omp parallel for schedule(static)
      for(int z=0;z<d;z++)for(int y=0;y<d;y++)for(int x=0;x<d;x++){
        size_t i=((size_t)z*d+y)*d+x;
        lw[i]= lm[i]? trilin(prevw,pz,py,px,(z0+z)/s,(y0+y)/s,(x0+x)/s):0; }
      umbilicus lu; lu.n=2; lu.z=malloc(8);lu.y=malloc(8);lu.x=malloc(8);
      lu.z[0]=0;lu.z[1]=(f32)(d-1); lu.y[0]=lu.y[1]=(f32)(fy*0.5-y0); lu.x[0]=lu.x[1]=(f32)(fx*0.5-x0);
      wfield_params lp=winding_default_params(); lp.warm_start=1; lp.iters=30;
      fprintf(stderr,"refining winding at native LOD0...\n");
      winding_field_solve(lm,d,d,d,&lu,&lp,NULL,NULL,lw);
      double lmvf=winding_mvf(lm,lw,d,d,d,&lu,64,d/2);
      float wmin=1e30f,wmax=-1e30f; for(size_t i=0;i<dn;i++) if(lm[i]){if(lw[i]<wmin)wmin=lw[i];if(lw[i]>wmax)wmax=lw[i];}
      printf("region global winding %.2f..%.2f (~%.1f wraps), native MVF=%.4f\n",wmin,wmax,wmax-wmin,lmvf);
      int W=(int)((wmax-wmin)*ppw)+1; if(W<1)W=1; if(W>60000)W=60000; int H=d;
      double *acc=calloc((size_t)W*H,sizeof(double)); int *cnt=calloc((size_t)W*H,sizeof(int));
      #pragma omp parallel for schedule(static)
      for(int z=0;z<d;z++)for(int y=0;y<d;y++)for(int x=0;x<d;x++){
        size_t i=((size_t)z*d+y)*d+x; if(!lm[i])continue;
        int u=(int)((lw[i]-wmin)*ppw); if(u<0||u>=W)continue;
        size_t pidx=(size_t)z*W+u;
        #pragma omp atomic
        acc[pidx]+=img[i];
        #pragma omp atomic
        cnt[pidx]++;
      }
      free(lm);free(lw);free(lu.z);free(lu.y);free(lu.x);
      u8 *flat=calloc((size_t)W*H,1); size_t fill=0;
      for(size_t i=0;i<(size_t)W*H;i++){ if(cnt[i]){flat[i]=(u8)(acc[i]/cnt[i]);fill++;} }
      tiff_save_u8(outp,flat,1,H,W);
      printf("wrote native-res unroll %s (%dx%d), %.0f%% filled\n",outp,W,H,100.0*fill/((size_t)W*H));
      free(img);free(acc);free(cnt);free(flat);
    }
  }
  free(prevw);
  return 0;
}
