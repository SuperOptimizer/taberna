/* unroll_wind — flatten a region into an unrolled papyrus image using a precomputed winding field.
 *
 * Reads a merged winding volume (the _vol.f32 / merged.f32 format: int32 {dz,dy,dx,z0,y0,x0} then
 * dz*dy*dx float32, NaN off-material, in the archive's `lod` grid) and the matching CT intensity
 * from the archive. Vertical axis is always z. Two horizontal modes:
 *
 *  RADIAL (default): horizontal = winding w directly. sheet_sep3d's winding is ~RADIAL (constant
 *    along each sheet), so this collapses the azimuth -> a sheet-vs-z QUALITY reslice, not a text
 *    image. Max-composited.
 *
 *  SPIRAL (pass CY CX = umbilicus center in merged-local y,x): reconstructs the true Archimedean
 *    unroll coordinate  w_spiral = w + atan2(y-cy, x-cx)/2pi , spreading each wrap across the
 *    azimuth (the along-sheet / text direction). Cut at theta=+/-pi. Mean-composited (matches
 *    tools/unroll). NOTE at L2 (4x downsampled) ink is not resolvable -- this shows flattened
 *    SHEET structure, not letters.
 *
 * The inner umbilicus core (winding ill-defined as r->0) is trimmed via the 1..99.5 percentile of
 * the render coordinate.
 *
 *   unroll_wind ARCHIVE lod WINDING.f32 OUT.pgm [SAMP=160] [CY CX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"

static int cmp_f32(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return (x>y)-(x<y); }

// NaN-aware trilinear sample of a winding volume (averages only finite corners; NAN if none).
static double wtrilin(const float*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0,acc=0,wt=0;
  int zs[2]={z0,z1},ys[2]={y0,y1},xs[2]={x0,x1}; double zw[2]={1-dz,dz},yw[2]={1-dy,dy},xw[2]={1-dx,dx};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){
    double v=w[((size_t)zs[a]*ny+ys[b])*nx+xs[c]]; if(!isfinite(v))continue;
    double g=zw[a]*yw[b]*xw[c]; acc+=g*v; wt+=g; }
  return wt>1e-9?acc/wt:NAN;
}

int main(int argc,char**argv){
  if(argc<5){ fprintf(stderr,"usage: %s ARCHIVE lod WINDING.f32 OUT.pgm [SAMP=160] [CY CX]\n",argv[0]); return 2; }
  const char*arc=argv[1]; int lod=atoi(argv[2]); const char*wf=argv[3],*outp=argv[4];
  int SAMP=argc>5?atoi(argv[5]):160;
  int spiral = (argc>7);
  double cy=spiral?atof(argv[6]):0, cx=spiral?atof(argv[7]):0;
  FILE*f=fopen(wf,"rb"); if(!f){ fprintf(stderr,"open %s fail\n",wf); return 1; }
  int hdr[6]; if(fread(hdr,sizeof(int),6,f)!=6){ fprintf(stderr,"bad header\n"); return 1; }
  int dz=hdr[0],dy=hdr[1],dx=hdr[2],z0=hdr[3],y0=hdr[4],x0=hdr[5];
  size_t nn=(size_t)dz*dy*dx; float*w=malloc(nn*sizeof(float));
  if(fread(w,sizeof(float),nn,f)!=nn){ fprintf(stderr,"short winding volume\n"); return 1; } fclose(f);
  fprintf(stderr,"winding %dx%dx%d @ lod%d (z%d y%d x%d) mode=%s\n",dz,dy,dx,lod,z0,y0,x0,spiral?"SPIRAL":"radial");

  const double TWO_PI=6.283185307179586;
  // FINE mode (argv[8]=finer archive, [9]=scale, [10]=ZA, [11]=ZN): sample the winding
  // (trilinear-upsampled) against a FINER CT archive over a thin z-band [ZA,ZA+ZN) of the winding,
  // for ink-scale resolution the L2 winding doesn't carry. The finer archive is ==scale x the
  // winding archive (padded), so finevoxel = scale * windingvoxel. Requires spiral (CY CX).
  int fine = (argc>11);
  const char*ctarc = fine?argv[8]:arc; int scale=fine?atoi(argv[9]):1;
  int za=fine?atoi(argv[10]):0, zn=fine?atoi(argv[11]):dz;
  if(fine && !spiral){ fprintf(stderr,"fine mode needs CY CX (spiral)\n"); return 2; }
  if(za<0)za=0; if(za+zn>dz)zn=dz-za;

  // winding range = 1..99.5 percentile of finite render coords over the rendered z-band
  double wmin,wmax;
  {
    int zlo=fine?za:0, zhi=fine?za+zn:dz;
    size_t cap=(size_t)(zhi-zlo)*dy*dx; float*tmp=malloc((cap?cap:1)*sizeof(float)); size_t k=0;
    for(int z=zlo;z<zhi;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){
      size_t p=((size_t)z*dy+y)*dx+x; if(!isfinite(w[p]))continue;
      double rc=w[p]; if(spiral) rc+=atan2((double)y-cy,(double)x-cx)/TWO_PI; tmp[k++]=(float)rc; }
    qsort(tmp,k,sizeof(float),cmp_f32);
    wmin=k?tmp[(size_t)(0.010*k)]:0; wmax=k?tmp[(size_t)(0.995*(k-1))]:1; free(tmp);
  }
  if(wmax<=wmin){ fprintf(stderr,"degenerate winding range\n"); return 1; }
  int UW=(int)((wmax-wmin)*SAMP)+1; if(UW>120000)UW=120000;

  mca_handle*h=mca_open(ctarc); if(!h){ fprintf(stderr,"open arc %s fail\n",ctarc); return 1; }

  if(fine){
    int Fz0=(z0+za)*scale, Fy0=y0*scale, Fx0=x0*scale, Fdz=zn*scale, Fdy=dy*scale, Fdx=dx*scale;
    fprintf(stderr,"FINE render: CT %s scale %d, fine region z%d y%d x%d + %dx%dx%d -> %d cols x %d rows\n",
            ctarc,scale,Fz0,Fy0,Fx0,Fdz,Fdy,Fdx,UW,Fdz);
    u8*ct=mca_read(h,0,Fz0,Fy0,Fx0,Fdz,Fdy,Fdx); if(!ct){ fprintf(stderr,"fine CT read fail\n"); return 1; }
    double*acc=calloc((size_t)Fdz*UW,sizeof(double)); int*cnt=calloc((size_t)Fdz*UW,sizeof(int));
    // parallel over zf: each zf writes a disjoint output row [zf*UW,(zf+1)*UW) -> no atomics needed
    #pragma omp parallel for schedule(static)
    for(int zf=0;zf<Fdz;zf++){ double wlz=za+(double)zf/scale;
      for(int yf=0;yf<Fdy;yf++){ double wly=(double)yf/scale;
        for(int xf=0;xf<Fdx;xf++){ double wlx=(double)xf/scale;
          double wv=wtrilin(w,dz,dy,dx,wlz,wly,wlx); if(!isfinite(wv))continue;
          double rc=wv+atan2(wly-cy,wlx-cx)/TWO_PI;
          int u=(int)((rc-wmin)*SAMP); if(u<0||u>=UW)continue;
          size_t o=(size_t)zf*UW+u; acc[o]+=ct[((size_t)zf*Fdy+yf)*Fdx+xf]; cnt[o]++; } } }
    u8*img=calloc((size_t)Fdz*UW,1);
    for(size_t i=0;i<(size_t)Fdz*UW;i++) if(cnt[i]) img[i]=(u8)(acc[i]/cnt[i]);
    FILE*o=fopen(outp,"wb"); if(!o){ fprintf(stderr,"write %s fail\n",outp); return 1; }
    fprintf(o,"P5\n%d %d\n255\n",UW,Fdz); fwrite(img,1,(size_t)Fdz*UW,o); fclose(o);
    fprintf(stderr,"wrote %s (%dx%d)\n",outp,UW,Fdz);
    free(w); free(img); free(acc); free(cnt); free(ct); mca_close(h); return 0;
  }

  u8*ct=mca_read(h,lod,z0,y0,x0,dz,dy,dx); if(!ct){ fprintf(stderr,"CT read fail\n"); return 1; }
  fprintf(stderr,"render-coord %.1f..%.1f (%.0f wraps) -> %d cols x %d rows\n",wmin,wmax,wmax-wmin,UW,dz);
  double*acc=calloc((size_t)dz*UW,sizeof(double)); int*cnt=calloc((size_t)dz*UW,sizeof(int));
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){
    size_t p=((size_t)z*dy+y)*dx+x; float wv=w[p]; if(!isfinite(wv))continue;
    double rc=wv; if(spiral) rc+=atan2((double)y-cy,(double)x-cx)/TWO_PI;
    int u=(int)((rc-wmin)*SAMP); if(u<0||u>=UW)continue;
    size_t o=(size_t)z*UW+u; acc[o]+=ct[p]; cnt[o]++; }
  u8*img=calloc((size_t)dz*UW,1);
  for(size_t i=0;i<(size_t)dz*UW;i++) if(cnt[i]) img[i]=(u8)(acc[i]/cnt[i]);

  FILE*o=fopen(outp,"wb"); if(!o){ fprintf(stderr,"write %s fail\n",outp); return 1; }
  fprintf(o,"P5\n%d %d\n255\n",UW,dz); fwrite(img,1,(size_t)dz*UW,o); fclose(o);
  fprintf(stderr,"wrote %s (%dx%d)\n",outp,UW,dz);
  free(w); free(img); free(acc); free(cnt); mca_close(h); return 0;
}
