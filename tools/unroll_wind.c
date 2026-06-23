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

  mca_handle*h=mca_open(arc); if(!h){ fprintf(stderr,"open arc fail\n"); return 1; }
  u8*ct=mca_read(h,lod,z0,y0,x0,dz,dy,dx); if(!ct){ fprintf(stderr,"CT read fail\n"); return 1; }

  // render coordinate per voxel: radial = w; spiral = w + atan2(dy,dx)/2pi
  const double TWO_PI=6.283185307179586;
  // winding range = 1..99.5 percentile of finite render coords (trims the ill-defined core)
  double wmin,wmax;
  {
    size_t nf=0; for(size_t p=0;p<nn;p++) if(isfinite(w[p])) nf++;
    float*tmp=malloc((nf?nf:1)*sizeof(float)); size_t k=0;
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){
      size_t p=((size_t)z*dy+y)*dx+x; if(!isfinite(w[p]))continue;
      double rc=w[p]; if(spiral) rc+=atan2((double)y-cy,(double)x-cx)/TWO_PI; tmp[k++]=(float)rc; }
    qsort(tmp,k,sizeof(float),cmp_f32);
    wmin=k?tmp[(size_t)(0.010*k)]:0; wmax=k?tmp[(size_t)(0.995*(k-1))]:1; free(tmp);
  }
  if(wmax<=wmin){ fprintf(stderr,"degenerate winding range\n"); return 1; }
  int UW=(int)((wmax-wmin)*SAMP)+1; if(UW>60000)UW=60000;
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
