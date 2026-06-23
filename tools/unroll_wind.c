/* unroll_wind — flatten a region into an unrolled papyrus image using a precomputed winding field.
 *
 * Reads a merged winding volume (the _vol.f32 / merged.f32 format: int32 {dz,dy,dx,z0,y0,x0} then
 * dz*dy*dx float32, NaN off-material, in the archive's `lod` grid) and the matching CT intensity
 * from the archive, then renders the unroll: horizontal axis = winding w (SAMP columns per wrap),
 * vertical axis = z. Every material voxel deposits its CT intensity at column round((w-wmin)*SAMP);
 * voxels sharing a (w,z) bin (the radial thickness of one sheet at one azimuth) are MAX-composited,
 * so each sheet lays flat as a vertical band and crossing wraps reads left-to-right.
 *
 * The inner umbilicus core (winding ill-defined as r->0, often <0 or huge) wastes columns; by
 * default we clip the winding range to the 1..99.5 percentile of finite voxels. Pass explicit
 * WMIN WMAX to override.
 *
 *   unroll_wind ARCHIVE lod WINDING.f32 OUT.pgm [SAMP=160] [WMIN WMAX]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"

static int cmp_f32(const void*a,const void*b){ float x=*(const float*)a,y=*(const float*)b; return (x>y)-(x<y); }

int main(int argc,char**argv){
  if(argc<5){ fprintf(stderr,"usage: %s ARCHIVE lod WINDING.f32 OUT.pgm [SAMP=160]\n",argv[0]); return 2; }
  const char*arc=argv[1]; int lod=atoi(argv[2]); const char*wf=argv[3],*outp=argv[4];
  int SAMP=argc>5?atoi(argv[5]):160;
  FILE*f=fopen(wf,"rb"); if(!f){ fprintf(stderr,"open %s fail\n",wf); return 1; }
  int hdr[6]; if(fread(hdr,sizeof(int),6,f)!=6){ fprintf(stderr,"bad header\n"); return 1; }
  int dz=hdr[0],dy=hdr[1],dx=hdr[2],z0=hdr[3],y0=hdr[4],x0=hdr[5];
  size_t nn=(size_t)dz*dy*dx; float*w=malloc(nn*sizeof(float));
  if(fread(w,sizeof(float),nn,f)!=nn){ fprintf(stderr,"short winding volume\n"); return 1; } fclose(f);
  fprintf(stderr,"winding %dx%dx%d @ lod%d (z%d y%d x%d)\n",dz,dy,dx,lod,z0,y0,x0);

  mca_handle*h=mca_open(arc); if(!h){ fprintf(stderr,"open arc fail\n"); return 1; }
  u8*ct=mca_read(h,lod,z0,y0,x0,dz,dy,dx); if(!ct){ fprintf(stderr,"CT read fail\n"); return 1; }

  // winding range: explicit WMIN WMAX, else 1..99.5 percentile of finite voxels (trims the core)
  double wmin,wmax;
  if(argc>7){ wmin=atof(argv[6]); wmax=atof(argv[7]); }
  else{
    size_t nf=0; for(size_t p=0;p<nn;p++) if(isfinite(w[p])) nf++;
    float*tmp=malloc((nf?nf:1)*sizeof(float)); size_t k=0;
    for(size_t p=0;p<nn;p++) if(isfinite(w[p])) tmp[k++]=w[p];
    qsort(tmp,nf,sizeof(float),cmp_f32);
    wmin=nf?tmp[(size_t)(0.010*nf)]:0; wmax=nf?tmp[(size_t)(0.995*(nf-1))]:1; free(tmp);
  }
  if(wmax<=wmin){ fprintf(stderr,"degenerate winding range\n"); return 1; }
  int UW=(int)((wmax-wmin)*SAMP)+1;
  fprintf(stderr,"winding %.1f..%.1f (%.0f wraps) -> %d cols x %d rows\n",wmin,wmax,wmax-wmin,UW,dz);
  u8*img=calloc((size_t)dz*UW,1);
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){
    size_t p=((size_t)z*dy+y)*dx+x; float wv=w[p]; if(!isfinite(wv))continue;
    int u=(int)((wv-wmin)*SAMP); if(u<0||u>=UW)continue;
    u8 c=ct[p]; size_t o=(size_t)z*UW+u; if(c>img[o])img[o]=c; }

  FILE*o=fopen(outp,"wb"); if(!o){ fprintf(stderr,"write %s fail\n",outp); return 1; }
  fprintf(o,"P5\n%d %d\n255\n",UW,dz); fwrite(img,1,(size_t)dz*UW,o); fclose(o);
  fprintf(stderr,"wrote %s (%dx%d)\n",outp,UW,dz);
  free(w); free(img); mca_close(h); return 0;
}
