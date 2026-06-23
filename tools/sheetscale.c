/* sheetscale — multiscale structure-tensor sheetness + argmax-sigma FUSED-WRAP detector.
 * Phase 0a of docs/touching-sheets-plan.md.
 *
 * The two-peak (Sparrow) limit: two parallel sheets a gap `d` apart give TWO resolved sheetness peaks
 * when the tensor integration scale sigma < d/2, but FUSE into one thick slab when sigma > d/2. So at a
 * fixed single sigma (sheet_sep3d uses ~1.5) two touching wraps may merge purely as a scale artifact.
 * Here we run st_sheet_detect over a sigma SET and, per voxel, record the sigma that MAXIMIZES sheetness:
 *   - a LONE sheet peaks at its natural (small) scale,
 *   - a FUSED pair peaks at a LARGER scale (it only looks planar once smoothed into one slab).
 * So the ARGMAX-sigma map is itself a fused-wrap detector: voxels selecting a sigma above the regional
 * mode (the lone-sheet scale) are touch candidates. Reports a ground-truth-free "touch fraction" and
 * writes a PPM colored by selected scale (cool=lone/small -> warm=fused/large).
 *
 *   sheetscale ARC OUT lod z0 y0 x0 dz dy dx [zc=mid] [sigmas=0.7,1.0,1.5,2.0]
 *
 * RESULT (2026-06-23, L2 pherc0332): the argmax-sigma TOUCH detector is WEAK -- touch fraction is flat
 * ~0.30 across the compressed core, centered, AND delaminated-outer regions (a real detector must score
 * the core >> the delaminated outer), and the warm (large-sigma) voxels track sheet THICKNESS/contrast,
 * not specifically fusions. That is the thickness<->fusion confound: a fused pair and a genuinely thick
 * lone sheet are the SAME slab, locally indistinguishable -- the touching-sheets wall, re-confirmed from
 * the scale-space angle. BUT the SCALE finding is solid and actionable: lone sheets peak at sigma<=0.7
 * (70% mode at the smallest sigma), so sheet_sep3d's default sigma_tensor~1.5 OVER-SMOOTHS most sheets.
 * Use small-sigma detection (Phase 2a) + small-sigma normals for supervoxel orientation (Phase 1).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"

#define MAXS 8

static int air_threshold(const u8*v,size_t n){ long h[256]={0},t=0; double s=0; long z=0;
  for(size_t i=0;i<n;i++){int x=v[i]; if(x){s+=x;z++;} if(x>=1&&x<=254){h[x]++;t++;}}
  if(t<256)return z?(int)(0.5*s/z+0.5):1; double sum=0; for(int i=1;i<=254;i++)sum+=(double)i*h[i];
  double sB=0,wB=0,b=-1; int thr=1; for(int k=1;k<=254;k++){wB+=h[k]; if(!wB)continue; double wF=t-wB; if(wF<=0)break;
    sB+=(double)k*h[k]; double mB=sB/wB,mF=(sum-sB)/wF,bb=wB*wF*(mB-mF)*(mB-mF); if(bb>b){b=bb;thr=k;}} return thr; }

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [zc=mid] [sigmas=0.7,1.0,1.5,2.0]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int zc=argc>10?atoi(argv[10]):dz/2;
  float sig[MAXS]={0.7f,1.0f,1.5f,2.0f}; int ns=4;
  if(argc>11){ ns=0; char buf[256]; strncpy(buf,argv[11],255); buf[255]=0; char*t=strtok(buf,",");
    while(t&&ns<MAXS){ sig[ns++]=atof(t); t=strtok(NULL,","); } }
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr; nmat+=mask[p]; }
  fprintf(stderr,"region %dx%dx%d air<%d, %ld material; sigmas:",dz,dy,dx,athr,nmat);
  for(int s=0;s<ns;s++) fprintf(stderr," %.2f",sig[s]); fprintf(stderr,"\n");

  f32*vf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vf[p]=v[p];
  f32*sh=malloc(nn*sizeof(f32)),*shmax=calloc(nn,sizeof(f32)); u8*amax=calloc(nn,1);
  for(int s=0;s<ns;s++){ st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=sig[s];
    st_sheet_detect(vf,dz,dy,dx,&sp,sh,NULL);
    for(size_t p=0;p<nn;p++) if(mask[p] && sh[p]>shmax[p]){ shmax[p]=sh[p]; amax[p]=(u8)s; }
    fprintf(stderr,"  sigma=%.2f done\n",sig[s]); }
  free(sh); free(vf);

  // argmax-sigma histogram over material; mode = lone-sheet scale; touch frac = above-mode fraction
  long hist[MAXS]={0}; for(size_t p=0;p<nn;p++) if(mask[p]&&shmax[p]>1e-4) hist[amax[p]]++;
  long tot=0; for(int s=0;s<ns;s++) tot+=hist[s];
  int mode=0; for(int s=1;s<ns;s++) if(hist[s]>hist[mode]) mode=s;
  long above=0; for(int s=mode+1;s<ns;s++) above+=hist[s];
  fprintf(stderr,"argmax-sigma histogram (material with sheetness>0):\n");
  for(int s=0;s<ns;s++) fprintf(stderr,"  sigma=%.2f : %8ld (%.1f%%)%s\n",sig[s],hist[s],
    tot?100.0*hist[s]/tot:0, s==mode?"  <- mode (lone-sheet scale)":"");
  printf("TOUCH FRACTION (argmax-sigma above lone-sheet mode %.2f): %.3f  (%ld/%ld voxels select a larger scale = fused-wrap candidates)\n",
    sig[mode],tot?(double)above/tot:0,above,tot);

  // PPM: color material by selected scale, cool(small/lone)->warm(large/fused); dim CT background
  u8*rgb=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
    if(mask[p]&&shmax[p]>1e-4){ double t=ns>1?(double)amax[p]/(ns-1):0;   // 0=small 1=large
      rgb[o]=(u8)(40+215*t); rgb[o+1]=(u8)(60*(1-fabs(2*t-1))); rgb[o+2]=(u8)(40+215*(1-t)); } }
  char fn[700]; snprintf(fn,sizeof fn,"%s_scale.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb);free(shmax);free(amax);free(mask);free(v); mca_close(arc); return 0;
}
