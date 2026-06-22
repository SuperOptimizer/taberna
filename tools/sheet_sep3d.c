/* sheet_sep3d.c -- 3D scroll winding over a sub-volume. ONE global winding-number solve.
 *
 * Hard-won architecture (see commit log): 3D connected-component segmentation PERCOLATES
 * across touching wraps (the touching-sheet frontier) -- any bridge anywhere in z merges
 * wraps into one node -> flat winding; thinning the ridge to avoid that SHATTERS each wrap
 * into fragments instead. The working approach detects sheets per-slice in 2D (clean), and
 * the WINDING NUMBER comes from RADIAL valley-counted step-edges (the reliable engine),
 * solved globally. Ties (same-wrap, valley-count 0) FLATTEN the radial gradient because
 * "same wrap" is unreliable at touches -- so ties are weak/off by default (WEQ, z-tie args).
 *
 *   1. read a sub-volume; bimodal air-cut; per-z umbilicus axis -> in-plane radial dir.
 *   2. 3D ridge (radial intensity max) + recto (material, air outward) detection.
 *   3. PER-SLICE 2D cc_label (no z-percolation) -> nodes.
 *   4. radial walk + VALLEY-COUNT winding step (keep walking through same-wrap fragments,
 *      stop at the first valley crossing) -> step-edges; optional same-wrap/z ties.
 *   5. ONE global IRLS least-squares winding solve + a WEAK RADIUS PRIOR. The radial
 *      step-edges give precise local structure, but the graph is many disconnected
 *      components (per-slice/per-angle) that float with arbitrary integer offsets and
 *      cancel when pooled. A weak pull toward winding ~= (radius-r_inner)/pitch fixes each
 *      component's offset from its (reliable) radius WITHOUT wrap-bridging z/same-wrap ties.
 *   6. dense anisotropic 3D propagation (radial-blocked barriers, z-coupled, soft anchors).
 *
 * THICK ridges (2D-quality) + per-slice 2D labels + tangential ties (same-wrap, can't bridge
 * wraps) + radius prior -> the node winding is MONOTONIC in radius and matches the pitch
 * geometry (span ~= radial-span/pitch). Recto is sparse/flat in compressed cores and flattens
 * the fit -> RIDGE-ONLY by default (argv[15]=1 to add recto for delaminated outer regions).
 *
 * usage: sheet_sep3d ARCHIVE OUTBASE lod z0 y0 x0 dz dy dx [minseg=40] [pitch=auto] [LAM=0.6]
 *        [WEQ=3] [z-tie=0] [use_recto=0] [LPRI=0.15]
 * STATUS: node winding CORRECT; dense backward-switch 316->139/1000 voxels (radial diffusion, barriers-only block). Outer region + 3D structure-tensor normals are the remaining polish.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "io/mca.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"
#include "segmentation/sheet_tensor.h"

typedef unsigned char u8; typedef unsigned int u32; typedef float f32;
#define IDX(z,y,x) (((size_t)(z)*dy+(y))*dx+(x))
static double tnow(void){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return t.tv_sec+t.tv_nsec*1e-9; }
static double g_tl;
#define PHASE(n) do{ fprintf(stderr,"[t] %-22s %.2fs\n",n,tnow()-g_tl); g_tl=tnow(); }while(0)

static int air_threshold(const u8 *v, size_t n){
  long h[256]={0}; long tot=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x>=1&&x<=254){ h[x]++; tot++; } }
  if(tot<256) return 1;
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sumB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){ wB+=h[t]; if(wB==0)continue; double wF=(double)tot-wB; if(wF<=0)break;
    sumB+=(double)t*h[t]; double mB=sumB/wB,mF=(sum-sumB)/wF,btw=wB*wF*(mB-mF)*(mB-mF);
    if(btw>best){best=btw;thr=t;} }
  return thr;
}

// aggregate autocorr pitch on one z-slice (non-square dy x dx)
static double measure_pitch(const u8 *sl, int dy, int dx, double cyf, double cxf, int athr){
  int R=(int)(0.45*(dy<dx?dy:dx)); if(R>800)R=800; if(R<16) return 8.0;
  double *prof=malloc(R*sizeof(double)),*ac=malloc(R*sizeof(double));
  double *agg=calloc(R,sizeof(double)); long *cnt=calloc(R,sizeof(long)); int nr=0;
  for(int a=0;a<360;a++){ double th=a*(3.14159265/180.0),uy=sin(th),ux=cos(th); int n=0;
    for(int k=0;k<R;k++){ int yy=(int)lround(cyf+uy*k),xx=(int)lround(cxf+ux*k);
      double val=(yy<0||yy>=dy||xx<0||xx>=dx)?0:sl[(size_t)yy*dx+xx]; prof[k]=val<athr?0:val-athr; if(prof[k]>0)n=k; }
    if(n<24)continue; double mean=0; for(int k=0;k<=n;k++)mean+=prof[k]; mean/=(n+1);
    for(int k=0;k<=n;k++)prof[k]-=mean;
    for(int L=0;L<=n/2;L++){ double s=0; for(int k=0;k+L<=n;k++)s+=prof[k]*prof[k+L]; ac[L]=s; }
    if(ac[0]<=0)continue; for(int L=0;L<=n/2;L++){ agg[L]+=ac[L]/ac[0]; cnt[L]++; } nr++; }
  free(prof);free(ac); if(nr<8){free(agg);free(cnt);return 8.0;}
  for(int L=0;L<R;L++) if(cnt[L]>0) agg[L]/=cnt[L];
  int best=-1; double bv=-1; for(int L=6;L<R/3;L++)
    if(agg[L]>=agg[L-1]&&agg[L]>=agg[L-2]&&agg[L]>=agg[L+1]&&agg[L]>=agg[L+2]&&agg[L]>bv){bv=agg[L];best=L;}
  free(agg);free(cnt); return best>=6?(double)best:8.0;
}

// count inter-sheet valleys along the IN-PLANE radial ray at fixed z (wrap-aware, pitch-gated)
// wdcal = total-variation per wrap (windingDistance calibration); >0 enables the rescue below.
// tvsm = radial smoothing radius for the TV signal (suppresses within-sheet texture so the raw-CT
// gradient integral measures only the pitch-scale ridge->valley swings = sheet crossings).
static int count_valleys(const f32 *vsz, int dy, int dx, double y, double x, double uy, double ux,
                         int kstop, double prom, double pitchmin, double wdcal, int tvsm){
  if(kstop<2) return 1;
  double peakprom=1.8*prom; double pos[260]; signed char typ[260]; int ne=0;
  pos[ne]=0; typ[ne]=1; ne++;
  double prof[800]; int np=kstop<799?kstop:799;   // buffer large enough for full-radius ground-truth rays too
  for(int t=0;t<=np;t++){ int iy=(int)lround(y+uy*t),ix=(int)lround(x+ux*t);
    prof[t]=(iy<0||iy>=dy||ix<0||ix>=dx)?0.0:(double)vsz[(size_t)iy*dx+ix]; }
  double cmax=-1e30,cmin=1e30,maxp=0,minp=0; int dir=0;
  for(int t=0;t<=np;t++){ double val=prof[t];
    if(t==0){cmax=cmin=val;maxp=minp=0;continue;}
    if(dir>=0){ if(val>cmax){cmax=val;maxp=t;} if(cmax-val>=peakprom){ if(ne<258&&typ[ne-1]!=1){pos[ne]=maxp;typ[ne]=1;ne++;} dir=-1;cmin=val;minp=t; } }
    if(dir<=0){ if(val<cmin){cmin=val;minp=t;} if(val-cmin>=prom){ if(ne<258&&typ[ne-1]!=-1){pos[ne]=minp;typ[ne]=-1;ne++;} dir=1;cmax=val;maxp=t; } } }
  if(typ[ne-1]!=1){ pos[ne]=np; typ[ne]=1; ne++; }
  int nv=0; for(int i=1;i<ne-1;i++) if(typ[i]==-1 && pos[i+1]-pos[i-1]>=pitchmin) nv++;
  // windingDistance RESCUE (adapted from VC3D FiberIntersections::windingDistance): at a touch the
  // inter-sheet valley can be too shallow for prominence detection (nv under-counts -> sheets get
  // merged -> backward switch). Integrate the gradient of the RADIALLY-SMOOTHED profile (texture
  // removed, only pitch-scale swings survive); if valley-count found nothing but smoothed-TV is
  // >=0.7 wrap's worth, the segment really did cross a (weak) boundary.
  if(wdcal>0 && nv==0){
    double tv=0,prev=prof[0];
    for(int t=1;t<=np;t++){ double s=0;int c=0; for(int d=-tvsm;d<=tvsm;d++){int u=t+d;if(u<0||u>np)continue;s+=prof[u];c++;}
      double sm=s/c; tv+=fabs(sm-prev); prev=sm; }
    if(tv>=0.7*wdcal) nv=1;
  }
  return nv;
}

// RANSAC/IRLS umbilicus-from-normals (adapted from VC3D normalgridtools::align_and_extract_umbilicus).
// Sheet normals point radially, so each material point's gradient defines a LINE that passes
// through the scroll center. The center is the robust least-squares intersection of those lines:
// minimise sum_i w_i (c-p_i)^T (I - n n^T) (c-p_i)  ->  (sum w M) c = (sum w M p), a 2x2 solve,
// IRLS-reweighted (Tukey on point-line distance, weight by gradient magnitude). Works even when
// the crop is OFF-CENTER (lines still intersect outside the crop), where the coarse-LOD5 estimate
// extrapolates poorly. Returns the per-line-distance inlier RMS for diagnostics.
static double refine_center_slice(const f32*vs,const u8*v,int athr,int dy,int dx,
                                  double*cy_io,double*cx_io,double pitch){
  double gthr=0.05*athr; if(gthr<2)gthr=2; double scale=0.5*pitch; if(scale<4)scale=4;
  double cyv=*cy_io,cxv=*cx_io,rms=0;
  for(int rep=0;rep<5;rep++){
    double A00=0,A01=0,A11=0,b0=0,b1=0,sd=0; long ni=0;
    for(int y=1;y<dy-1;y++)for(int x=1;x<dx-1;x++){ size_t p=(size_t)y*dx+x; if(v[p]<athr)continue;
      double gy=0.5*((double)vs[p+dx]-vs[p-dx]),gx=0.5*((double)vs[p+1]-vs[p-1]); double g=sqrt(gy*gy+gx*gx); if(g<gthr)continue;
      double ny=gy/g,nx=gx/g, M00=1-ny*ny,M01=-ny*nx,M11=1-nx*nx;
      double dyp=cyv-y,dxp=cxv-x, d2=dyp*(M00*dyp+M01*dxp)+dxp*(M01*dyp+M11*dxp), d=sqrt(d2>0?d2:0);
      double tw=(d<scale)?(1-(d/scale)*(d/scale))*(1-(d/scale)*(d/scale)):0; if(tw<=0)continue; double w=g*tw;
      A00+=w*M00;A01+=w*M01;A11+=w*M11; b0+=w*(M00*y+M01*x); b1+=w*(M01*y+M11*x); sd+=tw*d2; ni++; }
    double det=A00*A11-A01*A01; if(fabs(det)<1e-6||ni<50) break;
    cyv=(A11*b0-A01*b1)/det; cxv=(A00*b1-A01*b0)/det; rms=sqrt(sd/ni);
  }
  *cy_io=cyv; *cx_io=cxv; return rms;
}

static int cmp_long(const void*a,const void*b){ long x=*(const long*)a,y=*(const long*)b; return (x>y)-(x<y); }
static int cmp_dbl(const void*a,const void*b){ double x=*(const double*)a,y=*(const double*)b; return (x>y)-(x<y); }
// aggregate packed (key<<4|step) edge keys -> sparse (a,b,count,step) keeping run>=minc
static long agg_edges(long*raw,long nraw,int K,int**ao,int**bo,long**wo,long**so,long minc){
  if(nraw==0){ *ao=malloc(4);*bo=malloc(4);*wo=malloc(8); if(so)*so=malloc(8); return 0; }
  qsort(raw,nraw,sizeof(long),cmp_long);
  long ne=0; for(long i=0;i<nraw;){ long key=raw[i]>>4; long j=i; while(j<nraw&&(raw[j]>>4)==key)j++; if(j-i>=minc)ne++; i=j; }
  int*ea=malloc((size_t)(ne?ne:1)*sizeof(int)),*eb=malloc((size_t)(ne?ne:1)*sizeof(int));
  long*ew=malloc((size_t)(ne?ne:1)*sizeof(long)),*es=so?malloc((size_t)(ne?ne:1)*sizeof(long)):NULL;
  long e=0; for(long i=0;i<nraw;){ long key=raw[i]>>4; long j=i,ss=0; while(j<nraw&&(raw[j]>>4)==key){ ss+=raw[j]&15; j++; } long c=j-i;
    if(c>=minc){ ea[e]=(int)(key/K); eb[e]=(int)(key%K); ew[e]=c; if(es)es[e]=lround((double)ss/c); e++; } i=j; }
  *ao=ea;*bo=eb;*wo=ew; if(so)*so=es; return ne;
}

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARCHIVE OUTBASE lod z0 y0 x0 dz dy dx [minseg=40] [pitch=auto]\n",argv[0]); return 2; }
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]);
  int dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int minseg=argc>10?atoi(argv[10]):40; double pitch_arg=argc>11?atof(argv[11]):0;
  int use_recto=argc>15?atoi(argv[15]):0;   // recto (delamination rims) is sparse/flat in compressed cores; 0 = ridge-only
  g_tl=tnow();
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);

  // umbilicus axis from coarse LOD5
  int cl=5; double cs=(double)(1<<cl); int ccz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,ccz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse fail\n");return 1;}
  size_t ccn=(size_t)ccz*ccy*ccx; u8*ccm=malloc(ccn); for(size_t i=0;i<ccn;i++)ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,ccz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;} free(ccm);
  double ls=(double)(1<<lod);
  // per-z center in sub-volume coords
  double *cy=malloc(dz*sizeof(double)),*cx=malloc(dz*sizeof(double));
  for(int zz=0;zz<dz;zz++){ double coarse_z=(double)(z0+zz)*ls/cs; if(coarse_z<0)coarse_z=0; if(coarse_z>ccz-1)coarse_z=ccz-1;
    f32 ucy,ucx; umbilicus_center(&umb,(f32)coarse_z,&ucy,&ucx); cy[zz]=ucy*cs/ls - y0; cx[zz]=ucx*cs/ls - x0; }
  fprintf(stderr,"archive %dx%dx%d; subvol z%ld y%ld x%ld + %dx%dx%d; center z0=(%.0f,%.0f) zmid=(%.0f,%.0f)\n",
          fz,fy,fx,z0,y0,x0,dz,dy,dx,cy[0],cx[0],cy[dz/2],cx[dz/2]);

  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,dz,dy,dx); if(!v){fprintf(stderr,"subvol read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx;
  int athr=air_threshold(v,nn); fprintf(stderr,"bimodal air threshold v<%d\n",athr); PHASE("read+threshold");

  // smoothed vs (per-slice 2D box, 2 iters), air->0
  f32*vs=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vs[p]=v[p];
  { f32*t=malloc(nn*sizeof(f32));
    for(int it=0;it<2;it++){
      #pragma omp parallel for schedule(static)
      for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); double s=vs[p];int c=1;
        if(x>0){s+=vs[p-1];c++;}if(x<dx-1){s+=vs[p+1];c++;}if(y>0){s+=vs[p-dx];c++;}if(y<dy-1){s+=vs[p+dx];c++;}
        t[p]=(f32)(s/c);} memcpy(vs,t,nn*sizeof(f32)); } free(t); }
  for(size_t p=0;p<nn;p++) if(v[p]<athr) vs[p]=0;

  double pitch = pitch_arg>0? pitch_arg : measure_pitch(v+(size_t)(dz/2)*dy*dx,dy,dx,cy[dz/2],cx[dz/2],athr);
  int walkr=(int)(1.6*pitch+0.5); if(walkr<8)walkr=8; if(walkr>120)walkr=120;
  int kmin=(int)(0.5*pitch+0.5); if(kmin<2)kmin=2;
  int tier=(int)(0.7*pitch+0.5); if(tier<3)tier=3; if(tier>30)tier=30;
  double vprom=0.35*athr; if(vprom<6)vprom=6; double pitchmin=0.20*pitch;
  fprintf(stderr,"pitch=%.1f walk=%d kmin=%d tie=%d vprom=%.1f\n",pitch,walkr,kmin,tier,vprom);

  // CLASSICAL SURFACE-PROBABILITY FIELD (no ML): structure-tensor sheetness in [0,1], the non-ML
  // equivalent of VC3D's U-Net surface prediction. VC3D's gradient/TV methods failed on raw CT
  // because raw intensity has within-sheet texture; the sheetness field is high on sheets, low
  // between, texture suppressed -> windingDistance etc. should work on IT. argv[24]=usesheet.
  int usesheet=argc>24?atoi(argv[24]):1; f32*sh=NULL;   // default ON: drives detection from the clean
  // classical sheetness field -> backward-switch -16..-24% in ALL regimes (core AND delaminated),
  // better z-coherence + wrap-count accuracy. argv[24]=0 reverts to raw-intensity detection.
  if(usesheet){ sh=malloc(nn*sizeof(f32)); f32*vf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vf[p]=v[p];
    st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=(f32)(pitch*0.12); if(sp.sigma_tensor<1)sp.sigma_tensor=1; if(sp.sigma_tensor>3)sp.sigma_tensor=3;
    st_sheet_detect(vf,dz,dy,dx,&sp,sh,NULL); free(vf);
    for(size_t p=0;p<nn;p++) if(v[p]<athr) sh[p]=0;   // mask air
    double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p];
    for(size_t p=0;p<nn;p++) sh[p]=(f32)(sh[p]/(smx>1e-6?smx:1)*255.0);   // scale to ~[0,255] for the valley/peak detectors
    fprintf(stderr,"sheetness field computed (sigma_tensor=%.1f, scaled to 0..255)\n",sp.sigma_tensor); }
  // SIGNAL driving detection (ridge/bar) + valley-counting: the clean sheetness field when usesheet,
  // else raw smoothed intensity. Sheets are MAXIMA in both (bright papyrus / high sheetness),
  // inter-sheet gaps are MINIMA -> the radial max/min detector and the prominence valley-counter
  // work unchanged; only the prominence + ridge-floor thresholds adapt to the signal's scale.
  f32*sig = usesheet? sh : vs;
  double sigprom = usesheet? 0.25*255.0 : vprom;    // inter-sheet valley prominence on the signal
  double ridgethr= usesheet? 0.12*255.0 : athr;     // "is on a sheet" floor for the ridge

  // refine per-z umbilicus from sheet normals (robust line intersection); keeps coarse estimate
  // as init/fallback. Sanity-clamp: reject a refined center that bolts >0.6*dim from coarse.
  int umbref=argc>21?atoi(argv[21]):1;
  if(umbref){ double *ny0=malloc(dz*sizeof(double)),*nx0=malloc(dz*sizeof(double)); double mr=0; int nref=0;
    for(int z=0;z<dz;z++){ double cyz=cy[z],cxz=cx[z];
      double rms=refine_center_slice(vs+(size_t)z*dy*dx,v+(size_t)z*dy*dx,athr,dy,dx,&cyz,&cxz,pitch);
      double clampr=0.6*((dy<dx)?dy:dx);
      if(fabs(cyz-cy[z])<clampr && fabs(cxz-cx[z])<clampr && rms>0){ ny0[z]=cyz; nx0[z]=cxz; mr+=rms; nref++; }
      else { ny0[z]=cy[z]; nx0[z]=cx[z]; } }
    // smooth refined centers across z (window +/-4) -> robust to per-slice noise
    for(int z=0;z<dz;z++){ double sy=0,sx=0;int n=0; for(int d=-4;d<=4;d++){int zz=z+d;if(zz<0||zz>=dz)continue;sy+=ny0[zz];sx+=nx0[zz];n++;} cy[z]=sy/n;cx[z]=sx/n; }
    fprintf(stderr,"umbilicus refined from normals: %d/%d slices, mean line-RMS=%.2f, center z0=(%.0f,%.0f) zmid=(%.0f,%.0f)\n",
      nref,dz,nref?mr/nref:0,cy[0],cx[0],cy[dz/2],cx[dz/2]);
    free(ny0);free(nx0); }

  // windingDistance calibration: TV (gradient integral) accrued over one wrap. Sample radial
  // segments of length `pitch` at many material points (mid-z) and take the median TV -> the
  // per-wrap gradient budget, so count_valleys can rescue prominence-missed crossings at touches.
  int wdrescue=argc>22?atoi(argv[22]):0; double wdcal=0;
  int tvsm=(int)(pitch/6); if(tvsm<1)tvsm=1; if(tvsm>10)tvsm=10;   // radial smoothing radius for the TV signal
  // windingDistance on RAW CT: integrate the gradient of the radially-SMOOTHED profile (removes
  // within-sheet texture, leaves pitch-scale ridge->valley swings). Calibrate per-wrap budget as
  // the median smoothed-TV over one-pitch radial segments at mid-z.
  if(wdrescue){ const f32*vsz=vs+(size_t)(dz/2)*dy*dx; int zm=dz/2; double cyz=cy[zm],cxz=cx[zm];
    double *tvs=malloc((size_t)dy*dx/16*sizeof(double)); int nt=0; int seg=(int)(pitch+0.5);
    for(int y=2;y<dy-2;y+=4)for(int x=2;x<dx-2;x+=4){ size_t p=(size_t)y*dx+x; if(v[(size_t)zm*dy*dx+p]<athr)continue;
      double ddy=y-cyz,ddx=x-cxz,r=sqrt(ddy*ddy+ddx*ddx); if(r<pitch)continue; double uy=ddy/r,ux=ddx/r;
      double prof[130]; int sp=seg<129?seg:129;
      for(int t=0;t<=sp;t++){ int iy=(int)lround(y+uy*t),ix=(int)lround(x+ux*t); prof[t]=(iy<0||iy>=dy||ix<0||ix>=dx)?0.0:(double)vsz[(size_t)iy*dx+ix]; }
      double tv=0,prev=prof[0];
      for(int t=1;t<=sp;t++){ double s=0;int c=0; for(int d=-tvsm;d<=tvsm;d++){int u=t+d;if(u<0||u>sp)continue;s+=prof[u];c++;} double sm=s/c; tv+=fabs(sm-prev); prev=sm; }
      if(nt<(int)((size_t)dy*dx/16)) tvs[nt++]=tv; }
    if(nt>0){ qsort(tvs,nt,sizeof(double),cmp_dbl); wdcal=tvs[nt/2]; }
    free(tvs);
    fprintf(stderr,"windingDistance calib: smoothed-TV/wrap=%.1f (radial-sm=%d, %d samples)\n",wdcal,tvsm,nt); }
  double wdedge=(wdrescue==1)?wdcal:0;   // per-edge rescue only in mode 1 (proven harmful); mode 2 = prior only

  // radial unit dir per voxel (in-plane, outward from per-z center)
  // ridge = radial intensity max; recto = material with air outward; bar = radial min
  u8*ridge=calloc(nn,1),*recto=calloc(nn,1),*bar=calloc(nn,1);
  #pragma omp parallel for schedule(static)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr)continue;
    double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<2)continue; double uy=ddy/r,ux=ddx/r;
    int iy=(int)lround(y-uy),ix=(int)lround(x-ux),oy=(int)lround(y+uy),ox=(int)lround(x+ux);
    if(iy<0||iy>=dy||ix<0||ix>=dx||oy<0||oy>=dy||ox<0||ox>=dx)continue;
    // THICK ridge (2D-style): a radial local max. Per-slice labeling means this can't
    // percolate across z; within a slice the thick crest is TANGENTIALLY CONNECTED into one
    // big wrap-arc node (not fragments), so the graph is connected -- the whole point of B.
    f32 c=sig[p],in=sig[IDX(z,iy,ix)],ou=sig[IDX(z,oy,ox)];
    if(c>=in&&c>=ou&&c>=ridgethr) ridge[p]=1;
    if(c<=in&&c<=ou) bar[p]=1;                            // radial local min = inter-sheet boundary
    for(int k=1;k<=3;k++){ int ry=(int)lround(y+uy*k),rx=(int)lround(x+ux*k);
      if(ry>=0&&ry<dy&&rx>=0&&rx<dx&&v[IDX(z,ry,rx)]<athr){ recto[p]=1; break; } } }
  if(!use_recto) memset(recto,0,nn);   // ridge-only mode: drop the unreliable recto graph
  // BARRIER COMPLETION: the dense fill leaks across tightly-packed wraps wherever the
  // inter-sheet wall (bar) has a gap. Per-slice in-plane morphological CLOSE (dilate r,
  // erode r) seals tangential gaps in the walls without net thickening, so isotropic
  // diffusion is actually stopped between wraps. argv[18] radius (0=off).
  int barclose=argc>18?atoi(argv[18]):1;
  if(barclose>0){ u8*tb=malloc((size_t)dy*dx);
    for(int z=0;z<dz;z++){ u8*bz=bar+(size_t)z*dy*dx;
      for(int it=0;it<barclose;it++){ memcpy(tb,bz,(size_t)dy*dx);
        for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t q=(size_t)y*dx+x; if(tb[q]||v[IDX(z,y,x)]<athr)continue;
          if((x>0&&tb[q-1])||(x<dx-1&&tb[q+1])||(y>0&&tb[q-dx])||(y<dy-1&&tb[q+dx])) bz[q]=1; } }
      for(int it=0;it<barclose;it++){ memcpy(tb,bz,(size_t)dy*dx);
        for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t q=(size_t)y*dx+x; if(!tb[q])continue;
          if((x>0&&!tb[q-1])||(x<dx-1&&!tb[q+1])||(y>0&&!tb[q-dx])||(y<dy-1&&!tb[q+dx])) bz[q]=0; } } }
    free(tb); }
  PHASE("3D detect");

  // PER-SLICE 2D labeling (NOT 3D): 3D connected components percolate across touching wraps
  // (the touching-sheet frontier) -- any bridge anywhere in z merges wraps into one node ->
  // flat winding. So label each z-slice in 2D (clean, like the 2D tool), giving per-slice
  // nodes, then tie them across z with overlap edges in the GLOBAL solve below. This yields
  // real 3D winding WITHOUT relying on 3D connectivity to separate sheets.
  u32*rl=calloc(nn,sizeof(u32)),*sl=calloc(nn,sizeof(u32)); u32 nrid=0,nrec=0;
  { u32*tl=malloc((size_t)dy*dx*sizeof(u32));
    for(int z=0;z<dz;z++){ u32 nc=cc_label(ridge+(size_t)z*dy*dx,1,dy,dx,TOPO_CONN26,tl);
      for(size_t i=0;i<(size_t)dy*dx;i++) rl[(size_t)z*dy*dx+i]= tl[i]? nrid+tl[i] : 0; nrid+=nc; }
    for(int z=0;z<dz;z++){ u32 nc=cc_label(recto+(size_t)z*dy*dx,1,dy,dx,TOPO_CONN26,tl);
      for(size_t i=0;i<(size_t)dy*dx;i++) sl[(size_t)z*dy*dx+i]= tl[i]? nrec+tl[i] : 0; nrec+=nc; }
    free(tl); }
  // sizes + mean radius
  size_t *ra=calloc((size_t)nrid+1,sizeof(size_t)),*sa=calloc((size_t)nrec+1,sizeof(size_t));
  double *rr=calloc((size_t)nrid+1,sizeof(double)),*sr=calloc((size_t)nrec+1,sizeof(double));
  double *rz=calloc((size_t)nrid+1,sizeof(double)),*sz=calloc((size_t)nrec+1,sizeof(double));  // per-node mean z (labels are per-slice)
  double *ry=calloc((size_t)nrid+1,sizeof(double)),*rx=calloc((size_t)nrid+1,sizeof(double));  // per-node centroid (y,x) for the wd-prior ray
  double *syc=calloc((size_t)nrec+1,sizeof(double)),*sxc=calloc((size_t)nrec+1,sizeof(double));
  double *rmn=malloc(((size_t)nrid+1)*sizeof(double)),*rmx=malloc(((size_t)nrid+1)*sizeof(double));
  for(u32 L=0;L<=nrid;L++){ rmn[L]=1e30; rmx[L]=-1e30; }
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x);
    double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx);
    u32 L=rl[p]; if(L){ra[L]++;rr[L]+=r;rz[L]+=z;ry[L]+=y;rx[L]+=x; if(r<rmn[L])rmn[L]=r; if(r>rmx[L])rmx[L]=r;} u32 M=sl[p]; if(M){sa[M]++;sr[M]+=r;sz[M]+=z;syc[M]+=y;sxc[M]+=x;} }
  // DIAGNOSTIC: ridge nodes whose radial extent >> pitch are MERGED multi-wrap blobs
  { int merged=0,nbig=0; size_t bigp=0; u32 bigL=0; double bigext=0;
    for(u32 L=1;L<=nrid;L++){ if(ra[L]<(size_t)minseg)continue; nbig++; double ext=rmx[L]-rmn[L];
      if(ext>1.5*pitch) merged++; if(ra[L]>bigp){bigp=ra[L];bigL=L;} if(ext>bigext)bigext=ext; }
    fprintf(stderr,"ridge nodes radial-extent: %d of %d span >1.5*pitch (MERGED); biggest=%zu vox extent=%.0f (pitch=%.0f); max extent=%.0f\n",
      merged, nbig, bigp, bigL?rmx[bigL]-rmn[bigL]:0, pitch, bigext); }
  int *rcid=malloc((size_t)(nrid+1)*sizeof(int)),*cid=malloc((size_t)(nrec+1)*sizeof(int)); int RK=0,K=0;
  for(u32 L=1;L<=nrid;L++) rcid[L]=(ra[L]>=(size_t)minseg)?RK++:-1;
  for(u32 L=1;L<=nrec;L++) cid[L]=(sa[L]>=(size_t)minseg)?K++:-1;
  if(K<2&&RK<2){ fprintf(stderr,"too few 3D nodes (recto=%d ridge=%d)\n",K,RK); return 1; }
  double *rmeanr=calloc(RK?RK:1,sizeof(double)),*meanr=calloc(K?K:1,sizeof(double));
  double *nodez=calloc((size_t)(K+RK?K+RK:1),sizeof(double));   // per-node mean z, indexed like the solve (recto 0..K-1, ridge K..)
  double *nodeY=calloc((size_t)(K+RK?K+RK:1),sizeof(double)),*nodeX=calloc((size_t)(K+RK?K+RK:1),sizeof(double));
  for(u32 L=1;L<=nrid;L++) if(rcid[L]>=0){ rmeanr[rcid[L]]=rr[L]/ra[L]; nodez[K+rcid[L]]=rz[L]/ra[L]; nodeY[K+rcid[L]]=ry[L]/ra[L]; nodeX[K+rcid[L]]=rx[L]/ra[L]; }
  for(u32 L=1;L<=nrec;L++) if(cid[L]>=0){ meanr[cid[L]]=sr[L]/sa[L]; nodez[cid[L]]=sz[L]/sa[L]; nodeY[cid[L]]=syc[L]/sa[L]; nodeX[cid[L]]=sxc[L]/sa[L]; }
  fprintf(stderr,"3D nodes: recto=%d ridge=%d (>=%d vox)\n",K,RK,minseg); PHASE("3D cc-label");

  // ---- 3D radial-adjacency graph (recto idx 0..K-1, ridge idx K..K+RK-1) ----
  int NT=K+RK;
  long cap1=1<<16,nr1=0; long*raw1=malloc(cap1*sizeof(long));
  long cap0=1<<16,nr0=0; long*raw0=malloc(cap0*sizeof(long));
  #define PUSH(arr,n,cap,key) do{ if((n)>=(cap)){(cap)*=2;(arr)=realloc((arr),(cap)*sizeof(long));} (arr)[(n)++]=(key); }while(0)
  // Walk outward, but DON'T stop at the first different node -- per-slice nodes are
  // fragmented, so the first hit is usually a SAME-WRAP fragment (valley-count 0). Keep
  // walking: nv==0 hits become same-wrap TIES (reconnect the wrap's fragments); stop at
  // the first hit that has actually CROSSED a valley (nv>=1) = the next wrap (+nv edge).
  #define RADWALK(MARK,LBL,CID,BASE) \
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(!(MARK)[p])continue; \
      int La=(CID)[(LBL)[p]]; if(La<0)continue; double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<1)continue; \
      double uy=ddy/r,ux=ddx/r; const f32*vsz=sig+(size_t)z*dy*dx; int prevL=La; \
      for(int k=2;k<walkr;k++){ int ny=(int)lround(y+uy*k),nx=(int)lround(x+ux*k); if(ny<0||ny>=dy||nx<0||nx>=dx)break; size_t q=IDX(z,ny,nx); \
        if((MARK)[q]){ int Lb=(CID)[(LBL)[q]]; if(Lb>=0&&Lb!=prevL&&k>=kmin){ int nv=count_valleys(vsz,dy,dx,y,x,uy,ux,k,sigprom,pitchmin,wdedge,tvsm); if(nv>15)nv=15; \
          if(nv==0){ prevL=Lb; continue; }   /* same-wrap fragment: walk THROUGH (no tie -- ties at touches are unreliable) */ \
          PUSH(raw1,nr1,cap1,(((long)((BASE)+La)*NT+(BASE)+Lb)<<4)|nv); break; } } } }
  RADWALK(recto,sl,cid,0);
  RADWALK(ridge,rl,rcid,K);
  #undef RADWALK
  // same-wrap ties: ridge centerline -> recto rim of the SAME sheet, a SHORT outward walk
  // (the recto edge is only ~half a sheet thickness out of the centerline). tier over-reaches
  // to the NEXT wrap's recto and flattens; cap at ~0.3*pitch so the tie is same-sheet-reliable.
  int rrtie=(int)(0.3*pitch); if(rrtie<2)rrtie=2; if(rrtie>6)rrtie=6; (void)tier;
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(!ridge[p])continue;
    int Ld=rcid[rl[p]]; if(Ld<0)continue; double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<1)continue; double uy=ddy/r,ux=ddx/r;
    for(int k=1;k<=rrtie;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k); if(oy<0||oy>=dy||ox<0||ox>=dx)break; size_t q=IDX(z,oy,ox);
      if(recto[q]){ int Lr=cid[sl[q]]; if(Lr>=0) PUSH(raw0,nr0,cap0,(long)(K+Ld)*NT+Lr); break; } } }
  // TANGENTIAL ties: connect ridge fragments along the TANGENT (perp to radial) = same wrap,
  // same radius, adjacent angle. These join the angular islands into per-wrap rings so the
  // graph is ONE connected component (else each radial chain floats with its own offset ->
  // huge winding spread at fixed radius). Tangential motion stays on the same wrap (wraps are
  // radially separated), so unlike radial/z ties they CANNOT bridge different wraps.
  int rtang=(int)(0.5*pitch); if(rtang<3)rtang=3; if(rtang>14)rtang=14;
  #define TANGTIE(MARK,LBL,CID,BASE) \
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(!(MARK)[p])continue; \
      int La=(CID)[(LBL)[p]]; if(La<0)continue; double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<1)continue; \
      double ty=-ddx/r,tx=ddy/r; \
      for(int sgn=-1;sgn<=1;sgn+=2)for(int k=1;k<=rtang;k++){ int ny=(int)lround(y+ty*sgn*k),nx=(int)lround(x+tx*sgn*k); \
        if(ny<0||ny>=dy||nx<0||nx>=dx)break; size_t q=IDX(z,ny,nx); \
        if((MARK)[q]){ int Lb=(CID)[(LBL)[q]]; if(Lb>=0&&Lb!=La){ PUSH(raw0,nr0,cap0,(long)((BASE)+La)*NT+(BASE)+Lb); } break; } } }
  TANGTIE(ridge,rl,rcid,K);
  if(use_recto) TANGTIE(recto,sl,cid,0);
  #undef TANGTIE
  // Z-TIE edges (THE 3D link): a per-slice node and the spatially-OVERLAPPING node in the
  // next slice are the same wrap -> step-0 tie. This propagates one global winding through
  // z WITHOUT 3D connectivity percolating across touching wraps. Overlap = shared (y,x).
  int do_ztie=argc>14?atoi(argv[14]):0;
  if(do_ztie) for(int z=0;z<dz-1;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x), pz=p+(size_t)dy*dx;
    if(ridge[p]&&ridge[pz]){ int A=rcid[rl[p]],B=rcid[rl[pz]]; if(A>=0&&B>=0&&A!=B) PUSH(raw0,nr0,cap0,(long)(K+A)*NT+(K+B)); }
    if(recto[p]&&recto[pz]){ int A=cid[sl[p]],B=cid[sl[pz]]; if(A>=0&&B>=0&&A!=B) PUSH(raw0,nr0,cap0,(long)A*NT+B); } }
  #undef PUSH
  int *a1,*b1,*a0,*b0; long *w1,*w0,*s1;
  long n1=agg_edges(raw1,nr1,NT,&a1,&b1,&w1,&s1,3);
  long n0=agg_edges(raw0,nr0,NT,&a0,&b0,&w0,NULL,2); free(raw1);free(raw0);
  { long sh[5]={0}; for(long e=0;e<n1;e++){ int s=s1[e]; if(s>4)s=4; sh[s]++; }
    fprintf(stderr,"3D graph: %ld step-edges, %ld ties; step hist [0..4+]: %ld %ld %ld %ld %ld\n",n1,n0,sh[0],sh[1],sh[2],sh[3],sh[4]); }
  PHASE("3D graph");
  double WEQ_ARG=argc>13?atof(argv[13]):3.0;

  // ---- ONE global winding solve (IRLS) with a WEAK RADIUS PRIOR ----
  // The radial step-edges give precise LOCAL structure but the graph is many disconnected
  // components (per-slice, per-angular-sector); each floats with an arbitrary integer
  // offset, and pooling them (no reliable z-tie) cancels the radial gradient. A weak pull
  // toward winding ~= (radius - r_inner)/pitch fixes each component's offset from its radius
  // (reliable) WITHOUT the wrap-bridging z/same-wrap ties. Step-edges still dominate locally.
  double rmin=1e30; for(int i=0;i<K;i++) if(meanr[i]<rmin)rmin=meanr[i]; for(int i=0;i<RK;i++) if(rmeanr[i]<rmin)rmin=rmeanr[i];
  // per-node radius indexed like the solve (recto 0..K-1, ridge K..) -- pairs with nodez[]
  double *nodeR=malloc(NT*sizeof(double));
  for(int i=0;i<K;i++) nodeR[i]=meanr[i]; for(int i=0;i<RK;i++) nodeR[K+i]=rmeanr[i];
  double *prior=malloc(NT*sizeof(double));
  for(int i=0;i<NT;i++) prior[i]=(nodeR[i]-rmin)/pitch;   // initial: constant-pitch Archimedean prior
  // WINDING-DISTANCE PRIOR (wdrescue==2): replace the Archimedean prior with the cumulative count
  // of sheets crossed along the ray from the per-z center to each node (smoothed-TV / per-wrap
  // budget). Unlike (r-rmin)/pitch this COUNTS ACTUAL SHEETS -> handles variable pitch and the
  // off-center case where radius is a poor predictor of winding (the radius prior's worst regime).
  if(wdrescue==2 && wdcal>0){
    for(int i=0;i<NT;i++){ int z=(int)lround(nodez[i]); if(z<0)z=0;if(z>=dz)z=dz-1; const f32*vsz=vs+(size_t)z*dy*dx;
      double ddy=nodeY[i]-cy[z],ddx=nodeX[i]-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<2){prior[i]=0;continue;} double uy=ddy/r,ux=ddx/r;
      int R=(int)(r+0.5); double prof[700]; int sp=R<699?R:699;
      for(int t=0;t<=sp;t++){ int iy=(int)lround(cy[z]+uy*t),ix=(int)lround(cx[z]+ux*t); prof[t]=(iy<0||iy>=dy||ix<0||ix>=dx)?0.0:(double)vsz[(size_t)iy*dx+ix]; }
      double tv=0,prev=prof[0]; for(int t=1;t<=sp;t++){ double s=0;int c=0; for(int d=-tvsm;d<=tvsm;d++){int u=t+d;if(u<0||u>sp)continue;s+=prof[u];c++;} double sm=s/c; tv+=fabs(sm-prev); prev=sm; }
      prior[i]=tv/wdcal; }
    fprintf(stderr,"using windingDistance prior (cumulative sheet-crossings from center)\n");
  }
  const double LPRI=argc>16?atof(argv[16]):0.15;        // radius-prior strength
  // VARIABLE-PITCH SPIRAL FIT (adapted from Henderson spiral-fitting's per-slice affine +
  // gap-scaling field, done classically): after a winding solve, fit a SMOOTH per-z spiral
  // r ~= r0(z) + pitch(z)*w (robust Tukey line per slice, then smoothed across z), and feed
  // it back as the prior. This replaces the single global constant `pitch`/`rmin` -> models
  // umbilicus radial offset + locally-varying winding spacing that the constant prior can't.
  const int SP=argc>20?atoi(argv[20]):0;                 // # spiral-fit refits; default 0 -- per-z variable
                                                         // pitch is a no-op on L2 (pitch ~z-invariant, prior weak by design). kept for L1/full-volume.
  double *r0z=malloc(dz*sizeof(double)),*ptz=malloc(dz*sizeof(double));
  double *wind=calloc(NT,sizeof(double)),*nw=malloc(NT*sizeof(double)),*den=malloc(NT*sizeof(double));
  for(int i=0;i<NT;i++) wind[i]=prior[i];                // warm-start from the radius prior
  double *rw=malloc((size_t)(n1?n1:1)*sizeof(double));
  const double WEQ=WEQ_ARG;
  for(int outer=0;outer<=SP;outer++){
    if(outer>0){
      // robust per-z line fit  r = r0[z] + ptz[z]*w  over nodes at that z
      double *sa0=calloc(dz,sizeof(double)); int *valid=calloc(dz,sizeof(int));
      for(int z=0;z<dz;z++){ (void)sa0;
        double a=rmin,b=pitch;                            // start from global
        for(int rep=0;rep<3;rep++){ double Sww=0,Sw=0,S1=0,Swr=0,Sr=0;
          for(int i=0;i<NT;i++){ if((int)lround(nodez[i])!=z)continue; double w=wind[i],r=nodeR[i];
            double e=r-(a+b*w),s=0.5*pitch,tw=(fabs(e)<s)?(1-(e/s)*(e/s))*(1-(e/s)*(e/s)):0;  // Tukey
            Sww+=tw*w*w;Sw+=tw*w;S1+=tw;Swr+=tw*w*r;Sr+=tw*r; }
          double det=S1*Sww-Sw*Sw; if(S1>=6&&fabs(det)>1e-9){ b=(S1*Swr-Sw*Sr)/det; a=(Sr-b*Sw)/S1; valid[z]=1; } }
        if(b<0.3*pitch)b=0.3*pitch; if(b>3*pitch)b=3*pitch; r0z[z]=a; ptz[z]=b; }
      // fill invalid z from nearest valid, then smooth across z (window +/-3) -> per-slice-affine smoothness
      for(int z=0;z<dz;z++) if(!valid[z]){ int lo=z,hi=z; while(lo>=0&&!valid[lo])lo--; while(hi<dz&&!valid[hi])hi++;
        if(lo<0&&hi>=dz){r0z[z]=rmin;ptz[z]=pitch;} else if(lo<0){r0z[z]=r0z[hi];ptz[z]=ptz[hi];}
        else if(hi>=dz){r0z[z]=r0z[lo];ptz[z]=ptz[lo];} else { double f=(double)(z-lo)/(hi-lo); r0z[z]=r0z[lo]+f*(r0z[hi]-r0z[lo]); ptz[z]=ptz[lo]+f*(ptz[hi]-ptz[lo]); } }
      { double *t0=malloc(dz*sizeof(double)),*t1=malloc(dz*sizeof(double));
        for(int z=0;z<dz;z++){ double s0=0,s1=0;int n=0; for(int d=-3;d<=3;d++){int zz=z+d; if(zz<0||zz>=dz)continue; s0+=r0z[zz];s1+=ptz[zz];n++;} t0[z]=s0/n;t1[z]=s1/n; }
        memcpy(r0z,t0,dz*sizeof(double)); memcpy(ptz,t1,dz*sizeof(double)); free(t0);free(t1); }
      for(int i=0;i<NT;i++){ int z=(int)lround(nodez[i]); if(z<0)z=0;if(z>=dz)z=dz-1; prior[i]=(nodeR[i]-r0z[z])/ptz[z]; }
      // report spiral-fit residual (median |r - model| in wraps)
      double sr2=0;long ns=0; for(int i=0;i<NT;i++){ int z=(int)lround(nodez[i]);if(z<0)z=0;if(z>=dz)z=dz-1; double e=(nodeR[i]-r0z[z])/ptz[z]-wind[i]; sr2+=fabs(e);ns++; }
      double pmn=1e30,pmx=-1e30; for(int z=0;z<dz;z++){ if(ptz[z]<pmn)pmn=ptz[z]; if(ptz[z]>pmx)pmx=ptz[z]; }
      fprintf(stderr,"spiral-fit refit %d: pitch(z) %.1f..%.1f, mean|r-model|=%.3f wraps\n",outer,pmn,pmx,ns?sr2/ns:0);
      free(sa0);free(valid);
    }
    for(long e=0;e<n1;e++)rw[e]=1.0;
    for(int round=0;round<4;round++){
      for(int it=0;it<800;it++){ memset(nw,0,NT*sizeof(double)); memset(den,0,NT*sizeof(double));
        for(long e=0;e<n1;e++){ int a=a1[e],b=b1[e]; double w=w1[e]*rw[e],s=(double)s1[e]; nw[a]+=w*(wind[b]-s);den[a]+=w; nw[b]+=w*(wind[a]+s);den[b]+=w; }
        for(long e=0;e<n0;e++){ int a=a0[e],b=b0[e]; double w=WEQ*w0[e]; nw[a]+=w*wind[b];den[a]+=w; nw[b]+=w*wind[a];den[b]+=w; }
        for(int i=0;i<NT;i++){ nw[i]+=LPRI*prior[i]; den[i]+=LPRI; }   // weak (variable-pitch) radius prior
        for(int i=0;i<NT;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
      if(round<3) for(long e=0;e<n1;e++){ double res=wind[b1[e]]-wind[a1[e]]-s1[e]; double z=res/0.4; rw[e]=1.0/(1.0+z*z); }
    }
  }
  free(prior);free(nodeR);free(r0z);free(ptz);
  double wmin=1e30,wmax=-1e30,sres=0; long nin=0;
  for(int i=0;i<NT;i++){ if(wind[i]<wmin)wmin=wind[i]; if(wind[i]>wmax)wmax=wind[i]; }
  for(long e=0;e<n1;e++) if(rw[e]>0.5){ sres+=fabs(wind[b1[e]]-wind[a1[e]]-s1[e]); nin++; }
  fprintf(stderr,"GLOBAL winding: %d nodes, %.1f..%.1f (%.0f wraps), inlier-resid=%.3f\n",NT,wmin,wmax,wmax-wmin,nin?sres/nin:0);
  // GROUND-TRUTH wrap count: cast rays from the center at mid-z and count ACTUAL sheet crossings
  // (full prominence valley-count) out to the max material radius. Median over angles = true wraps
  // crossed. The solved wmax-wmin should match this; whichever prior matches is the correct scale.
  { int zm2=dz/2; const f32*vsz=vs+(size_t)zm2*dy*dx; double cyz=cy[zm2],cxz=cx[zm2];
    double Rmx=0; for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ if(v[(size_t)zm2*dy*dx+(size_t)y*dx+x]<athr)continue;
      double rr2=hypot(y-cyz,x-cxz); if(rr2>Rmx)Rmx=rr2; }
    int NA=72; double cnts[72]; int nc=0; double shc[72]; int nsh=0;
    const f32*shz = sh? sh+(size_t)zm2*dy*dx : NULL;
    for(int a=0;a<NA;a++){ double th=2*M_PI*a/NA,uy=sin(th),ux=cos(th); int Rr=(int)(Rmx+0.5); if(Rr>799)Rr=799;
      int nv=count_valleys(vsz,dy,dx,cyz,cxz,uy,ux,Rr,vprom,pitchmin,0,tvsm); if(nv>0)cnts[nc++]=nv;
      // sheetness crossings: a sheet is a sheetness PEAK -> count valleys of (255-sheetness).
      // The sheetness field has texture suppressed, so this is the clean windingDistance signal.
      if(shz){ int np2=Rr<799?Rr:799; double prof[800]; for(int t=0;t<=np2;t++){ int iy=(int)lround(cyz+uy*t),ix=(int)lround(cxz+ux*t);
          prof[t]=(iy<0||iy>=dy||ix<0||ix>=dx)?0.0:255.0-(double)shz[(size_t)iy*dx+ix]; }
        // inline prominence valley-count on the inverted sheetness (prom ~ 0.25*255)
        double prom=60,pk=1.8*prom,cmax=-1e30,cmin=1e30; int dir=0,nv2=0,st=0;
        for(int t=0;t<=np2;t++){ double val=prof[t]; if(!st){cmax=cmin=val;st=1;continue;}
          if(dir>=0){ if(val>cmax)cmax=val; if(cmax-val>=pk){dir=-1;cmin=val;} }
          if(dir<=0){ if(val<cmin)cmin=val; if(val-cmin>=prom){ if(dir==-1)nv2++; dir=1;cmax=val; } } }
        if(nv2>0)shc[nsh++]=nv2; } }
    if(nc){ qsort(cnts,nc,sizeof(double),cmp_dbl); double shm=0; if(nsh){qsort(shc,nsh,sizeof(double),cmp_dbl);shm=shc[nsh/2];}
      if(shz) fprintf(stderr,"GROUND-TRUTH radial crossings (median/%d rays): raw-valley=%.0f, SHEETNESS=%.0f, solved=%.0f\n",nc,cnts[nc/2],shm,wmax-wmin);
      else    fprintf(stderr,"GROUND-TRUTH radial crossings (median/%d rays): raw-valley=%.0f, solved=%.0f\n",nc,cnts[nc/2],wmax-wmin); } }
  // does NODE winding climb with NODE radius? (isolates graph vs dense). bin by mean radius.
  { int NB=12; double bsum[12]={0};long bc[12]={0}; double Rmax=0;
    for(int i=0;i<K;i++) if(meanr[i]>Rmax)Rmax=meanr[i]; for(int i=0;i<RK;i++) if(rmeanr[i]>Rmax)Rmax=rmeanr[i];
    for(int i=0;i<K;i++){ int b=(int)(meanr[i]/Rmax*NB); if(b<0)b=0;if(b>=NB)b=NB-1; bsum[b]+=wind[i];bc[b]++; }
    for(int i=0;i<RK;i++){ int b=(int)(rmeanr[i]/Rmax*NB); if(b<0)b=0;if(b>=NB)b=NB-1; bsum[b]+=wind[K+i];bc[b]++; }
    fprintf(stderr,"NODE winding vs radius (inner->outer):"); for(int b=0;b<NB;b++) fprintf(stderr," %.1f",bc[b]?bsum[b]/bc[b]:0); fprintf(stderr,"\n"); }
  PHASE("3D solve");

  // ---- dense anisotropic 3D propagation: radial-blocked in-plane + z-coupled, soft anchors ----
  f32*wd=calloc(nn,sizeof(f32)); u8*anc=calloc(nn,1); f32*ancw=calloc(nn,sizeof(f32));
  for(size_t p=0;p<nn;p++){ int idx=-1; u32 Lr=sl[p],Ld=rl[p];
    if(Lr&&cid[Lr]>=0) idx=cid[Lr]; else if(Ld&&rcid[Ld]>=0) idx=K+rcid[Ld];
    if(idx>=0){ wd[p]=(f32)wind[idx]; ancw[p]=(f32)wind[idx]; anc[p]=1; } }
  // DIAGNOSTIC: dump the raw NODE winding at mid-z (anchors only, BEFORE dense) -- isolates
  // graph quality from dense-propagation quality.
  { int zmm=dz/2; size_t zbb=(size_t)zmm*dy*dx; u8*rg=malloc((size_t)dy*dx*3);
    double lo=wmin,hi=wmax,spn=(hi>lo)?hi-lo:1;
    for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=zbb+(size_t)y*dx+x,o=((size_t)y*dx+x)*3; int g=v[p]/6; rg[o]=rg[o+1]=rg[o+2]=(u8)g;
      if(!anc[p])continue; double t=(wd[p]-lo)/spn; if(t<0)t=0;if(t>1)t=1; double ss=t*4;int sg=(int)ss;double fr=ss-sg;u8 R,G,B;
      switch(sg){case 0:R=0;G=(u8)(255*fr);B=255;break;case 1:R=0;G=255;B=(u8)(255*(1-fr));break;case 2:R=(u8)(255*fr);G=255;B=0;break;case 3:R=255;G=(u8)(255*(1-fr));B=0;break;default:R=255;G=0;B=0;}
      rg[o]=R;rg[o+1]=G;rg[o+2]=B; }
    char fn[700]; snprintf(fn,sizeof fn,"%s_nodes.ppm",base); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rg,1,(size_t)dy*dx*3,f); fclose(f);} free(rg); fprintf(stderr,"wrote %s_nodes.ppm (graph node winding, mid-z)\n",base); }
  #define MAT(q) (v[q]>=athr && !bar[q])
  // in-plane coupling is HARD radial-blocked: weight = (1-(e.radial)^2)^4 so a radially-
  // aligned neighbour gets ~0 (winding must NOT diffuse across wraps). A small floor would
  // leak across many wraps over the iterations and collapse the radial gradient to the
  // anchor mean. Anchors held firmly (LAM) so each sheet keeps its solved winding.
  // ISOTROPIC diffusion, blocked ONLY at barriers (bar = inter-sheet radial minima). The
  // old hard radial tangent-gate blocked WITHIN a wrap too -> the ridge winding couldn't
  // fill the rest of its own wrap's radial band (the unanchored sheet body) -> most voxels
  // were poorly determined -> ~30% backward switches. Within a wrap, winding flows freely
  // (it's ~constant there); the barriers stop it at the inter-sheet boundary. A mild radial
  // down-weight (radw) still discourages leaking through barrier GAPS without starving the fill.
  const double omega=1.6, LAM=(argc>12?atof(argv[12]):0.6), radw=(argc>17?atof(argv[17]):1.0);
  #define INPW(comp) (1.0 - (1.0-radw)*(comp)*(comp))   // 1 tangential, radw radial
  // ROBUST dense update (adapted from VC3D GrowSurface robust_affine / drop-worst-neighbour): a
  // barrier GAP lets a wrong-wrap neighbour (~1 wrap off) leak into the mean -> backward switch.
  // Take the neighbour MEDIAN as consensus and Tukey-downweight neighbours far from it. Safe:
  // winding has NO legitimate sharp gradient tangentially or in z (only radially, barrier-blocked),
  // so a strongly-disagreeing neighbour is a leak, not signal. argv[23]=robustdiff sigma (0=off).
  // Default: ON (sigma 0.3) in the tight-pitch COMPRESSED CORE where touches/leaks dominate
  // (backward-switch -43..-60%, z-coherence ~3x better); OFF in the wide-pitch DELAMINATED regime
  // (pitch>=50) where sheets are sparse and the median rejects real signal (it hurt there).
  // Explicit argv[23] overrides the gate entirely (for testing).
  const double robsig=argc>23?atof(argv[23]):(pitch<50?0.3:0.0);
  for(int it=0;it<200;it++) for(int color=0;color<2;color++){
    #pragma omp parallel for schedule(static)
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ if(((x+y+z)&1)!=color)continue; size_t p=IDX(z,y,x);
      if(v[p]<athr||bar[p])continue;
      double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); double uy=r>1e-6?ddy/r:0,ux=r>1e-6?ddx/r:1;
      double nval[6],nwt[6]; int nn2=0;
      if(x>0   &&MAT(p-1)){ nval[nn2]=wd[p-1];nwt[nn2]=INPW(ux);nn2++; }
      if(x<dx-1&&MAT(p+1)){ nval[nn2]=wd[p+1];nwt[nn2]=INPW(ux);nn2++; }
      if(y>0   &&MAT(p-dx)){ nval[nn2]=wd[p-dx];nwt[nn2]=INPW(uy);nn2++; }
      if(y<dy-1&&MAT(p+dx)){ nval[nn2]=wd[p+dx];nwt[nn2]=INPW(uy);nn2++; }
      if(z>0   &&MAT(p-(size_t)dy*dx)){ nval[nn2]=wd[p-(size_t)dy*dx];nwt[nn2]=1.0;nn2++; }
      if(z<dz-1&&MAT(p+(size_t)dy*dx)){ nval[nn2]=wd[p+(size_t)dy*dx];nwt[nn2]=1.0;nn2++; }
      double s=0,ws=0;
      if(robsig>0 && nn2>=3){
        double sv[6]; for(int j=0;j<nn2;j++)sv[j]=nval[j];   // median of neighbour values = consensus
        for(int i=1;i<nn2;i++){ double t=sv[i];int j=i-1; while(j>=0&&sv[j]>t){sv[j+1]=sv[j];j--;} sv[j+1]=t; }
        double cons=(nn2&1)?sv[nn2/2]:0.5*(sv[nn2/2-1]+sv[nn2/2]);
        for(int j=0;j<nn2;j++){ double d=(nval[j]-cons)/robsig, rob=(d*d<1)?(1-d*d)*(1-d*d):0; double w=nwt[j]*rob; s+=w*nval[j]; ws+=w; }
      } else { for(int j=0;j<nn2;j++){ s+=nwt[j]*nval[j]; ws+=nwt[j]; } }
      if(anc[p]){ double wa=LAM*ws+1e-6; s+=wa*ancw[p]; ws+=wa; }
      if(ws>1e-9){ double tgt=s/ws; wd[p]=(f32)(wd[p]+omega*(tgt-wd[p])); } }
  }
  #undef MAT
  #undef INPW
  #undef TANW
  // POST-SOLVE Z-MEDIAN: firmly-held per-slice anchors bake a whole-slice integer offset that
  // the local z-diffusion can't overcome (it only blends boundaries). A slice offset by +1 sits
  // between z-neighbours at 0, so a robust z-median over +/-zmedw removes the integer-jump
  // outlier (wind_couple's winning move: median voting beats linear coupling) while a monotonic
  // real z-tilt passes through the median unchanged. Off-material voxels excluded from the window.
  int zmedw=argc>19?atoi(argv[19]):0;   // off: field is smooth in z (jump<<1); enable only for noisy crops
  if(zmedw>0 && dz>2*zmedw){
    f32 *wm=malloc((size_t)dz*dy*dx*sizeof(f32));
    #pragma omp parallel for schedule(static)
    for(int z=0;z<dz;z++){ f32 win[33]; for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x);
      if(v[p]<athr||bar[p]){ wm[p]=wd[p]; continue; }
      int n=0; for(int dzz=-zmedw;dzz<=zmedw;dzz++){ int zz=z+dzz; if(zz<0||zz>=dz)continue; size_t q=IDX(zz,y,x); if(v[q]<athr||bar[q])continue; win[n++]=wd[q]; }
      if(n<3){ wm[p]=wd[p]; continue; }
      for(int i=1;i<n;i++){ f32 t=win[i]; int j=i-1; while(j>=0&&win[j]>t){win[j+1]=win[j];j--;} win[j+1]=t; }
      wm[p]=win[n/2]; } }
    memcpy(wd,wm,(size_t)dz*dy*dx*sizeof(f32)); free(wm);
    fprintf(stderr,"post-solve z-median (+/-%d slices) applied\n",zmedw);
  }
  PHASE("3D dense");

  // ---- outputs: mid-z xy slice (banded), (z,radius) reslice, raw volume ----
  int zm=dz/2; size_t zb=(size_t)zm*dy*dx;
  double sp=(wmax>wmin)?wmax-wmin:1.0;
  { u8*rgb=malloc((size_t)dy*dx*3);
    for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=zb+(size_t)y*dx+x; int g=v[p]/5; size_t o=((size_t)y*dx+x)*3; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
      if(v[p]<athr||bar[p])continue; double t=(wd[p]-wmin)/sp; if(t<0)t=0;if(t>1)t=1;
      double ss=t*4; int seg=(int)ss; double fr=ss-seg; u8 R,G,B;
      switch(seg){case 0:R=0;G=(u8)(255*fr);B=255;break;case 1:R=0;G=255;B=(u8)(255*(1-fr));break;
        case 2:R=(u8)(255*fr);G=255;B=0;break;case 3:R=255;G=(u8)(255*(1-fr));B=0;break;default:R=255;G=0;B=0;}
      rgb[o]=R;rgb[o+1]=G;rgb[o+2]=B; }
    char fn[700]; snprintf(fn,sizeof fn,"%s_midz.ppm",base); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); } free(rgb);
    fprintf(stderr,"wrote %s (mid-z xy winding)\n",fn); }
  // (z,radius) reslice at fixed angle through the axis -> smooth vertical colour = z-coherent
  { double ang=0.6,uy=sin(ang),ux=cos(ang); int R=(int)(0.45*(dy<dx?dy:dx)); u8*img=calloc((size_t)dz*R*3,1);
    for(int z=0;z<dz;z++)for(int k=0;k<R;k++){ int yy=(int)lround(cy[z]+uy*k),xx=(int)lround(cx[z]+ux*k);
      if(yy<0||yy>=dy||xx<0||xx>=dx)continue; size_t p=IDX(z,yy,xx); if(v[p]<athr||bar[p])continue;
      double t=(wd[p]-wmin)/sp; if(t<0)t=0;if(t>1)t=1; size_t o=((size_t)z*R+k)*3; img[o]=(u8)(255*t); img[o+2]=(u8)(255*(1-t)); }
    char fn[700]; snprintf(fn,sizeof fn,"%s_reslice.ppm",base); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P6\n%d %d\n255\n",R,dz); fwrite(img,1,(size_t)dz*R*3,f); fclose(f); } free(img);
    fprintf(stderr,"wrote %s_reslice.ppm (rows=z cols=radius; smooth vertical=coherent)\n",base); }
  // VOXEL-FAIR switch metric: each material voxel counted ONCE (the per-ray version
  // over-samples the umbilicus-centered core ~360x and inflates it). Outward radial
  // neighbour should have HIGHER winding; a drop is a backward sheet switch. Split
  // inner-third vs outer to localise.
  { long sw=0,matc=0,sw_in=0,mat_in=0;
    double Rmax=0; for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr||bar[p])continue;
      double r=hypot(y-cy[z],x-cx[z]); if(r>Rmax)Rmax=r; }
    double rcore=Rmax*0.34;
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr||bar[p])continue;
      double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<2)continue; double uy=ddy/r,ux=ddx/r;
      int oy=(int)lround(y+uy*2),ox=(int)lround(x+ux*2); if(oy<0||oy>=dy||ox<0||ox>=dx)continue; size_t o=IDX(z,oy,ox);
      if(v[o]<athr||bar[o])continue; matc++; int inr=(r<rcore); if(inr)mat_in++;
      if((double)wd[o]-(double)wd[p]<-0.3){ sw++; if(inr)sw_in++; } }
    printf("3D winding: %.0f wraps, inlier-resid=%.3f, backward-switch=%.1f/1000 voxels (inner %.1f / outer %.1f)\n",
      wmax-wmin,nin?sres/nin:0,1000.0*sw/(matc?matc:1),1000.0*sw_in/(mat_in?mat_in:1),1000.0*(sw-sw_in)/((matc-mat_in)?(matc-mat_in):1)); }
  // radial winding profile at mid-z (does winding actually climb with radius?)
  { int R=(int)(0.45*(dy<dx?dy:dx)),NB=12; double bsum[12]={0};long bc[12]={0};
    for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=zb+(size_t)y*dx+x; if(v[p]<athr||bar[p])continue;
      double ddy=y-cy[zm],ddx=x-cx[zm],r=sqrt(ddy*ddy+ddx*ddx); int b=(int)(r/R*NB); if(b<0||b>=NB)continue; bsum[b]+=wd[p];bc[b]++; }
    fprintf(stderr,"radial winding profile (mid-z, r-bin -> mean wind):");
    for(int b=0;b<NB;b++) fprintf(stderr," %.1f",bc[b]?bsum[b]/bc[b]:0); fprintf(stderr,"\n"); }
  // Z-COHERENCE (the whole point of a 3D solve): the radial winding profile should be the
  // SAME at every z. Compute the per-z profile and report the std across z per radius bin.
  // Low = the winding number is consistent through the volume (a coherent 3D field).
  { int R=(int)(0.45*(dy<dx?dy:dx)),NB=12;
    double *pm=calloc((size_t)dz*NB,sizeof(double)); long *pc=calloc((size_t)dz*NB,sizeof(long));
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr||bar[p])continue;
      double r=hypot(y-cy[z],x-cx[z]); int b=(int)(r/R*NB); if(b<0||b>=NB)continue; pm[(size_t)z*NB+b]+=wd[p]; pc[(size_t)z*NB+b]++; }
    for(int z=0;z<dz;z++)for(int b=0;b<NB;b++) if(pc[(size_t)z*NB+b]) pm[(size_t)z*NB+b]/=pc[(size_t)z*NB+b];
    double tot=0,tdet=0,tjump=0; int nbn=0;
    for(int b=0;b<NB;b++){ double m=0,m2=0; int n=0; for(int z=0;z<dz;z++) if(pc[(size_t)z*NB+b]){ double q=pm[(size_t)z*NB+b]; m+=q; m2+=q*q; n++; }
      if(n>1){ m/=n; double var=m2/n-m*m; tot+=sqrt(var>0?var:0); nbn++;
        // DETREND (quadratic): fit wind ~= a + c*z + e*z^2 and report residual std. A smooth
        // z-drift OR curve (real scroll-axis tilt/curvature, NOT an error) fits the parabola;
        // only non-smooth jumps (solver defects) survive. JUMP = mean |profile(z+1)-profile(z)|:
        // tiny => smooth (geometry); large => integer hops (solver). detrended+jump both small
        // => the apparent incoherence is real geometry and the 3D solve is z-consistent.
        double S[3][3]={{0}},rhs[3]={0};
        for(int z=0;z<dz;z++) if(pc[(size_t)z*NB+b]){ double zz=z,b0=1,b1=zz,b2=zz*zz,q=pm[(size_t)z*NB+b];
          double bb[3]={b0,b1,b2}; for(int i=0;i<3;i++){ for(int j=0;j<3;j++)S[i][j]+=bb[i]*bb[j]; rhs[i]+=bb[i]*q; } }
        double cf[3]={0}; { double A[3][4]; for(int i=0;i<3;i++){ for(int j=0;j<3;j++)A[i][j]=S[i][j]; A[i][3]=rhs[i]; }
          for(int c=0;c<3;c++){ int pv=c; for(int r=c+1;r<3;r++) if(fabs(A[r][c])>fabs(A[pv][c]))pv=r; for(int k=0;k<4;k++){double t=A[c][k];A[c][k]=A[pv][k];A[pv][k]=t;}
            if(fabs(A[c][c])<1e-9)continue; for(int r=0;r<3;r++) if(r!=c){ double f=A[r][c]/A[c][c]; for(int k=0;k<4;k++)A[r][k]-=f*A[c][k]; } }
          for(int i=0;i<3;i++) cf[i]=fabs(A[i][i])>1e-9?A[i][3]/A[i][i]:0; }
        double rss=0; for(int z=0;z<dz;z++) if(pc[(size_t)z*NB+b]){ double zz=z,fit=cf[0]+cf[1]*zz+cf[2]*zz*zz,e=pm[(size_t)z*NB+b]-fit; rss+=e*e; }
        tdet+=sqrt(rss/n);
        double js=0; int jn=0,pz=-1; for(int z=0;z<dz;z++) if(pc[(size_t)z*NB+b]){ if(pz>=0){ js+=fabs(pm[(size_t)z*NB+b]-pm[(size_t)pz*NB+b]); jn++; } pz=z; }
        if(jn) tjump+=js/jn; } }
    printf("3D z-coherence: std across z = %.2f wraps (raw), %.2f (quad-detrended=solver jitter), adj-z jump %.3f wraps/slice (small=smooth real tilt)\n",
      nbn?tot/nbn:0, nbn?tdet/nbn:0, nbn?tjump/nbn:0);
    free(pm);free(pc); }
  // winding VOLUME dump for downstream use (header {dz,dy,dx,z0,y0,x0} + dz*dy*dx float, NAN off-material)
  { char fn[700]; snprintf(fn,sizeof fn,"%s_vol.f32",base); FILE*f=fopen(fn,"wb");
    if(f){ int hdr[6]={dz,dy,dx,(int)z0,(int)y0,(int)x0}; fwrite(hdr,sizeof(int),6,f);
      float*wo=malloc(nn*sizeof(float)); for(size_t p=0;p<nn;p++) wo[p]=(v[p]>=athr&&!bar[p])?wd[p]:NAN;
      fwrite(wo,sizeof(float),nn,f); free(wo); fclose(f); fprintf(stderr,"wrote %s_vol.f32 (winding volume %dx%dx%d)\n",base,dz,dy,dx); } }
  // UNROLL: the payoff. winding w advances by exactly 1 per revolution, so it IS the flattened
  // arc-length coordinate. Lay the scroll out flat: horizontal u = winding (SAMP cols/wrap sweeps
  // azimuth, one revolution per wrap-width), vertical = z, pixel = mean source intensity of the
  // material at that (winding,z). A clean layered papyrus image validates the whole chain.
  { const int SAMP=120; int UW=(int)((wmax-wmin)*SAMP)+1; if(UW<1)UW=1; if(UW>20000)UW=20000;
    double *acc=calloc((size_t)dz*UW,sizeof(double)); long *cnt=calloc((size_t)dz*UW,sizeof(long));
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr||bar[p])continue;
      int u=(int)((wd[p]-wmin)/(wmax-wmin>1e-9?wmax-wmin:1)*(UW-1)); if(u<0||u>=UW)continue;
      acc[(size_t)z*UW+u]+=v[p]; cnt[(size_t)z*UW+u]++; }
    u8*img=malloc((size_t)dz*UW); for(size_t i=0;i<(size_t)dz*UW;i++) img[i]=cnt[i]?(u8)(acc[i]/cnt[i]):0;
    char fn[700]; snprintf(fn,sizeof fn,"%s_unroll.pgm",base); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P5\n%d %d\n255\n",UW,dz); fwrite(img,1,(size_t)dz*UW,f); fclose(f);
      fprintf(stderr,"wrote %s_unroll.pgm (unrolled scroll %dx%d, %d cols/wrap; horiz=winding, vert=z)\n",base,UW,dz,SAMP); }
    free(acc);free(cnt);free(img); }
  // ARC-LENGTH UNROLL (adapted from spiral-fitting's area-correct flatten): mapping horizontal
  // linearly to winding w stretches every wrap to the same width, so inner wraps (small radius,
  // small circumference) are horizontally stretched ~N x vs outer wraps -> text distorted. The
  // fitted Archimedean spiral gives r(w) ~= rmin + pitch*w, so physical arc length traversed is
  // S(w) = integral 2*pi*r dw = 2*pi*(rmin*w + pitch*w^2/2). Map horizontal = S(w): each wrap
  // gets pixels proportional to its true circumference -> aspect-correct, undistorted layers.
  { const double Sw0=rmin*wmin+0.5*pitch*wmin*wmin, Sw1=rmin*wmax+0.5*pitch*wmax*wmax;
    double Sspan=Sw1-Sw0; if(Sspan<1e-6)Sspan=1;
    // resolution: keep the OUTER (densest) wrap at ~SAMP cols; dS/dw at wmax = rmin+pitch*wmax
    const int SAMP=120; double dSmax=rmin+pitch*wmax; int UW=(int)(Sspan/(dSmax>1e-6?dSmax:1)*SAMP)+1;
    if(UW<1)UW=1; if(UW>30000)UW=30000;
    double *acc=calloc((size_t)dz*UW,sizeof(double)); long *cnt=calloc((size_t)dz*UW,sizeof(long));
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr||bar[p])continue;
      double w=wd[p],Sw=rmin*w+0.5*pitch*w*w; int u=(int)((Sw-Sw0)/Sspan*(UW-1)); if(u<0||u>=UW)continue;
      acc[(size_t)z*UW+u]+=v[p]; cnt[(size_t)z*UW+u]++; }
    u8*img=malloc((size_t)dz*UW); for(size_t i=0;i<(size_t)dz*UW;i++) img[i]=cnt[i]?(u8)(acc[i]/cnt[i]):0;
    char fn[700]; snprintf(fn,sizeof fn,"%s_unroll_arc.pgm",base); FILE*f=fopen(fn,"wb");
    if(f){ fprintf(f,"P5\n%d %d\n255\n",UW,dz); fwrite(img,1,(size_t)dz*UW,f); fclose(f);
      fprintf(stderr,"wrote %s_unroll_arc.pgm (arc-length unroll %dx%d; horiz=physical arc length, aspect-correct)\n",base,UW,dz); }
    free(acc);free(cnt);free(img); }
  PHASE("output");
  return 0;
}
