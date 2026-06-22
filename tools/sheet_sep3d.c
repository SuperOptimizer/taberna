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
 *   5. ONE global IRLS least-squares winding solve (anchor innermost=0).
 *   6. dense anisotropic 3D propagation (radial-blocked barriers, z-coupled, soft anchors).
 *
 * usage: sheet_sep3d ARCHIVE OUTBASE lod z0 y0 x0 dz dy dx [minseg=40] [pitch=auto] [LAM=0.6] [WEQ=1] [z-tie=0]
 * STATUS: WIP -- radial winding structure is correct (span ~matches geometry); dense-field
 * switch-rate and reliable z/recto-ridge ties are the open work.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "io/mca.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"

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
static int count_valleys(const f32 *vsz, int dy, int dx, double y, double x, double uy, double ux,
                         int kstop, double prom, double pitchmin){
  if(kstop<2) return 1;
  double peakprom=1.8*prom; double pos[260]; signed char typ[260]; int ne=0;
  pos[ne]=0; typ[ne]=1; ne++;
  double cmax=-1e30,cmin=1e30,maxp=0,minp=0; int dir=0;
  for(int t=0;t<=kstop;t++){ int iy=(int)lround(y+uy*t),ix=(int)lround(x+ux*t);
    double val=(iy<0||iy>=dy||ix<0||ix>=dx)?0.0:(double)vsz[(size_t)iy*dx+ix];
    if(t==0){cmax=cmin=val;maxp=minp=0;continue;}
    if(dir>=0){ if(val>cmax){cmax=val;maxp=t;} if(cmax-val>=peakprom){ if(ne<258&&typ[ne-1]!=1){pos[ne]=maxp;typ[ne]=1;ne++;} dir=-1;cmin=val;minp=t; } }
    if(dir<=0){ if(val<cmin){cmin=val;minp=t;} if(val-cmin>=prom){ if(ne<258&&typ[ne-1]!=-1){pos[ne]=minp;typ[ne]=-1;ne++;} dir=1;cmax=val;maxp=t; } } }
  if(typ[ne-1]!=1){ pos[ne]=kstop; typ[ne]=1; ne++; }
  int nv=0; for(int i=1;i<ne-1;i++) if(typ[i]==-1 && pos[i+1]-pos[i-1]>=pitchmin) nv++;
  return nv;
}

static int cmp_long(const void*a,const void*b){ long x=*(const long*)a,y=*(const long*)b; return (x>y)-(x<y); }
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

  // radial unit dir per voxel (in-plane, outward from per-z center)
  // ridge = radial intensity max; recto = material with air outward; bar = radial min
  u8*ridge=calloc(nn,1),*recto=calloc(nn,1),*bar=calloc(nn,1);
  #pragma omp parallel for schedule(static)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(v[p]<athr)continue;
    double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<2)continue; double uy=ddy/r,ux=ddx/r;
    int iy=(int)lround(y-uy),ix=(int)lround(x-ux),oy=(int)lround(y+uy),ox=(int)lround(x+ux);
    if(iy<0||iy>=dy||ix<0||ix>=dx||oy<0||oy>=dy||ox<0||ox>=dx)continue;
    // STRICT radial peak (margin) so ridges are THIN -- thick plateaus 26-connect adjacent
    // wraps into one giant percolated 3D blob (one node -> flat winding). Check +/-1 and +/-2.
    int iy2=(int)lround(y-uy*2),ix2=(int)lround(x-ux*2),oy2=(int)lround(y+uy*2),ox2=(int)lround(x+ux*2);
    f32 c=vs[p],in=vs[IDX(z,iy,ix)],ou=vs[IDX(z,oy,ox)];
    f32 in2=(iy2>=0&&iy2<dy&&ix2>=0&&ix2<dx)?vs[IDX(z,iy2,ix2)]:in, ou2=(oy2>=0&&oy2<dy&&ox2>=0&&ox2<dx)?vs[IDX(z,oy2,ox2)]:ou;
    if(c>in&&c>ou&&c>=in2&&c>=ou2&&c>=athr) ridge[p]=1;   // strict local max over +/-2 radial
    if(c<in&&c<ou&&c<=in2&&c<=ou2) bar[p]=1;              // strict local min
    for(int k=1;k<=3;k++){ int ry=(int)lround(y+uy*k),rx=(int)lround(x+ux*k);
      if(ry>=0&&ry<dy&&rx>=0&&rx<dx&&v[IDX(z,ry,rx)]<athr){ recto[p]=1; break; } } }
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
  double *rmn=malloc(((size_t)nrid+1)*sizeof(double)),*rmx=malloc(((size_t)nrid+1)*sizeof(double));
  for(u32 L=0;L<=nrid;L++){ rmn[L]=1e30; rmx[L]=-1e30; }
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x);
    double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx);
    u32 L=rl[p]; if(L){ra[L]++;rr[L]+=r; if(r<rmn[L])rmn[L]=r; if(r>rmx[L])rmx[L]=r;} u32 M=sl[p]; if(M){sa[M]++;sr[M]+=r;} }
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
  for(u32 L=1;L<=nrid;L++) if(rcid[L]>=0) rmeanr[rcid[L]]=rr[L]/ra[L];
  for(u32 L=1;L<=nrec;L++) if(cid[L]>=0) meanr[cid[L]]=sr[L]/sa[L];
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
      double uy=ddy/r,ux=ddx/r; const f32*vsz=vs+(size_t)z*dy*dx; int prevL=La; \
      for(int k=2;k<walkr;k++){ int ny=(int)lround(y+uy*k),nx=(int)lround(x+ux*k); if(ny<0||ny>=dy||nx<0||nx>=dx)break; size_t q=IDX(z,ny,nx); \
        if((MARK)[q]){ int Lb=(CID)[(LBL)[q]]; if(Lb>=0&&Lb!=prevL&&k>=kmin){ int nv=count_valleys(vsz,dy,dx,y,x,uy,ux,k,vprom,pitchmin); if(nv>15)nv=15; \
          if(nv==0){ PUSH(raw0,nr0,cap0,(long)((BASE)+La)*NT+(BASE)+Lb); prevL=Lb; continue; } \
          PUSH(raw1,nr1,cap1,(((long)((BASE)+La)*NT+(BASE)+Lb)<<4)|nv); break; } } } }
  RADWALK(recto,sl,cid,0);
  RADWALK(ridge,rl,rcid,K);
  #undef RADWALK
  // same-wrap ties: ridge -> first recto outward (step 0)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=IDX(z,y,x); if(!ridge[p])continue;
    int Ld=rcid[rl[p]]; if(Ld<0)continue; double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); if(r<1)continue; double uy=ddy/r,ux=ddx/r;
    for(int k=1;k<=tier;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k); if(oy<0||oy>=dy||ox<0||ox>=dx)break; size_t q=IDX(z,oy,ox);
      if(recto[q]){ int Lr=cid[sl[q]]; if(Lr>=0) PUSH(raw0,nr0,cap0,(long)(K+Ld)*NT+Lr); break; } } }
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
  double WEQ_ARG=argc>13?atof(argv[13]):1.0;

  // ---- ONE global winding solve (IRLS), anchor innermost recto (or ridge) = 0 ----
  int inner=0; double bestr=1e30;
  for(int i=0;i<K;i++) if(meanr[i]<bestr){bestr=meanr[i];inner=i;}
  if(K<1){ for(int i=0;i<RK;i++) if(rmeanr[i]<bestr){bestr=rmeanr[i];inner=K+i;} }
  double *wind=calloc(NT,sizeof(double)),*nw=malloc(NT*sizeof(double)),*den=malloc(NT*sizeof(double));
  double *rw=malloc((size_t)(n1?n1:1)*sizeof(double)); for(long e=0;e<n1;e++)rw[e]=1.0;
  const double WEQ=WEQ_ARG;
  for(int round=0;round<4;round++){
    for(int it=0;it<800;it++){ memset(nw,0,NT*sizeof(double)); memset(den,0,NT*sizeof(double));
      for(long e=0;e<n1;e++){ int a=a1[e],b=b1[e]; double w=w1[e]*rw[e],s=(double)s1[e]; nw[a]+=w*(wind[b]-s);den[a]+=w; nw[b]+=w*(wind[a]+s);den[b]+=w; }
      for(long e=0;e<n0;e++){ int a=a0[e],b=b0[e]; double w=WEQ*w0[e]; nw[a]+=w*wind[b];den[a]+=w; nw[b]+=w*wind[a];den[b]+=w; }
      den[inner]+=5;
      for(int i=0;i<NT;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
    if(round<3) for(long e=0;e<n1;e++){ double res=wind[b1[e]]-wind[a1[e]]-s1[e]; double z=res/0.4; rw[e]=1.0/(1.0+z*z); }
  }
  double wmin=1e30,wmax=-1e30,sres=0; long nin=0;
  for(int i=0;i<NT;i++){ if(wind[i]<wmin)wmin=wind[i]; if(wind[i]>wmax)wmax=wind[i]; }
  for(long e=0;e<n1;e++) if(rw[e]>0.5){ sres+=fabs(wind[b1[e]]-wind[a1[e]]-s1[e]); nin++; }
  fprintf(stderr,"GLOBAL winding: %d nodes, %.1f..%.1f (%.0f wraps), inlier-resid=%.3f\n",NT,wmin,wmax,wmax-wmin,nin?sres/nin:0);
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
  #define MAT(q) (v[q]>=athr && !bar[q])
  // in-plane coupling is HARD radial-blocked: weight = (1-(e.radial)^2)^4 so a radially-
  // aligned neighbour gets ~0 (winding must NOT diffuse across wraps). A small floor would
  // leak across many wraps over the iterations and collapse the radial gradient to the
  // anchor mean. Anchors held firmly (LAM) so each sheet keeps its solved winding.
  const double omega=1.6, LAM=(argc>12?atof(argv[12]):0.6);
  #define TANW(comp) ({ double tt=1.0-(comp)*(comp); tt*tt*tt*tt; })   // sharp tangent gate
  for(int it=0;it<200;it++) for(int color=0;color<2;color++){
    #pragma omp parallel for schedule(static)
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ if(((x+y+z)&1)!=color)continue; size_t p=IDX(z,y,x);
      if(v[p]<athr||bar[p])continue;
      double ddy=y-cy[z],ddx=x-cx[z],r=sqrt(ddy*ddy+ddx*ddx); double uy=r>1e-6?ddy/r:0,ux=r>1e-6?ddx/r:1;
      double s=0,ws=0,w;
      if(x>0   &&MAT(p-1)){ w=TANW(ux); s+=w*wd[p-1];ws+=w; }
      if(x<dx-1&&MAT(p+1)){ w=TANW(ux); s+=w*wd[p+1];ws+=w; }
      if(y>0   &&MAT(p-dx)){ w=TANW(uy); s+=w*wd[p-dx];ws+=w; }
      if(y<dy-1&&MAT(p+dx)){ w=TANW(uy); s+=w*wd[p+dx];ws+=w; }
      if(z>0   &&MAT(p-(size_t)dy*dx)){ s+=1.0*wd[p-(size_t)dy*dx];ws+=1.0; }   // z-link (genuine 3D)
      if(z<dz-1&&MAT(p+(size_t)dy*dx)){ s+=1.0*wd[p+(size_t)dy*dx];ws+=1.0; }
      if(anc[p]){ double wa=LAM*ws+1e-6; s+=wa*ancw[p]; ws+=wa; }
      if(ws>1e-9){ double tgt=s/ws; wd[p]=(f32)(wd[p]+omega*(tgt-wd[p])); } }
  }
  #undef MAT
  #undef TANW
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
  // 3D coherence metric: along outward radial rays per z, fraction of monotone steps
  { long steps=0,back=0;
    for(int z=0;z<dz;z++)for(int a=0;a<360;a+=2){ double th=a*M_PI/180.0,uy=sin(th),ux=cos(th); double pw=0;int hp=0;
      for(int k=2;k<(int)(0.45*(dy<dx?dy:dx));k++){ int yy=(int)lround(cy[z]+uy*k),xx=(int)lround(cx[z]+ux*k);
        if(yy<0||yy>=dy||xx<0||xx>=dx)break; size_t p=IDX(z,yy,xx); if(v[p]<athr||bar[p])continue; double w=wd[p];
        if(hp){ steps++; if(w-pw<-0.3)back++; } pw=w;hp=1; } }
    printf("3D winding: %.0f wraps, inlier-resid=%.3f, backward-rate=%.1f/1000 radial steps\n",wmax-wmin,nin?sres/nin:0,1000.0*back/(steps?steps:1)); }
  // radial winding profile at mid-z (does winding actually climb with radius?)
  { int R=(int)(0.45*(dy<dx?dy:dx)),NB=12; double bsum[12]={0};long bc[12]={0};
    for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=zb+(size_t)y*dx+x; if(v[p]<athr||bar[p])continue;
      double ddy=y-cy[zm],ddx=x-cx[zm],r=sqrt(ddy*ddy+ddx*ddx); int b=(int)(r/R*NB); if(b<0||b>=NB)continue; bsum[b]+=wd[p];bc[b]++; }
    fprintf(stderr,"radial winding profile (mid-z, r-bin -> mean wind):");
    for(int b=0;b<NB;b++) fprintf(stderr," %.1f",bc[b]?bsum[b]/bc[b]:0); fprintf(stderr,"\n"); }
  PHASE("output");
  return 0;
}
