/* sheet_sep.c — separate touching papyrus sheets via the DELAMINATION GRAPH (no ML).
 *
 * The .mca masks air as 0. Where wraps delaminate there is real air between them; where
 * they touch there is not. This tool uses the reliable delaminated surfaces to fit the
 * global spiral and assign a consistent RELATIVE winding, even across touching regions:
 *
 *   1. recto/verso: each papyrus voxel (v!=0) with air a few voxels OUTWARD along the
 *      radius from the umbilicus is a RECTO face (outer-wrap side); air inward = VERSO.
 *   2. label the recto rims into delamination SEGMENTS (cc_label).
 *   3. radial-adjacency GRAPH: from each recto voxel, walk outward to the next recto
 *      voxel of a different segment -> edge "next wrap outward" (+1).
 *   4. fit the winding: least-squares so wind[B]=wind[A]+1 over all edges (the segments
 *      lie on one spiral), solved by iterative relaxation; anchor the innermost to 0.
 *
 * Validated in prototype: edge residual ~0.06 (the delaminations ARE globally spiral-
 * consistent). Output OUT_wind.ppm: each delamination segment coloured by its winding.
 *
 * Usage: sheet_sep ARCHIVE.mca OUTBASE lod z0 y0 x0 d [minseg=15]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "io/mca.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"

static double tnow(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+ts.tv_nsec*1e-9; }
static double g_tl;
#define PHASE(name) do{ fprintf(stderr,"[t] %-18s %.2fs\n",name,tnow()-g_tl); g_tl=tnow(); }while(0)

/* RADIAL PITCH (vox/wrap): cast rays from the umbilicus, record radial positions of
 * material-intensity local maxima (sheet centerlines), and take the median gap between
 * consecutive maxima. This is the TRUE wrap pitch (it scales correctly with LOD), unlike
 * the segment-radius-gap proxy (which shrinks as segment count grows). */
static double measure_pitch(const u8 *v, int d, double cyf, double cxf, int athr){
  // AUTOCORRELATION of the radial intensity profile: raw peak-picking locks onto
  // noise (~2-3 vox), but the wrap periodicity shows as the first dominant autocorr
  // lag. Per ray: build the profile, autocorrelate, find the first prominent local-max
  // lag in [2, R/2]; median over rays = pitch. Robust + scales correctly with LOD.
  int R=(int)(0.45*d); if(R>800)R=800; if(R<16) return 8.0;
  double *prof=malloc(R*sizeof(double)), *ac=malloc(R*sizeof(double));
  double lags[256]; int nl=0;
  for(int a=0;a<180 && nl<256;a++){
    double th=a*(3.14159265/90.0), uy=sin(th), ux=cos(th);
    int n=0; for(int k=0;k<R;k++){ int yy=(int)lround(cyf+uy*k),xx=(int)lround(cxf+ux*k);
      double val = (yy<0||yy>=d||xx<0||xx>=d)?0 : v[(size_t)yy*d+xx];
      prof[k]= val<athr?0:val-athr; if(prof[k]>0)n=k; }       // profile, air->0, last material at n
    if(n<24) continue;                                         // ray too short
    double mean=0; for(int k=0;k<=n;k++) mean+=prof[k]; mean/=(n+1);
    for(int k=0;k<=n;k++) prof[k]-=mean;
    for(int L=0;L<=n/2;L++){ double s=0; for(int k=0;k+L<=n;k++) s+=prof[k]*prof[k+L]; ac[L]=s; }
    if(ac[0]<=0) continue;
    int best=-1;                                               // first prominent local max, lag>=2
    for(int L=2;L<n/2-1;L++) if(ac[L]>ac[L-1]&&ac[L]>=ac[L+1]&&ac[L]>0.15*ac[0]){ best=L; break; }
    if(best>=2 && best<200) lags[nl++]=best;
  }
  free(prof); free(ac);
  if(nl<8) return 8.0;
  for(int i=0;i<nl;i++)for(int j=i+1;j<nl;j++) if(lags[j]<lags[i]){double t=lags[i];lags[i]=lags[j];lags[j]=t;}
  return lags[nl/2];
}

/* BIMODAL air threshold: the airless archive no longer stores air as 0 (only the
 * SAM2 exterior is 0); internal delamination air is low-but-nonzero. Find the valley
 * between the air mode (near 0) and the material mode (the big hump) on the region
 * histogram. A voxel is AIR iff v < thr, MATERIAL iff v >= thr. Replaces the old
 * v==0 / v!=0 tests so the pipeline is bimodal-aware (no reliance on baked zeros). */
static int air_threshold(const u8 *v, size_t n){
  // The airless export no longer pre-zeros INTRA-scroll air (the old archive baked
  // the fysics air-cut in; we now do it post-decode). The exterior is hard-zeroed by
  // the SAM2 mask -> a huge spike at v==0 that must be EXCLUDED, else it dominates the
  // histogram and the threshold lands just above 0 (missing the real internal air,
  // which is partial-volume bright). So: histogram the IN-VOLUME values 1..254 only,
  // then Otsu (parameter-free) to cut in the air/papyrus valley. (Otsu is robust to
  // the air hump being TALLER than the papyrus hump, which it is on this data.)
  long h[256]={0}; long tot=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x>=1 && x<=254){ h[x]++; tot++; } }
  if(tot<256) return 1;
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sumB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){
    wB+=h[t]; if(wB==0) continue; double wF=(double)tot-wB; if(wF<=0) break;
    sumB+=(double)t*h[t];
    double mB=sumB/wB, mF=(sum-sumB)/wF, btw=wB*wF*(mB-mF)*(mB-mF);
    if(btw>best){ best=btw; thr=t; }
  }
  return thr;
}

/* Build the radial-adjacency graph over marker voxels (recto rims OR intensity ridges)
 * and fit a consistent winding (wind[B]=wind[A]+1 over edges), anchoring the innermost
 * segment to 0. Returns malloc'd wind[K]; sets *residual and *nedges. */
static double *graph_winding(const u8 *marker, const u32 *lbl, const int *cid, int K,
                             const double *meanr, int d, double cyf, double cxf, int wmax, int kmin,
                             double *residual, long *nedges){
  long *edge=calloc((size_t)K*K,sizeof(long));
  // A +1 edge must reach the IMMEDIATE next wrap (~1 pitch out). Different-segment
  // markers closer than kmin are same-wrap FRAGMENTS (would be a false +1); skip them
  // and keep walking. Beyond wmax would be a 2-ring SKIP; the loop ends there.
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!marker[p])continue;
    int La=cid[lbl[p]]; if(La<0)continue;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1)continue; double uy=dy/r,ux=dx/r;
    for(int k=2;k<wmax;k++){ int ny_=(int)lround(y+uy*k),nx_=(int)lround(x+ux*k);
      if(ny_<0||ny_>=d||nx_<0||nx_>=d)break; size_t q=(size_t)ny_*d+nx_;
      if(marker[q]){ int Lb=cid[lbl[q]]; if(Lb>=0&&Lb!=La&&k>=kmin){ edge[(size_t)La*K+Lb]++; break; } } } }
  double *wind=calloc(K,sizeof(double));
  int inner=0; for(int i=1;i<K;i++) if(meanr[i]<meanr[inner]) inner=i;
  // SPARSE edge list: the K*K matrix is ~99.9% zero; relaxing over the nonzero edges
  // is O(iters*E) not O(iters*K*K) (an ~800x speedup at K~1000).
  long ne=0; for(int a=0;a<K;a++)for(int b=0;b<K;b++) if(edge[(size_t)a*K+b]>=3) ne++;
  int *ea=malloc((size_t)(ne?ne:1)*sizeof(int)),*eb=malloc((size_t)(ne?ne:1)*sizeof(int));
  long *ew=malloc((size_t)(ne?ne:1)*sizeof(long)); long ei=0;
  for(int a=0;a<K;a++)for(int b=0;b<K;b++){ long w=edge[(size_t)a*K+b]; if(w>=3){ ea[ei]=a;eb[ei]=b;ew[ei]=w;ei++; } }
  free(edge);
  double *nw=malloc(K*sizeof(double)),*den=malloc(K*sizeof(double));
  double *rw=malloc((size_t)(ne?ne:1)*sizeof(double)); for(long e=0;e<ne;e++) rw[e]=1.0; // IRLS robust weights
  // IRLS: solve, then DOWN-WEIGHT edges whose residual is far from +1 (skips / fragment
  // bridges that survived the kmin filter), and re-solve. Tukey-like w = 1/(1+(res/0.4)^2).
  for(int round=0;round<4;round++){
    for(int it=0;it<600;it++){ memset(nw,0,K*sizeof(double)); memset(den,0,K*sizeof(double));
      for(long e=0;e<ne;e++){ int a=ea[e],b=eb[e]; double w=ew[e]*rw[e];
        nw[a]+=w*(wind[b]-1); den[a]+=w; nw[b]+=w*(wind[a]+1); den[b]+=w; }
      den[inner]+=5;
      for(int i=0;i<K;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
    if(round<3) for(long e=0;e<ne;e++){ double res=wind[eb[e]]-wind[ea[e]]-1; double z=res/0.4; rw[e]=1.0/(1.0+z*z); }
  }
  // report INLIER residual (robust weight > 0.5 == within ~0.4 of a clean +1 step)
  double sres=0; long nin=0; for(long e=0;e<ne;e++){ if(rw[e]>0.5){ sres+=fabs(wind[eb[e]]-wind[ea[e]]-1); nin++; } }
  free(nw);free(den);free(rw);free(ea);free(eb);free(ew);
  if(residual)*residual = nin? sres/nin : 0; if(nedges)*nedges=nin; return wind;
}

/* DELAMINATION-ANCHORED UNIFIED winding. Two node sets share one spiral:
 *   recto faces (Kr nodes, indices 0..Kr-1)  — reliable, only in delaminated regions
 *   ridge centerlines (Kd nodes, indices Kr..Kr+Kd-1) — reach the compressed core, noisy
 * Edges of two kinds:
 *   +1  : radial-adjacency between DIFFERENT wraps (recto<->recto, ridge<->ridge)
 *   ==  : SAME-WRAP tie between a ridge centerline and the recto face of the same sheet
 *         (found by a SHORT outward walk from the ridge: the outer face is ~half a sheet
 *         thickness outward of the centerline). These ties pin the ridge winding to the
 *         trustworthy recto winding, fixing the radius-offset bug that wrecked the naive
 *         merge. Anchor innermost recto = 0. Returns wind[Kr+Kd]. */
static int cmp_long(const void*a,const void*b){ long x=*(const long*)a,y=*(const long*)b; return (x>y)-(x<y); }
// aggregate a raw list of (a*K+b) edge keys into sparse arrays, keeping keys with run>=minc.
static long agg_edges(long*raw,long nraw,int K,int**ao,int**bo,long**wo,long minc){
  if(nraw==0){ *ao=malloc(sizeof(int)); *bo=malloc(sizeof(int)); *wo=malloc(sizeof(long)); return 0; }
  qsort(raw,nraw,sizeof(long),cmp_long);
  long ne=0; for(long i=0;i<nraw;){ long j=i; while(j<nraw&&raw[j]==raw[i])j++; if(j-i>=minc)ne++; i=j; }
  int*ea=malloc((size_t)(ne?ne:1)*sizeof(int)),*eb=malloc((size_t)(ne?ne:1)*sizeof(int)); long*ew=malloc((size_t)(ne?ne:1)*sizeof(long));
  long e=0; for(long i=0;i<nraw;){ long j=i; while(j<nraw&&raw[j]==raw[i])j++; long c=j-i;
    if(c>=minc){ ea[e]=(int)(raw[i]/K); eb[e]=(int)(raw[i]%K); ew[e]=c; e++; } i=j; }
  *ao=ea; *bo=eb; *wo=ew; return ne;
}
static double *unified_winding(const u8*recto,const u32*rsl,const int*rcid,int Kr,const double*rmeanr,
                               const u8*ridge,const u32*rdl,const int*rdcid,int Kd,
                               int d,double cyf,double cxf,int wmax,int tiemax,int kmin,
                               double*resid_p1,double*resid_eq,long*np1,long*neq){
  int K=Kr+Kd; if(K<2) return NULL;
  // SPARSE edge accumulation: push (a*K+b) keys into growable lists, then sort+aggregate
  // -- O(E log E) over the ~10^5 real edges, NOT O(K*K) (214M at L2). Was the 20s phase.
  long cap1=4096,nr1=0; long*raw1=malloc(cap1*sizeof(long));
  long cap0=4096,nr0=0; long*raw0=malloc(cap0*sizeof(long));
  #define PUSH1(key) do{ if(nr1>=cap1){cap1*=2;raw1=realloc(raw1,cap1*sizeof(long));} raw1[nr1++]=(key); }while(0)
  #define PUSH0(key) do{ if(nr0>=cap0){cap0*=2;raw0=realloc(raw0,cap0*sizeof(long));} raw0[nr0++]=(key); }while(0)
  #define RADWALK(MARK,LBL,CID,BASE) \
    for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!(MARK)[p])continue; \
      int La=(CID)[(LBL)[p]]; if(La<0)continue; double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); \
      if(r<1)continue; double uy=dy/r,ux=dx/r; \
      for(int k=2;k<wmax;k++){ int ny_=(int)lround(y+uy*k),nx_=(int)lround(x+ux*k); \
        if(ny_<0||ny_>=d||nx_<0||nx_>=d)break; size_t q=(size_t)ny_*d+nx_; \
        if((MARK)[q]){ int Lb=(CID)[(LBL)[q]]; if(Lb>=0&&Lb!=La&&k>=kmin){ PUSH1((long)((BASE)+La)*K+(BASE)+Lb); break; } } } }
  RADWALK(recto,rsl,rcid,0);
  RADWALK(ridge,rdl,rdcid,Kr);
  #undef RADWALK
  // same-wrap ties: from each ridge voxel walk OUTWARD; first recto hit = same wrap.
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!ridge[p])continue;
    int Ld=rdcid[rdl[p]]; if(Ld<0)continue; double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx);
    if(r<1)continue; double uy=dy/r,ux=dx/r;
    for(int k=1;k<=tiemax;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k);
      if(oy<0||oy>=d||ox<0||ox>=d)break; size_t q=(size_t)oy*d+ox;
      if(recto[q]){ int Lr=rcid[rsl[q]]; if(Lr>=0){ PUSH0((long)(Kr+Ld)*K+Lr); } break; } } }
  #undef PUSH1
  #undef PUSH0
  int inner=0; for(int i=1;i<Kr;i++) if(rmeanr[i]<rmeanr[inner]) inner=i;
  const double WEQ=3.0;
  double *wind=calloc(K,sizeof(double));
  int *a1,*b1,*a0,*b0; long *w1l,*w0l;
  long n1=agg_edges(raw1,nr1,K,&a1,&b1,&w1l,3);
  long n0=agg_edges(raw0,nr0,K,&a0,&b0,&w0l,2);
  free(raw1);free(raw0);
  double*nw=malloc(K*sizeof(double)),*den=malloc(K*sizeof(double));
  for(int it=0;it<3000;it++){ memset(nw,0,K*sizeof(double)); memset(den,0,K*sizeof(double));
    for(long e=0;e<n1;e++){ int a=a1[e],b=b1[e]; long w=w1l[e]; nw[a]+=w*(wind[b]-1); den[a]+=w; nw[b]+=w*(wind[a]+1); den[b]+=w; }
    for(long e=0;e<n0;e++){ int a=a0[e],b=b0[e]; double w=WEQ*w0l[e]; nw[a]+=w*wind[b]; den[a]+=w; nw[b]+=w*wind[a]; den[b]+=w; }
    den[inner]+=5;                               // pull anchor toward 0
    for(int i=0;i<K;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
  double s1=0,s0=0;
  for(long e=0;e<n1;e++) s1+=fabs(wind[b1[e]]-wind[a1[e]]-1);
  for(long e=0;e<n0;e++) s0+=fabs(wind[b0[e]]-wind[a0[e]]);
  free(nw);free(den);free(a1);free(b1);free(w1l);free(a0);free(b0);free(w0l);
  if(resid_p1)*resid_p1=n1?s1/n1:0; if(resid_eq)*resid_eq=n0?s0/n0:0; if(np1)*np1=n1; if(neq)*neq=n0;
  return wind;
}

int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: %s ARCHIVE OUTBASE lod z0 y0 x0 d [minseg=15]\n",argv[0]);return 2;}
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]); int d=atoi(argv[7]);
  int minseg=argc>8?atoi(argv[8]):15;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);

  // global umbilicus from a coarse LOD5 read, scaled into region coords
  int cl=5; double cs=(double)(1<<cl); int cz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,cz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse read fail\n");return 1;}
  size_t ccn=(size_t)cz*ccy*ccx; u8*ccm=malloc(ccn);
  for(size_t i=0;i<ccn;i++) ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,cz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;}
  free(ccm);
  double ls=(double)(1<<lod); f32 ucy,ucx; umbilicus_center(&umb,(f32)(cz/2),&ucy,&ucx);
  double cyf=ucy*cs/ls - y0, cxf=ucx*cs/ls - x0;
  fprintf(stderr,"umbilicus in region coords: (y=%.0f, x=%.0f)\n",cyf,cxf);

  // read the (single z) plane. Archives are NON-SQUARE (ROI trim makes ny!=nx), so read
  // the actual ny x nx extent CLAMPED to the archive and embed it in a SQUARE d x d buffer
  // (remainder padded 0 = air). d<=0 -> auto-size to cover the WHOLE plane from (y0,x0).
  int ply=fy>>lod, plx=fx>>lod;
  if(d<=0){ int sy=ply-(int)y0, sx=plx-(int)x0; d=sy>sx?sy:sx; }
  int dyr=d, dxr=d;
  if((int)y0+dyr>ply) dyr=ply-(int)y0; if(dyr<0)dyr=0;
  if((int)x0+dxr>plx) dxr=plx-(int)x0; if(dxr<0)dxr=0;
  u8*vr=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,1,dyr,dxr); if(!vr){fprintf(stderr,"region read fail\n");return 1;}
  u8*v=calloc((size_t)d*d,1);
  for(int yy=0;yy<dyr;yy++) memcpy(v+(size_t)yy*d, vr+(size_t)yy*dxr, (size_t)dxr);
  free(vr);
  size_t nn=(size_t)d*d;
  fprintf(stderr,"plane %dx%d @lod%d: read %dx%d embedded in %dx%d square\n",ply,plx,lod,dyr,dxr,d,d);
  g_tl=tnow();
  int athr=air_threshold(v,nn);              // bimodal air/material split (was v==0)
  fprintf(stderr,"bimodal air threshold: v<%d = air, v>=%d = material\n",athr,athr);

  // (1) recto faces: material with air OUTWARD along the radius
  u8*recto=calloc(nn,1);
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1)continue; double uy=dy/r,ux=dx/r;
    for(int k=1;k<=3;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k);
      if(oy>=0&&oy<d&&ox>=0&&ox<d&&v[(size_t)oy*d+ox]<athr){ recto[p]=1; break; } } }

  // (2) label recto segments, keep big ones, compact-index them
  u32*sl=calloc(nn,sizeof(u32)); u32 nseg=cc_label(recto,1,d,d,TOPO_CONN26,sl);
  size_t*sa=calloc((size_t)nseg+1,sizeof(size_t)); double*sry=calloc((size_t)nseg+1,sizeof(double));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; u32 L=sl[p]; if(!L)continue;
    sa[L]++; double dy=y-cyf,dx=x-cxf; sry[L]+=sqrt(dy*dy+dx*dx); }
  int *cid=malloc((size_t)(nseg+1)*sizeof(int)); int K=0;
  for(u32 L=1;L<=nseg;L++){ cid[L]= (sa[L]>=(size_t)minseg)? K++ : -1; }
  if(K<2){ fprintf(stderr,"too few segments (%d)\n",K); return 1; }
  double *meanr=calloc(K,sizeof(double));
  for(u32 L=1;L<=nseg;L++) if(cid[L]>=0) meanr[cid[L]]=sry[L]/sa[L];
  fprintf(stderr,"delamination segments (>=%d px): %d\n",minseg,K);

  // PITCH-AWARE walk radii: the radial-adjacency walk must reach ~1.5x the wrap
  // pitch (vox/wrap) to link a segment to its immediate radial neighbour -- a fixed
  // 40 over-reaches at coarse LODs (pitch~4) and under-reaches at fine ones.
  double pitch = measure_pitch(v,d,cyf,cxf,athr);
  int walkr=(int)(1.6*pitch+0.5); if(walkr<8)walkr=8; if(walkr>120)walkr=120;
  int kmin =(int)(0.5*pitch+0.5); if(kmin<2)kmin=2;            // reject same-wrap-fragment +1 edges
  int tier =(int)(0.7*pitch+0.5); if(tier<3)tier=3; if(tier>30)tier=30;
  fprintf(stderr,"radial pitch=%.1f vox/wrap -> walk=%d tie=%d\n",pitch,walkr,tier); PHASE("threshold+pitch");

  // (3+4) radial-adjacency graph over the recto rims + winding fit (the helper)
  double resid; long ne; double *wind=graph_winding(recto,sl,cid,K,meanr,d,cyf,cxf,walkr,kmin,&resid,&ne);
  double wmin=1e30,wmax=-1e30; for(int i=0;i<K;i++){ if(wind[i]<wmin)wmin=wind[i]; if(wind[i]>wmax)wmax=wind[i]; }
  printf("recto-graph: segments=%d edges=%ld  winding %.1f..%.1f (%.0f wraps)  edge-residual=%.3f\n",
         K,ne,wmin,wmax,wmax-wmin,resid); PHASE("recto-graph");

  // output: colour each delamination segment by its winding (hue cycle), on dim data
  u8*rgb=malloc(nn*3);
  for(size_t p=0;p<nn;p++){ int g=v[p]/4; rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=(u8)g; }
  for(size_t p=0;p<nn;p++){ u32 L=sl[p]; if(!L||cid[L]<0)continue;
    double h=fmod((wind[cid[L]]-wmin),6.0); if(h<0)h+=6; int hi=(int)h; double fr=h-hi;
    int V=255,P=30,Q=(int)(255*(1-fr)),T=(int)(255*fr); u8 R,G,B;
    switch(hi){case 0:R=V;G=T;B=P;break;case 1:R=Q;G=V;B=P;break;case 2:R=P;G=V;B=T;break;
      case 3:R=P;G=Q;B=V;break;case 4:R=T;G=P;B=V;break;default:R=V;G=P;B=Q;}
    rgb[3*p]=R;rgb[3*p+1]=G;rgb[3*p+2]=B; }
  char fn[600]; snprintf(fn,sizeof fn,"%s_wind.ppm",base);
  FILE*f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
  printf("wrote %s (delamination segments coloured by spiral-fit winding)\n",fn);

  // (5) VALLEY-EXTENSION: split partially-delaminated touches. Each delamination air
  // channel runs along a sheet boundary; extend it from its endpoints by FOLLOWING the
  // intensity valley (the faint inter-sheet minimum) until it bridges to another void.
  // Only commit cuts that reach another void (validated bridges, not arbitrary).
  f32 *vs=malloc(nn*sizeof(f32));                 // lightly smoothed intensity (air=0)
  for(size_t p=0;p<nn;p++) vs[p]=v[p];
  for(int it=0;it<2;it++){ f32*t=malloc(nn*sizeof(f32));
    for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; double s=vs[p];int c=1;
      if(x>0){s+=vs[p-1];c++;}if(x<d-1){s+=vs[p+1];c++;}if(y>0){s+=vs[p-d];c++;}if(y<d-1){s+=vs[p+d];c++;}
      t[p]=(f32)(s/c);} memcpy(vs,t,nn*sizeof(f32)); free(t);}
  for(size_t p=0;p<nn;p++) if(v[p]<athr) vs[p]=0;
  u8*air=malloc(nn); for(size_t p=0;p<nn;p++) air[p]=v[p]<athr;
  u32*al=calloc(nn,sizeof(u32)); u32 nch=cc_label(air,1,d,d,TOPO_CONN6,al);
  // per-channel moments for PCA
  size_t*ca=calloc((size_t)nch+1,sizeof(size_t)); double*my=calloc((size_t)nch+1,sizeof(double)),*mx=calloc((size_t)nch+1,sizeof(double));
  double*syy=calloc((size_t)nch+1,sizeof(double)),*sxx=calloc((size_t)nch+1,sizeof(double)),*sxy=calloc((size_t)nch+1,sizeof(double));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ u32 L=al[(size_t)y*d+x]; if(!L)continue; ca[L]++; my[L]+=y; mx[L]+=x; syy[L]+=(double)y*y; sxx[L]+=(double)x*x; sxy[L]+=(double)x*y; }
  // axis (major eigenvector) + aspect per channel
  double*axy=calloc((size_t)nch+1,sizeof(double)),*axx=calloc((size_t)nch+1,sizeof(double)); char*good=calloc((size_t)nch+1,1);
  for(u32 L=1;L<=nch;L++){ if(ca[L]<25||ca[L]>8000)continue; double N=ca[L];
    double cyy=syy[L]/N-(my[L]/N)*(my[L]/N), cxx=sxx[L]/N-(mx[L]/N)*(mx[L]/N), cxy=sxy[L]/N-(mx[L]/N)*(my[L]/N);
    double tr=cyy+cxx, det=cyy*cxx-cxy*cxy, disc=sqrt(fmax(0,tr*tr/4-det));
    double l1=tr/2+disc, l0=tr/2-disc; if(l1<1e-6)continue; if(sqrt(l1/fmax(l0,1e-6))<2.5)continue;
    // eigenvector for l1: (cxy, l1-cyy) normalized (in (x,y))
    double ex=cxy, ey=l1-cyy, n2=hypot(ex,ey); if(n2<1e-9){ex=1;ey=0;n2=1;} axx[L]=ex/n2; axy[L]=ey/n2; good[L]=1; }
  // endpoints = extreme projections along the axis
  double*pmin=malloc((size_t)(nch+1)*sizeof(double)),*pmax=malloc((size_t)(nch+1)*sizeof(double));
  int*emin=malloc((size_t)(nch+1)*sizeof(int)),*emax=malloc((size_t)(nch+1)*sizeof(int));
  for(u32 L=0;L<=nch;L++){pmin[L]=1e30;pmax[L]=-1e30;emin[L]=emax[L]=-1;}
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ u32 L=al[(size_t)y*d+x]; if(!L||!good[L])continue;
    double pr=x*axx[L]+y*axy[L]; if(pr<pmin[L]){pmin[L]=pr;emin[L]=(int)((size_t)y*d+x);} if(pr>pmax[L]){pmax[L]=pr;emax[L]=(int)((size_t)y*d+x);} }
  u8*cut=calloc(nn,1); long ncut=0,nbridge=0;
  for(u32 L=1;L<=nch;L++){ if(!good[L])continue;
    for(int s=0;s<2;s++){ int ep= s? emax[L]:emin[L]; if(ep<0)continue; double sg= s?1:-1;
      double yy=ep/d, xx=ep%d, dy=axy[L]*sg, dx=axx[L]*sg; int tmp[120]; int npath=0; int bridged=0;
      for(int t=0;t<80;t++){ double bv=1e30,by=0,bx=0,bdy=0,bdx=0; int found=0;
        for(int ai=0;ai<9;ai++){ double ang=-0.6+1.2*ai/8.0,ca2=cos(ang),sa2=sin(ang);
          double ndy=dy*ca2-dx*sa2, ndx=dy*sa2+dx*ca2, ty=yy+ndy*1.5, tx=xx+ndx*1.5;
          int iy=(int)lround(ty),ix=(int)lround(tx); if(iy<0||iy>=d||ix<0||ix>=d)continue;
          size_t q=(size_t)iy*d+ix; if(v[q]<athr){ bridged=1; goto done; }
          double val=vs[q]+fabs(ang)*8; if(val<bv){bv=val;by=ty;bx=tx;bdy=ndy;bdx=ndx;found=1;} }
        if(!found)break; yy=by;xx=bx; double m=hypot(bdy,bdx); dy=bdy/m;dx=bdx/m;
        int iy=(int)lround(yy),ix=(int)lround(xx); if(npath<120)tmp[npath++]=(int)((size_t)iy*d+ix); }
      done:; if(bridged){ for(int i=0;i<npath;i++){ if(!cut[tmp[i]]){cut[tmp[i]]=1;ncut++;} } nbridge++; } } }
  printf("valley-extension: %ld bridging cuts (%ld cut px)\n",nbridge,ncut);
  // segmentation = material minus the cuts (dilated 1)
  u8*sep=malloc(nn); for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x;
    int c=cut[p]; if(!c){ if(x>0&&cut[p-1])c=1; else if(x<d-1&&cut[p+1])c=1; else if(y>0&&cut[p-d])c=1; else if(y<d-1&&cut[p+d])c=1; }
    sep[p]= (v[p]>=athr && !c); }
  u32*pl=calloc(nn,sizeof(u32)); u32 npc=cc_label(sep,1,d,d,TOPO_CONN26,pl);
  size_t*pa=calloc((size_t)npc+1,sizeof(size_t)); for(size_t p=0;p<nn;p++) pa[pl[p]]++;
  long nbig=0; for(u32 L=1;L<=npc;L++) if(pa[L]>=150)nbig++;
  printf("sheet pieces after valley-cut: %ld (>=150px)\n",nbig);
  for(size_t p=0;p<nn;p++){ rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=0; u32 L=pl[p]; if(L&&pa[L]>=150){
    rgb[3*p]=(u8)(40+(L*97)%216); rgb[3*p+1]=(u8)(40+(L*53)%216); rgb[3*p+2]=(u8)(40+(L*191)%216);} }
  snprintf(fn,sizeof fn,"%s_seg.ppm",base);
  f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
  printf("wrote %s (sheet pieces after delamination + valley-extension cuts)\n",fn); PHASE("valley-extension");

  // (6) RIDGE-graph: extend the spiral fit into the COMPRESSED core. Intensity ridges
  // (sheet centerlines) exist even where there is no air, so use them as graph nodes:
  // a material voxel is a ridge if it is a local intensity max along the RADIAL direction
  // (sheets are ~perpendicular to the radius). Same radial-adjacency graph + winding fit,
  // now covering the air-free blocks too.
  u8*ridge=calloc(nn,1);
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<2)continue; double uy=dy/r,ux=dx/r;
    int iy=(int)lround(y-uy),ix=(int)lround(x-ux), oy=(int)lround(y+uy),ox=(int)lround(x+ux);
    if(iy<0||iy>=d||ix<0||ix>=d||oy<0||oy>=d||ox<0||ox>=d)continue;
    f32 c=vs[p],in=vs[(size_t)iy*d+ix],ou=vs[(size_t)oy*d+ox];
    if(c>=in&&c>=ou&&c>=athr) ridge[p]=1; }
  u32*rl=calloc(nn,sizeof(u32)); u32 nrseg=cc_label(ridge,1,d,d,TOPO_CONN26,rl);
  size_t*ra=calloc((size_t)nrseg+1,sizeof(size_t)); double*rry=calloc((size_t)nrseg+1,sizeof(double));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; u32 L=rl[p]; if(!L)continue;
    ra[L]++; double dy=y-cyf,dx=x-cxf; rry[L]+=sqrt(dy*dy+dx*dx); }
  int*rcid=malloc((size_t)(nrseg+1)*sizeof(int)); int RK=0;
  for(u32 L=1;L<=nrseg;L++) rcid[L]= (ra[L]>=(size_t)minseg)? RK++ : -1;
  if(RK>=2){
    double*rmeanr=calloc(RK,sizeof(double));
    for(u32 L=1;L<=nrseg;L++) if(rcid[L]>=0) rmeanr[rcid[L]]=rry[L]/ra[L];
    double rresid; long rne; double*rwind=graph_winding(ridge,rl,rcid,RK,rmeanr,d,cyf,cxf,walkr,kmin,&rresid,&rne);
    double rwmin=1e30,rwmax=-1e30; for(int i=0;i<RK;i++){ if(rwind[i]<rwmin)rwmin=rwind[i]; if(rwind[i]>rwmax)rwmax=rwind[i]; }
    // coverage: fraction of material within 6px of a ridge node (how much the core is reached)
    size_t mat=0,cov=0; for(size_t p=0;p<nn;p++) if(v[p]>=athr) mat++;
    printf("ridge-graph: nodes=%d edges=%ld  winding %.1f..%.1f (%.0f wraps)  edge-residual=%.3f  ridge-px=%zu\n",
           RK,rne,rwmin,rwmax,rwmax-rwmin,rresid,(size_t)0);
    for(size_t p=0;p<nn;p++){ int g=v[p]/4; rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=(u8)g; }
    // colour by winding with a FULL-RANGE continuous ramp (no hue cycling): map the
    // whole [rwmin,rwmax] span to blue->cyan->green->yellow->red so each winding is a
    // distinct colour and the spiral order is readable inner(blue)->outer(red).
    double rspan = (rwmax>rwmin)? rwmax-rwmin : 1.0;
    for(size_t p=0;p<nn;p++){ u32 L=rl[p]; if(!L||rcid[L]<0)continue; cov++;
      double t=(rwind[rcid[L]]-rwmin)/rspan; if(t<0)t=0; if(t>1)t=1;
      double s=t*4.0; int seg=(int)s; double fr=s-seg; u8 R,G,B;            // 5-stop ramp
      switch(seg){
        case 0: R=0;            G=(u8)(255*fr);   B=255;            break;  // blue->cyan
        case 1: R=0;            G=255;            B=(u8)(255*(1-fr));break; // cyan->green
        case 2: R=(u8)(255*fr); G=255;            B=0;             break;  // green->yellow
        case 3: R=255;          G=(u8)(255*(1-fr));B=0;            break;  // yellow->red
        default:R=255;          G=0;              B=0;             break;  // red
      }
      rgb[3*p]=R;rgb[3*p+1]=G;rgb[3*p+2]=B; }
    printf("  ridge nodes cover %.1f%% of material (vs delaminations only in the core)\n",100.0*cov/mat);
    snprintf(fn,sizeof fn,"%s_ridge.ppm",base);
    f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
    printf("wrote %s (ridge-centerline graph winding, covers compressed core)\n",fn); PHASE("ridge-graph");

    // (7) DELAMINATION-ANCHORED UNIFIED winding: tie the noisy core-reaching ridge
    // graph to the reliable recto (delamination) graph via same-wrap edges, so the
    // recto's spiral consistency propagates into the compressed core.
    double up1,ueq; long un1,uneq;
    double*uw=unified_winding(recto,sl,cid,K,meanr, ridge,rl,rcid,RK, d,cyf,cxf,walkr,tier,kmin,&up1,&ueq,&un1,&uneq);
    if(uw){
      double uwmin=1e30,uwmax=-1e30; for(int i=0;i<K+RK;i++){ if(uw[i]<uwmin)uwmin=uw[i]; if(uw[i]>uwmax)uwmax=uw[i]; }
      printf("UNIFIED: %d recto + %d ridge nodes, %ld +1-edges (resid=%.3f) %ld same-wrap ties (resid=%.3f)  winding %.1f..%.1f (%.0f wraps)\n",
             K,RK,un1,up1,uneq,ueq,uwmin,uwmax,uwmax-uwmin);
      for(size_t p=0;p<nn;p++){ int g=v[p]/4; rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=(u8)g; }
      // colour recto (index cid[L]) and ridge (index Kr+rcid[L]) nodes by unified winding
      for(size_t p=0;p<nn;p++){ int idx=-1; u32 Lr=sl[p],Ld=rl[p];
        if(Lr&&cid[Lr]>=0) idx=cid[Lr]; else if(Ld&&rcid[Ld]>=0) idx=K+rcid[Ld];
        if(idx<0)continue;
        double h=fmod(uw[idx]-uwmin,6.0); if(h<0)h+=6; int hi=(int)h; double fr=h-hi;
        int V=255,P=30,Q=(int)(255*(1-fr)),T=(int)(255*fr); u8 R,G,B;
        switch(hi){case 0:R=V;G=T;B=P;break;case 1:R=Q;G=V;B=P;break;case 2:R=P;G=V;B=T;break;
          case 3:R=P;G=Q;B=V;break;case 4:R=T;G=P;B=V;break;default:R=V;G=P;B=Q;}
        rgb[3*p]=R;rgb[3*p+1]=G;rgb[3*p+2]=B; }
      snprintf(fn,sizeof fn,"%s_unified.ppm",base);
      f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
      printf("wrote %s (delamination-anchored unified winding)\n",fn); PHASE("unified");

      // (8) DENSE PER-WRAP SEGMENTS: propagate the unified node winding to ALL material
      // voxels (Jacobi diffusion within material only; recto+ridge nodes are Dirichlet
      // anchors), then band by integer winding. Each band = ONE CONTINUOUS WRAP -- the
      // half-integer winding contour cuts BETWEEN wraps even where they touch (no air),
      // because the winding interpolates between the two wraps' anchors. The stack of
      // bands is the spiral. This turns thousands of short fragments into ~N long arms.
      // VALLEY BARRIER: a material voxel that is a radial intensity MINIMUM is an
      // inter-sheet boundary (the faint dark line between two touching wraps). Crossing
      // it = a SHEET SWITCH. Mark these so the winding diffusion will not flow across
      // them -- confining each wrap's winding to its own sheet.
      u8 *bar=calloc(nn,1);
      #pragma omp parallel for schedule(static)
      for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
        double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<2)continue; double uy=dy/r,ux=dx/r;
        int iy=(int)lround(y-uy*2),ix=(int)lround(x-ux*2), oy=(int)lround(y+uy*2),ox=(int)lround(x+ux*2);
        if(iy<0||iy>=d||ix<0||ix>=d||oy<0||oy>=d||ox<0||ox>=d)continue;
        float c=vs[p],in=vs[(size_t)iy*d+ix],ou=vs[(size_t)oy*d+ox];
        if(c<=in && c<=ou) bar[p]=1; }                  // radial local-min = inter-sheet boundary
      long nbar=0; for(size_t p=0;p<nn;p++) nbar+=bar[p];
      // ANISOTROPIC weighted diffusion: edge weight is HIGH along bright sheet centers,
      // ~0 across valley barriers -> winding follows the sheet, never jumps to the neighbour.
      float *wd=calloc(nn,sizeof(float)); u8 *anc=calloc(nn,1);
      for(size_t p=0;p<nn;p++){ int idx=-1; u32 Lr=sl[p],Ld=rl[p];
        if(Lr&&cid[Lr]>=0) idx=cid[Lr]; else if(Ld&&rcid[Ld]>=0) idx=K+rcid[Ld];
        if(idx>=0){ wd[p]=(float)uw[idx]; anc[p]=1; } }
      #define EW(q) ( (v[q]<athr || bar[q]) ? 0.0 : (double)(v[q]-athr+1) )   // HARD valley block
      // RED-BLACK GAUSS-SEIDEL with SOR (omega=1.7): in-place (no double buffer), uses the
      // newest neighbour values -> ~3-4x faster convergence than omega=1 Jacobi, so 250
      // sweeps match the old 800. Red/black colouring keeps each colour's cells independent
      // (4-neighbour stencil -> neighbours are the opposite colour) so it stays parallel-safe.
      const double omega=1.7;
      for(int it=0;it<250;it++) for(int color=0;color<2;color++){
        #pragma omp parallel for schedule(static)
        for(int y=0;y<d;y++)for(int x=0;x<d;x++){ if(((x+y)&1)!=color)continue; size_t p=(size_t)y*d+x;
          if(v[p]<athr || bar[p] || anc[p]) continue;     // air + valley walls + anchors fixed
          double s=0,wsum=0,w;
          if(x>0)  {w=EW(p-1);s+=w*wd[p-1];wsum+=w;} if(x<d-1){w=EW(p+1);s+=w*wd[p+1];wsum+=w;}
          if(y>0)  {w=EW(p-d);s+=w*wd[p-d];wsum+=w;} if(y<d-1){w=EW(p+d);s+=w*wd[p+d];wsum+=w;}
          if(wsum>0){ double avg=s/wsum; wd[p]=(float)(wd[p]+omega*(avg-wd[p])); } }
      }
      #undef EW
      printf("valley barriers: %ld inter-sheet-boundary px\n",nbar);
      // band by round(winding); count distinct wraps actually present
      int wlo=(int)floor(uwmin), whi=(int)ceil(uwmax); int nb=whi-wlo+1; if(nb<1)nb=1;
      long *bandpx=calloc(nb,sizeof(long)); long matn=0;
      for(size_t p=0;p<nn;p++) if(v[p]>=athr && !bar[p]){ int w=(int)lround(wd[p])-wlo; if(w<0)w=0; if(w>=nb)w=nb-1; bandpx[w]++; matn++; }
      int nwrap=0; for(int i=0;i<nb;i++) if(bandpx[i]>0) nwrap++;
      printf("DENSE per-wrap bands: %d wraps over %ld material px (avg %ld px/wrap)\n",nwrap,matn,nwrap?matn/nwrap:0);
      // colour each material voxel by its wrap band (full-range ramp, no cycling)
      double sp=(uwmax>uwmin)?uwmax-uwmin:1.0;
      for(size_t p=0;p<nn;p++){ int g=v[p]/5; rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=(u8)g;
        if(v[p]<athr || bar[p]) continue;
        double t=(lround(wd[p])-uwmin)/sp; if(t<0)t=0; if(t>1)t=1;
        double ss=t*4.0; int seg=(int)ss; double fr=ss-seg; u8 R,G,B;
        switch(seg){ case 0:R=0;G=(u8)(255*fr);B=255;break; case 1:R=0;G=255;B=(u8)(255*(1-fr));break;
          case 2:R=(u8)(255*fr);G=255;B=0;break; case 3:R=255;G=(u8)(255*(1-fr));B=0;break; default:R=255;G=0;B=0; }
        rgb[3*p]=R;rgb[3*p+1]=G;rgb[3*p+2]=B; }
      snprintf(fn,sizeof fn,"%s_wraps.ppm",base);
      f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
      printf("wrote %s (dense per-wrap segmentation = the spiral, banded by winding)\n",fn);
      // ALTERNATING discrete colours: consecutive wraps get distinct flat colours so the
      // per-wrap membership is visible as separate concentric rings (a voxel's colour ==
      // which wrap it belongs to). 6-colour cycle on the INTEGER wrap number.
      static const u8 PAL[6][3]={{230,40,40},{40,200,40},{60,90,255},{240,220,40},{220,60,220},{40,220,220}};
      for(size_t p=0;p<nn;p++){ int g=v[p]/6; rgb[3*p]=rgb[3*p+1]=rgb[3*p+2]=(u8)g;
        if(v[p]<athr || bar[p]) continue; int w=(int)lround(wd[p]); int c=((w%6)+6)%6;
        rgb[3*p]=PAL[c][0];rgb[3*p+1]=PAL[c][1];rgb[3*p+2]=PAL[c][2]; }
      snprintf(fn,sizeof fn,"%s_wrapsalt.ppm",base);
      f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nn*3,f); fclose(f); }
      printf("wrote %s (per-wrap, alternating discrete colours)\n",fn); PHASE("dense+propagate");
      free(wd);free(anc);free(bandpx);free(bar);
      free(uw);
    }
  } else printf("ridge-graph: too few ridge segments\n");
  return 0;
}
