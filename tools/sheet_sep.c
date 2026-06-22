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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static double tnow(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+ts.tv_nsec*1e-9; }
// separable box smooth (radius r) of a float field, in place via a scratch buffer.
static void box2d(float *f, float *tmp, int d, int r){
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++){ for(int x=0;x<d;x++){ double s=0;int c=0;
    for(int k=-r;k<=r;k++){ int xx=x+k; if(xx>=0&&xx<d){s+=f[(size_t)y*d+xx];c++;} } tmp[(size_t)y*d+x]=(float)(s/c); } }
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++){ for(int x=0;x<d;x++){ double s=0;int c=0;
    for(int k=-r;k<=r;k++){ int yy=y+k; if(yy>=0&&yy<d){s+=tmp[(size_t)yy*d+x];c++;} } f[(size_t)y*d+x]=(float)(s/c); } }
}

// ACROSS-WRAP DIRECTION field: the unit direction pointing from one wrap to the NEXT,
// used by ridge/valley/recto detection and the graph walk INSTEAD of the global radius.
// NULL g_dirx => fall back to radial. Two non-radial modes let the tracing follow the
// real DEFORMED sheets (flat tops, bulges, concavities) the radius rounds toward a circle.
static float *g_dirx=NULL,*g_diry=NULL;
#define WDIRV(P,Y,X) double dy=(double)(Y)-cyf,dx=(double)(X)-cxf,r=sqrt(dy*dy+dx*dx),uy,ux; \
  if(g_dirx&&(g_dirx[P]!=0.0f||g_diry[P]!=0.0f)){ ux=g_dirx[P]; uy=g_diry[P]; } \
  else { double rr=(r<1e-6)?1e-6:r; uy=dy/rr; ux=dx/rr; }

// MODE "normal": local sheet normal = dominant eigenvector of the smoothed gradient
// structure tensor, oriented OUTWARD (away from the umbilicus). Follows local deformation.
static void dir_normal(const float*vs,const u8*v,int athr,int d,double cyf,double cxf,float*dx,float*dy){
  size_t nn=(size_t)d*d;
  float *Jxx=calloc(nn,sizeof(float)),*Jyy=calloc(nn,sizeof(float)),*Jxy=calloc(nn,sizeof(float)),*tmp=malloc(nn*sizeof(float));
  #pragma omp parallel for schedule(static)
  for(int y=1;y<d-1;y++)for(int x=1;x<d-1;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
    float gx=0.5f*(vs[p+1]-vs[p-1]),gy=0.5f*(vs[p+d]-vs[p-d]); Jxx[p]=gx*gx;Jyy[p]=gy*gy;Jxy[p]=gx*gy; }
  box2d(Jxx,tmp,d,3); box2d(Jyy,tmp,d,3); box2d(Jxy,tmp,d,3);
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr){dx[p]=dy[p]=0;continue;}
    double a=Jxx[p],b=Jxy[p],c=Jyy[p],tr=a+c,det=a*c-b*b,disc=sqrt(fmax(0,tr*tr/4-det)),l1=tr/2+disc;
    double ex=b,ey=l1-a,nm=hypot(ex,ey); if(nm<1e-9){ex=1;ey=0;nm=1;} ex/=nm;ey/=nm;
    double rx=x-cxf,ry=y-cyf; if(ex*rx+ey*ry<0){ex=-ex;ey=-ey;}   // orient outward
    dx[p]=(float)ex; dy[p]=(float)ey; }
  free(Jxx);free(Jyy);free(Jxy);free(tmp);
}
// MODE "envelope": warp the radius by the per-angle outer envelope R(theta) (the egg
// shape), then the across-wrap direction = normalized gradient of the warped radius.
static void dir_envelope(const u8*v,int athr,int d,double cyf,double cxf,float*dx,float*dy){
  size_t nn=(size_t)d*d; int NB=360;
  double *R=calloc(NB,sizeof(double));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ if(v[(size_t)y*d+x]<athr)continue;
    double ry=y-cyf,rx=x-cxf,r=sqrt(rx*rx+ry*ry),th=atan2(ry,rx); int b=(int)((th+M_PI)/(2*M_PI)*NB); b=((b%NB)+NB)%NB;
    if(r>R[b])R[b]=r; }
  double *Rs=calloc(NB,sizeof(double)); for(int b=0;b<NB;b++){ double s=0;int c=0;
    for(int k=-8;k<=8;k++){int bb=((b+k)%NB+NB)%NB; if(R[bb]>0){s+=R[bb];c++;}} Rs[b]=c?s/c:R[b]; }
  double Rref=0;int cc=0; for(int b=0;b<NB;b++) if(Rs[b]>0){Rref+=Rs[b];cc++;} Rref=cc?Rref/cc:1;
  float *rw=malloc(nn*sizeof(float));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; double ry=y-cyf,rx=x-cxf,r=sqrt(rx*rx+ry*ry),th=atan2(ry,rx);
    int b=(int)((th+M_PI)/(2*M_PI)*NB); b=((b%NB)+NB)%NB; double Rb=Rs[b]>1?Rs[b]:1; rw[p]=(float)(r*Rref/Rb); }
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr){dx[p]=dy[p]=0;continue;}
    float gx=(x>0&&x<d-1)?0.5f*(rw[p+1]-rw[p-1]):0, gy=(y>0&&y<d-1)?0.5f*(rw[p+d]-rw[p-d]):0;
    double nm=hypot(gx,gy);
    if(nm<1e-6){ double rx=x-cxf,ry=y-cyf,r=hypot(rx,ry); if(r<1){dx[p]=dy[p]=0;} else {dx[p]=(float)(rx/r);dy[p]=(float)(ry/r);} }
    else { dx[p]=(float)(gx/nm); dy[p]=(float)(gy/nm); } }
  free(R);free(Rs);free(rw);
}
// MODE "hybrid": envelope (robust global egg) refined by the local normal where the
// structure tensor is COHERENT (clear sheet). dir = normalize((1-w)*env + w*normal),
// w = 0.6*coherence -- noisy regions stay on the envelope, clear regions get the local
// wiggle. The blended field is box-smoothed + renormalized to suppress residual jitter.
static void dir_hybrid(const float*vs,const u8*v,int athr,int d,double cyf,double cxf,float*dx,float*dy){
  size_t nn=(size_t)d*d;
  float *ex=calloc(nn,sizeof(float)),*ey=calloc(nn,sizeof(float));
  dir_envelope(v,athr,d,cyf,cxf,ex,ey);                       // global egg base
  float *Jxx=calloc(nn,sizeof(float)),*Jyy=calloc(nn,sizeof(float)),*Jxy=calloc(nn,sizeof(float)),*tmp=malloc(nn*sizeof(float));
  #pragma omp parallel for schedule(static)
  for(int y=1;y<d-1;y++)for(int x=1;x<d-1;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
    float gx=0.5f*(vs[p+1]-vs[p-1]),gy=0.5f*(vs[p+d]-vs[p-d]); Jxx[p]=gx*gx;Jyy[p]=gy*gy;Jxy[p]=gx*gy; }
  box2d(Jxx,tmp,d,3); box2d(Jyy,tmp,d,3); box2d(Jxy,tmp,d,3);
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr){dx[p]=dy[p]=0;continue;}
    double a=Jxx[p],b=Jxy[p],c=Jyy[p],tr=a+c,det=a*c-b*b,disc=sqrt(fmax(0,tr*tr/4-det)),l1=tr/2+disc,l0=tr/2-disc;
    double nx=b,ny=l1-a,nm=hypot(nx,ny); if(nm<1e-9){nx=ex[p];ny=ey[p];nm=hypot(nx,ny);if(nm<1e-9){nm=1;}} nx/=nm;ny/=nm;
    if(nx*ex[p]+ny*ey[p]<0){nx=-nx;ny=-ny;}                    // orient normal to the envelope
    double coh=(l1+l0>1e-9)?(l1-l0)/(l1+l0):0.0, w=0.6*coh;
    double bx=(1-w)*ex[p]+w*nx, by=(1-w)*ey[p]+w*ny, bn=hypot(bx,by); if(bn<1e-9){bx=ex[p];by=ey[p];bn=hypot(bx,by);if(bn<1e-9)bn=1;}
    dx[p]=(float)(bx/bn); dy[p]=(float)(by/bn); }
  box2d(dx,tmp,d,2); box2d(dy,tmp,d,2);                       // smooth the blended direction
  #pragma omp parallel for schedule(static)
  for(size_t p=0;p<nn;p++){ if(v[p]<athr){dx[p]=dy[p]=0;continue;} double n=hypot(dx[p],dy[p]); if(n>1e-9){dx[p]/=n;dy[p]/=n;} }
  free(ex);free(ey);free(Jxx);free(Jyy);free(Jxy);free(tmp);
}
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

/* VALLEY-CROSSING COUNT: walk the smoothed profile vs() outward from (y,x) along the
 * unit direction (uy,ux) for t in [1,kstop] and count the inter-sheet intensity minima
 * (valleys) crossed -- this is the TRUE winding step between two radial-adjacent markers,
 * replacing the blanket "+1" assumption:
 *   1 valley  -> immediate next wrap        (+1, the common case)
 *   2+ valleys-> a wrap was MISSED between   (+2.., fixes skip under-counting)
 *   0 valleys -> touching sheets / same-wrap fragment, no separating gap (tie, NOT +1;
 *                this is the sheet-switch fix: don't assert a jump without evidence)
 * A zigzag/hysteresis detector with prominence `prom` is robust to noise; air gaps
 * (vs==0) are deep minima always counted. */
static int count_valleys(const float*vs,int d,double y,double x,double uy,double ux,
                         int kstop,double prom,double pitchmin){
  if(prom<0) return 1;            // sentinel: baseline mode (assume +1, no valley counting)
  if(kstop<2) return 1;
  // Extract the alternating sequence of PROMINENT extrema (peaks=ridges, valleys=gaps)
  // along the ray. A and B are themselves ridges -> bracket with implicit peaks at the
  // endpoints. Then a valley counts as a wrap boundary ONLY if its two FLANKING ridges
  // are >= pitchmin apart: a true inter-wrap gap spans ~1 pitch, whereas an intra-wrap
  // delamination slit has its flanking ridges much closer -> rejected (the fix for the
  // false +2 over-count from cracks inside a single wrap).
  // Ridges (peaks) are confirmed with a LARGER prominence than valleys: a real sheet
  // centre stands well above the gaps, so requiring a bigger rise to call something a
  // ridge stops intra-gap noise bumps from splitting one inter-wrap gap into two sub-
  // valleys (which the pitchmin gate would then both reject -> under-count).
  double peakprom = 1.8*prom;
  double pos[260]; signed char typ[260]; int ne=0;             // typ: +1 ridge, -1 valley
  pos[ne]=0; typ[ne]=1; ne++;                                  // A is a ridge
  double cmax=-1e30,cmin=1e30,maxp=0,minp=0; int dir=0;
  for(int t=0;t<=kstop;t++){
    int iy=(int)lround(y+uy*t),ix=(int)lround(x+ux*t);
    double val = (iy<0||iy>=d||ix<0||ix>=d)?0.0:(double)vs[(size_t)iy*d+ix];
    if(t==0){ cmax=cmin=val; maxp=minp=0; continue; }
    if(dir>=0){ if(val>cmax){cmax=val;maxp=t;}
      if(cmax-val>=peakprom){ if(ne<258&&typ[ne-1]!=1){pos[ne]=maxp;typ[ne]=1;ne++;} dir=-1; cmin=val; minp=t; } }
    if(dir<=0){ if(val<cmin){cmin=val;minp=t;}
      if(val-cmin>=prom){ if(ne<258&&typ[ne-1]!=-1){pos[ne]=minp;typ[ne]=-1;ne++;} dir=1; cmax=val; maxp=t; } }
  }
  if(typ[ne-1]!=1){ pos[ne]=kstop; typ[ne]=1; ne++; }          // B is a ridge
  int nv=0;
  for(int i=1;i<ne-1;i++) if(typ[i]==-1){                      // valley flanked by peaks i-1, i+1
    if(pos[i+1]-pos[i-1] >= pitchmin) nv++; }
  return nv;
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
                             const float *vs, double vprom, double pitchmin, double *residual, long *nedges){
  long *edge=calloc((size_t)K*K,sizeof(long));     // vote count per (a,b)
  long *step=calloc((size_t)K*K,sizeof(long));      // sum of valley-counts per (a,b)
  // Walk outward to the IMMEDIATE different-segment neighbour, then COUNT the inter-sheet
  // valleys crossed to get the true winding step (0=touch/fragment, 1=adjacent, 2+=skip)
  // instead of assuming +1. kmin still rejects sub-half-pitch adjacency noise.
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!marker[p])continue;
    int La=cid[lbl[p]]; if(La<0)continue;
    WDIRV(p,y,x); if(r<1)continue;
    for(int k=2;k<wmax;k++){ int ny_=(int)lround(y+uy*k),nx_=(int)lround(x+ux*k);
      if(ny_<0||ny_>=d||nx_<0||nx_>=d)break; size_t q=(size_t)ny_*d+nx_;
      if(marker[q]){ int Lb=cid[lbl[q]]; if(Lb>=0&&Lb!=La&&k>=kmin){
        int nv=count_valleys(vs,d,y,x,uy,ux,k,vprom,pitchmin);
        edge[(size_t)La*K+Lb]++; step[(size_t)La*K+Lb]+=nv; break; } } } }
  double *wind=calloc(K,sizeof(double));
  int inner=0; for(int i=1;i<K;i++) if(meanr[i]<meanr[inner]) inner=i;
  // SPARSE edge list: the K*K matrix is ~99.9% zero; relaxing over the nonzero edges
  // is O(iters*E) not O(iters*K*K) (an ~800x speedup at K~1000). es[] = majority step.
  long ne=0; for(int a=0;a<K;a++)for(int b=0;b<K;b++) if(edge[(size_t)a*K+b]>=3) ne++;
  int *ea=malloc((size_t)(ne?ne:1)*sizeof(int)),*eb=malloc((size_t)(ne?ne:1)*sizeof(int));
  long *ew=malloc((size_t)(ne?ne:1)*sizeof(long)); double *es=malloc((size_t)(ne?ne:1)*sizeof(double)); long ei=0;
  for(int a=0;a<K;a++)for(int b=0;b<K;b++){ long w=edge[(size_t)a*K+b];
    if(w>=3){ ea[ei]=a;eb[ei]=b;ew[ei]=w; es[ei]=(double)lround((double)step[(size_t)a*K+b]/w); ei++; } }
  free(edge); free(step);
  double *nw=malloc(K*sizeof(double)),*den=malloc(K*sizeof(double));
  double *rw=malloc((size_t)(ne?ne:1)*sizeof(double)); for(long e=0;e<ne;e++) rw[e]=1.0; // IRLS robust weights
  // IRLS: solve wind[b]=wind[a]+es[e], then DOWN-WEIGHT edges whose residual is far from
  // their integer step (mis-counts / stray bridges), and re-solve. w = 1/(1+(res/0.4)^2).
  for(int round=0;round<4;round++){
    for(int it=0;it<600;it++){ memset(nw,0,K*sizeof(double)); memset(den,0,K*sizeof(double));
      for(long e=0;e<ne;e++){ int a=ea[e],b=eb[e]; double w=ew[e]*rw[e], s=es[e];
        nw[a]+=w*(wind[b]-s); den[a]+=w; nw[b]+=w*(wind[a]+s); den[b]+=w; }
      den[inner]+=5;
      for(int i=0;i<K;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
    if(round<3) for(long e=0;e<ne;e++){ double res=wind[eb[e]]-wind[ea[e]]-es[e]; double z=res/0.4; rw[e]=1.0/(1.0+z*z); }
  }
  // report INLIER residual (robust weight > 0.5 == within ~0.4 of its integer step)
  double sres=0; long nin=0; for(long e=0;e<ne;e++){ if(rw[e]>0.5){ sres+=fabs(wind[eb[e]]-wind[ea[e]]-es[e]); nin++; } }
  free(nw);free(den);free(rw);free(ea);free(eb);free(ew);free(es);
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
// Aggregate a raw list of edge keys into sparse arrays, keeping keys with run>=minc.
// When packbits>0 each raw entry is (key<<packbits | step); entries are grouped by key
// and the per-edge step (es, rounded mean) is returned in *so. packbits==0 => plain keys.
static long agg_edges(long*raw,long nraw,int K,int**ao,int**bo,long**wo,long**so,long minc,int packbits){
  long M = packbits>0 ? (1L<<packbits) : 1;
  if(nraw==0){ *ao=malloc(sizeof(int)); *bo=malloc(sizeof(int)); *wo=malloc(sizeof(long)); if(so)*so=malloc(sizeof(long)); return 0; }
  qsort(raw,nraw,sizeof(long),cmp_long);
  long ne=0; for(long i=0;i<nraw;){ long key=raw[i]/M; long j=i; while(j<nraw&&raw[j]/M==key)j++; if(j-i>=minc)ne++; i=j; }
  int*ea=malloc((size_t)(ne?ne:1)*sizeof(int)),*eb=malloc((size_t)(ne?ne:1)*sizeof(int));
  long*ew=malloc((size_t)(ne?ne:1)*sizeof(long)); long*es= so?malloc((size_t)(ne?ne:1)*sizeof(long)):NULL;
  long e=0; for(long i=0;i<nraw;){ long key=raw[i]/M; long j=i; long ss=0; while(j<nraw&&raw[j]/M==key){ ss+=raw[j]%M; j++; } long c=j-i;
    if(c>=minc){ ea[e]=(int)(key/K); eb[e]=(int)(key%K); ew[e]=c; if(es) es[e]=lround((double)ss/c); e++; } i=j; }
  *ao=ea; *bo=eb; *wo=ew; if(so)*so=es; return ne;
}
static double *unified_winding(const u8*recto,const u32*rsl,const int*rcid,int Kr,const double*rmeanr,
                               const u8*ridge,const u32*rdl,const int*rdcid,int Kd,
                               int d,double cyf,double cxf,int wmax,int tiemax,int kmin,
                               const float*vs,double vprom,double pitchmin,
                               double*resid_p1,double*resid_eq,long*np1,long*neq){
  int K=Kr+Kd; if(K<2) return NULL;
  // SPARSE edge accumulation: push (key<<4 | valley-step) into growable lists, then
  // sort+aggregate -- O(E log E) over the ~10^5 real edges, NOT O(K*K) (214M at L2).
  // The valley-step (count_valleys) replaces the assumed +1: 0=touch/fragment tie,
  // 1=adjacent, 2+=skip across a missed wrap.
  long cap1=4096,nr1=0; long*raw1=malloc(cap1*sizeof(long));
  long cap0=4096,nr0=0; long*raw0=malloc(cap0*sizeof(long));
  #define PUSH1(key) do{ if(nr1>=cap1){cap1*=2;raw1=realloc(raw1,cap1*sizeof(long));} raw1[nr1++]=(key); }while(0)
  #define PUSH0(key) do{ if(nr0>=cap0){cap0*=2;raw0=realloc(raw0,cap0*sizeof(long));} raw0[nr0++]=(key); }while(0)
  #define RADWALK(MARK,LBL,CID,BASE) \
    for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!(MARK)[p])continue; \
      int La=(CID)[(LBL)[p]]; if(La<0)continue; WDIRV(p,y,x); if(r<1)continue; \
      for(int k=2;k<wmax;k++){ int ny_=(int)lround(y+uy*k),nx_=(int)lround(x+ux*k); \
        if(ny_<0||ny_>=d||nx_<0||nx_>=d)break; size_t q=(size_t)ny_*d+nx_; \
        if((MARK)[q]){ int Lb=(CID)[(LBL)[q]]; if(Lb>=0&&Lb!=La&&k>=kmin){ \
          int nv=count_valleys(vs,d,y,x,uy,ux,k,vprom,pitchmin); if(nv>15)nv=15; \
          PUSH1((((long)((BASE)+La)*K+(BASE)+Lb)<<4)|nv); break; } } } }
  RADWALK(recto,rsl,rcid,0);
  RADWALK(ridge,rdl,rdcid,Kr);
  #undef RADWALK
  // same-wrap ties: from each ridge voxel walk OUTWARD; first recto hit = same wrap.
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!ridge[p])continue;
    int Ld=rdcid[rdl[p]]; if(Ld<0)continue; WDIRV(p,y,x); if(r<1)continue;
    for(int k=1;k<=tiemax;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k);
      if(oy<0||oy>=d||ox<0||ox>=d)break; size_t q=(size_t)oy*d+ox;
      if(recto[q]){ int Lr=rcid[rsl[q]]; if(Lr>=0){ PUSH0((long)(Kr+Ld)*K+Lr); } break; } } }
  #undef PUSH1
  #undef PUSH0
  int inner=0; for(int i=1;i<Kr;i++) if(rmeanr[i]<rmeanr[inner]) inner=i;
  const double WEQ=3.0;
  double *wind=calloc(K,sizeof(double));
  int *a1,*b1,*a0,*b0; long *w1l,*w0l,*s1l;
  long n1=agg_edges(raw1,nr1,K,&a1,&b1,&w1l,&s1l,3,4);   // packed step in low 4 bits
  long n0=agg_edges(raw0,nr0,K,&a0,&b0,&w0l,NULL,2,0);   // ties: step always 0
  free(raw1);free(raw0);
  double*nw=malloc(K*sizeof(double)),*den=malloc(K*sizeof(double));
  for(int it=0;it<3000;it++){ memset(nw,0,K*sizeof(double)); memset(den,0,K*sizeof(double));
    for(long e=0;e<n1;e++){ int a=a1[e],b=b1[e]; long w=w1l[e]; double s=(double)s1l[e]; nw[a]+=w*(wind[b]-s); den[a]+=w; nw[b]+=w*(wind[a]+s); den[b]+=w; }
    for(long e=0;e<n0;e++){ int a=a0[e],b=b0[e]; double w=WEQ*w0l[e]; nw[a]+=w*wind[b]; den[a]+=w; nw[b]+=w*wind[a]; den[b]+=w; }
    den[inner]+=5;                               // pull anchor toward 0
    for(int i=0;i<K;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]); }
  double s1=0,s0=0;
  for(long e=0;e<n1;e++) s1+=fabs(wind[b1[e]]-wind[a1[e]]-(double)s1l[e]);
  for(long e=0;e<n0;e++) s0+=fabs(wind[b0[e]]-wind[a0[e]]);
  free(nw);free(den);free(a1);free(b1);free(w1l);free(s1l);free(a0);free(b0);free(w0l);
  if(resid_p1)*resid_p1=n1?s1/n1:0; if(resid_eq)*resid_eq=n0?s0/n0:0; if(np1)*np1=n1; if(neq)*neq=n0;
  return wind;
}

int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: %s ARCHIVE OUTBASE lod z0 y0 x0 d [minseg=15] [dirmode=envelope] [valley-prom=0.35*athr]\n",argv[0]);return 2;}
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]); int d=atoi(argv[7]);
  int minseg=argc>8?atoi(argv[8]):15;
  const char*dirmode=argc>9?argv[9]:"auto"; // auto | radial | normal | envelope | hybrid (across-wrap direction)
                                                 // envelope wins on real data: global egg shape is robust,
                                                 // local structure-tensor normals too noisy at these LODs.
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

  // smoothed intensity (air->0), used by ridge/valley detection + the normal direction.
  f32 *vs=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vs[p]=v[p];
  { f32*t=malloc(nn*sizeof(f32)); for(int it=0;it<2;it++){
      for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; double s=vs[p];int c=1;
        if(x>0){s+=vs[p-1];c++;}if(x<d-1){s+=vs[p+1];c++;}if(y>0){s+=vs[p-d];c++;}if(y<d-1){s+=vs[p+d];c++;}
        t[p]=(f32)(s/c);} memcpy(vs,t,nn*sizeof(f32)); } free(t); }
  for(size_t p=0;p<nn;p++) if(v[p]<athr) vs[p]=0;

  // PITCH (vox/wrap) -- measured up-front; also the resolution proxy for AUTO dir mode.
  double pitch = measure_pitch(v,d,cyf,cxf,athr);

  // ACROSS-WRAP DIRECTION (radial | normal | envelope | hybrid | auto) -- lets detection
  // and the graph walk follow the real deformed sheets instead of the global radius.
  // AUTO: structure-tensor normals are reliable only when sheets are resolved over
  // several voxels, so pick by pitch -- fine data (pitch>=16, ~L2 and below) uses HYBRID
  // (local-normal deformation wins on real data); coarse data uses ENVELOPE (normals too
  // noisy when downscaled, global egg-shape wins). Verified by metrics: L2 hybrid geom
  // 0.05 / backward 15.7 vs envelope 0.14 / 19.0; L3 envelope wins.
  char dmode[16]; snprintf(dmode,sizeof dmode,"%s", strcmp(dirmode,"auto")? dirmode : (pitch>=16.0?"hybrid":"envelope"));
  if(strcmp(dmode,"normal")==0){ g_dirx=calloc(nn,sizeof(float)); g_diry=calloc(nn,sizeof(float)); dir_normal(vs,v,athr,d,cyf,cxf,g_dirx,g_diry); }
  else if(strcmp(dmode,"envelope")==0){ g_dirx=calloc(nn,sizeof(float)); g_diry=calloc(nn,sizeof(float)); dir_envelope(v,athr,d,cyf,cxf,g_dirx,g_diry); }
  else if(strcmp(dmode,"hybrid")==0){ g_dirx=calloc(nn,sizeof(float)); g_diry=calloc(nn,sizeof(float)); dir_hybrid(vs,v,athr,d,cyf,cxf,g_dirx,g_diry); }
  fprintf(stderr,"across-wrap direction mode: %s%s\n",dmode,strcmp(dirmode,"auto")?"":" (auto)");

  // (1) recto faces: material with air OUTWARD along the radius
  u8*recto=calloc(nn,1);
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
    WDIRV(p,y,x); if(r<1)continue;
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
  int walkr=(int)(1.6*pitch+0.5); if(walkr<8)walkr=8; if(walkr>120)walkr=120;
  int kmin =(int)(0.5*pitch+0.5); if(kmin<2)kmin=2;            // reject same-wrap-fragment +1 edges
  int tier =(int)(0.7*pitch+0.5); if(tier<3)tier=3; if(tier>30)tier=30;
  // valley-prominence for the winding-step counter: a dip must drop this many u8 below a
  // flanking ridge to count as an inter-sheet boundary. Scaled to the air cut (faint
  // intra-core minima sit a fraction of athr below the sheet centres); air gaps (vs==0)
  // always clear it. Tunable via argv (default 0.35*athr).
  double vprom = argc>10? atof(argv[10]) : 0.35*athr; if(vprom>=0 && vprom<6) vprom=6;  // <0 = baseline +1 mode
  double pitchmin = (argc>11? atof(argv[11]) : 0.20)*pitch;   // a counted valley's flanking ridges must span >= this*pitch (reject intra-wrap slits);
                                                              // 0.20 = geom-optimal on L3 (wraps match span/pitch); higher trims switches at the cost of under-count
  fprintf(stderr,"radial pitch=%.1f vox/wrap -> walk=%d tie=%d  valley-prom=%.1f\n",pitch,walkr,tier,vprom); PHASE("threshold+pitch");

  // (3+4) radial-adjacency graph over the recto rims + winding fit (the helper)
  double resid; long ne; double *wind=graph_winding(recto,sl,cid,K,meanr,d,cyf,cxf,walkr,kmin,vs,vprom,pitchmin,&resid,&ne);
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
    WDIRV(p,y,x); if(r<2)continue;
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
    double rresid; long rne; double*rwind=graph_winding(ridge,rl,rcid,RK,rmeanr,d,cyf,cxf,walkr,kmin,vs,vprom,pitchmin,&rresid,&rne);
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
    double*uw=unified_winding(recto,sl,cid,K,meanr, ridge,rl,rcid,RK, d,cyf,cxf,walkr,tier,kmin,vs,vprom,pitchmin,&up1,&ueq,&un1,&uneq);
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
        WDIRV(p,y,x); if(r<2)continue;
        int iy=(int)lround(y-uy*2),ix=(int)lround(x-ux*2), oy=(int)lround(y+uy*2),ox=(int)lround(x+ux*2);
        if(iy<0||iy>=d||ix<0||ix>=d||oy<0||oy>=d||ox<0||ox>=d)continue;
        float c=vs[p],in=vs[(size_t)iy*d+ix],ou=vs[(size_t)oy*d+ox];
        if(c<=in && c<=ou) bar[p]=1; }                  // radial local-min = inter-sheet boundary
      // BARRIER CLOSE: seal azimuthal GAPS in the inter-sheet walls (winding leaks through a
      // gap and contaminates an arc -> the coherent tangential backward-switch streaks seen
      // near barriers). Morphological close (dilate r then erode r) over MATERIAL only, so
      // it bridges gaps without thickening into the sheets. argv[14] radius (0 = off).
      int barclose=argc>14?atoi(argv[14]):(pitch>=16.0?1:0);  // fine res: seal salty barrier gaps (cuts switches);
                                                              // coarse res: barriers already thick, closing over-merges -> off
      if(barclose>0){ u8*tb=malloc(nn);
        for(int it2=0;it2<barclose;it2++){ memcpy(tb,bar,nn);
          for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t q=(size_t)y*d+x; if(tb[q]||v[q]<athr)continue;
            if((x>0&&tb[q-1])||(x<d-1&&tb[q+1])||(y>0&&tb[q-d])||(y<d-1&&tb[q+d])) bar[q]=1; } }
        for(int it2=0;it2<barclose;it2++){ memcpy(tb,bar,nn);
          for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t q=(size_t)y*d+x; if(!tb[q])continue;
            if((x>0&&!tb[q-1])||(x<d-1&&!tb[q+1])||(y>0&&!tb[q-d])||(y<d-1&&!tb[q+d])) bar[q]=0; } }
        free(tb); }
      long nbar=0; for(size_t p=0;p<nn;p++) nbar+=bar[p];
      // ANISOTROPIC weighted diffusion: edge weight is HIGH along bright sheet centers,
      // ~0 across valley barriers -> winding follows the sheet, never jumps to the neighbour.
      float *wd=calloc(nn,sizeof(float)); u8 *anc=calloc(nn,1); float *ancw=calloc(nn,sizeof(float));
      for(size_t p=0;p<nn;p++){ int idx=-1; u32 Lr=sl[p],Ld=rl[p];
        if(Lr&&cid[Lr]>=0) idx=cid[Lr]; else if(Ld&&rcid[Ld]>=0) idx=K+rcid[Ld];
        if(idx>=0){ wd[p]=(float)uw[idx]; ancw[p]=(float)uw[idx]; anc[p]=1; } }
      // SHEET-NORMAL field (2D structure tensor): make the propagation follow the LOCAL
      // deformed sheet instead of the global radius. n = dominant eigenvector of the
      // smoothed gradient tensor (points ACROSS the sheet). Diffusing ALONG the tangent
      // (perp to n) and NOT across n preserves the real wrap shape (flat tops, bulges,
      // concavities) that isotropic diffusion rounds toward a circle.
      float *nrx=calloc(nn,sizeof(float)),*nry=calloc(nn,sizeof(float));
      { float *Jxx=calloc(nn,sizeof(float)),*Jyy=calloc(nn,sizeof(float)),*Jxy=calloc(nn,sizeof(float)),*tmp=malloc(nn*sizeof(float));
        #pragma omp parallel for schedule(static)
        for(int y=1;y<d-1;y++)for(int x=1;x<d-1;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr)continue;
          float gx=0.5f*(vs[p+1]-vs[p-1]), gy=0.5f*(vs[p+d]-vs[p-d]); Jxx[p]=gx*gx;Jyy[p]=gy*gy;Jxy[p]=gx*gy; }
        box2d(Jxx,tmp,d,3); box2d(Jyy,tmp,d,3); box2d(Jxy,tmp,d,3);   // integration scale ~3 vox
        #pragma omp parallel for schedule(static)
        for(size_t p=0;p<nn;p++){ if(v[p]<athr)continue;
          double a=Jxx[p],b=Jxy[p],c=Jyy[p], tr=a+c, det=a*c-b*b, disc=sqrt(fmax(0,tr*tr/4-det)), l1=tr/2+disc;
          double ex=b, ey=l1-a, nn2=hypot(ex,ey); if(nn2<1e-9){ex=1;ey=0;nn2=1;} nrx[p]=(float)(ex/nn2); nry[p]=(float)(ey/nn2); }
        free(Jxx);free(Jyy);free(Jxy);free(tmp); }
      #define EW(q) ( (v[q]<athr || bar[q]) ? 0.0 : (double)(v[q]-athr+1) )   // HARD valley block
      // RED-BLACK GAUSS-SEIDEL with SOR (omega=1.7), ANISOTROPIC: a neighbour in direction
      // e is weighted by (1-(e.n)^2) -- ~1 along the sheet tangent, ~0 across the normal.
      // SOFT ANCHORS: the graph nodes are NOT hard-fixed (that pinned the field to their
      // per-node winding NOISE -> round(wd) flickered azimuthally -> wraps shattered into
      // ~30 fragments). Instead each node adds a DATA TERM pulling wd[p] toward its winding
      // with weight LAM, while still participating in the tangential smoothing. The field
      // tracks the nodes globally but averages out their azimuthal jitter -> continuous wraps.
      const double omega=1.7, LAM=(argc>12?atof(argv[12]):0.1);  // soft-anchor strength; 0.1 = switch-rate parity with hard +1 baseline at full geom
      // OUTWARD-DRIFT (screened Poisson): winding MUST climb outward at ~1/pitch across the
      // wraps (it's a spiral). A pure harmonic (Laplace) solve has no slope -> wherever
      // anchors locally disagree the field sits flat and wiggles, producing backward sheet
      // switches. Each neighbour contributes its winding MINUS the expected step along the
      // across-wrap direction, so the field is pushed to climb at the right rate. DR scales
      // the drift (argv[13]; 0 = old harmonic).
      const double DR=(argc>13?atof(argv[13]):0.0)/fmax(pitch,1e-6);  // outward drift: 0 (field already has correct slope; forcing it hurt)
      for(int it=0;it<400;it++) for(int color=0;color<2;color++){
        #pragma omp parallel for schedule(static)
        for(int y=0;y<d;y++)for(int x=0;x<d;x++){ if(((x+y)&1)!=color)continue; size_t p=(size_t)y*d+x;
          if(v[p]<athr || bar[p]) continue;               // air + valley walls fixed; anchors now SOFT
          double nx_=nrx[p],ny_=nry[p];
          double ax=0.12+0.88*(1.0-nx_*nx_), ay=0.12+0.88*(1.0-ny_*ny_);  // tangent-biased, small isotropic floor
          WDIRV(p,y,x); double gx=DR*ux, gy=DR*uy;          // expected d(wind) per +1 step in x / y (outward)
          double s=0,wsum=0,w;
          if(x>0)  {w=EW(p-1)*ax;s+=w*(wd[p-1]+gx);wsum+=w;} if(x<d-1){w=EW(p+1)*ax;s+=w*(wd[p+1]-gx);wsum+=w;}
          if(y>0)  {w=EW(p-d)*ay;s+=w*(wd[p-d]+gy);wsum+=w;} if(y<d-1){w=EW(p+d)*ay;s+=w*(wd[p+d]-gy);wsum+=w;}
          if(anc[p]){ double wa=LAM*wsum+1e-6; s+=wa*ancw[p]; wsum+=wa; }   // soft data term toward node winding
          if(wsum>1e-9){ double avg=s/wsum; wd[p]=(float)(wd[p]+omega*(avg-wd[p])); } }
      }
      #undef EW
      free(nrx);free(nry);free(ancw);
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

      // ===== QUALITY METRICS: grade the spiral fit from the dense winding field wd[] =====
      // These measure the things the edge-residual does NOT: skip/under-count (geometry),
      // sheet switches (radial monotonicity), and per-wrap integrity (coverage/fragments/holes).
      {
        // (a) GEOMETRIC CONSISTENCY: a correct fit has wraps ~= radial-span / pitch.
        long *rh=calloc(4096,sizeof(long)); long rn=0;
        for(size_t p=0;p<nn;p++) if(v[p]>=athr){ int yy=(int)(p/d),xx=(int)(p%d); double rr=hypot(yy-cyf,xx-cxf); int ri=(int)rr; if(ri>=0&&ri<4096){rh[ri]++;rn++;} }
        double Rin=0,Rout=0; long acc=0; int gotin=0;
        for(int i=0;i<4096;i++){ acc+=rh[i]; if(!gotin && acc>=rn*0.02){Rin=i;gotin=1;} if(acc>=rn*0.98){Rout=i;break;} }
        free(rh);
        double span=Rout-Rin, geomw=span/fmax(pitch,1e-6);
        printf("METRIC geom-consistency  wraps=%d  pitch=%.1f  radial-span=%.0f  implied-wraps=%.1f  rel.err=%.2f\n",
               nwrap,pitch,span,geomw,fabs(nwrap-geomw)/fmax(geomw,1));

        // (b) SWITCH RATE -- VOXEL-FAIR (each material voxel counted once; radial rays
        // over-sample the core ~360x and inflate it). For each voxel, its outward neighbour
        // along the across-wrap direction should have HIGHER winding; a DROP is a backward
        // sheet switch. Split inner-third vs outer to localise; flag at-barrier (touch).
        long sw=0,jumps=0,matc=0, sw_in=0,mat_in=0, sw_bar=0; int Rr=(int)Rout+4; (void)Rr;
        double rcore=Rin+(Rout-Rin)*0.34;
        u8 *swmap=calloc(nn,1);                  // heatmap: 1=backward switch, 2=wrap jump
        for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(v[p]<athr||bar[p])continue;
          WDIRV(p,y,x); if(r<2)continue;
          int oy=(int)lround(y+uy*2),ox=(int)lround(x+ux*2); if(oy<0||oy>=d||ox<0||ox>=d)continue;
          size_t o=(size_t)oy*d+ox; if(v[o]<athr||bar[o])continue;
          matc++; int inr=(r<rcore); if(inr)mat_in++;
          double diff=(double)wd[o]-(double)wd[p];
          if(diff<-0.3){ sw++; if(inr)sw_in++; swmap[p]=1;
            int nb=0; for(int dyy=-2;dyy<=2&&!nb;dyy++)for(int dxx=-2;dxx<=2&&!nb;dxx++){ int yz=y+dyy,xz=x+dxx; if(yz>=0&&yz<d&&xz>=0&&xz<d&&bar[(size_t)yz*d+xz])nb=1; }
            if(nb) sw_bar++; }
          else if(diff>1.5){ jumps++; swmap[p]=2; } }
        // LOCAL PITCH UNIFORMITY via outward rays (integer-crossing spacings -> CoV).
        double psum=0,psum2=0; long pn=0;
        for(int a=0;a<360;a++){ double th=a*M_PI/180.0, uy=sin(th),ux=cos(th);
          double pw=0,pr=-1; int havep=0;
          for(int k=2;k<Rr;k++){ int yy=(int)lround(cyf+uy*k),xx=(int)lround(cxf+ux*k); if(yy<0||yy>=d||xx<0||xx>=d)break;
            size_t q=(size_t)yy*d+xx; if(v[q]<athr||bar[q])continue; double w=wd[q];
            if(havep && (int)floor(pw)!=(int)floor(w) && fabs(w-pw)<1.5){ if(pr>=0){ double sp2=k-pr; if(sp2>1&&sp2<5*pitch){psum+=sp2;psum2+=sp2*sp2;pn++;} } pr=k; }
            pw=w; havep=1; } }
        printf("METRIC switch-rate       %.1f backward + %.1f wrap-jump  per 1000 material voxels  (%ld+%ld over %ld)\n",
               1000.0*sw/fmax(matc,1),1000.0*jumps/fmax(matc,1),sw,jumps,matc);
        printf("METRIC switch-locale     inner-third %.1f vs outer %.1f backward/1000  |  %.0f%% of switches at a barrier (touch)\n",
               1000.0*sw_in/fmax(mat_in,1), 1000.0*(sw-sw_in)/fmax(matc-mat_in,1), 100.0*sw_bar/fmax(sw,1));
        { u8 *hm=malloc(nn*3); for(size_t p=0;p<nn;p++){ int g=v[p]/5; hm[3*p]=hm[3*p+1]=hm[3*p+2]=(u8)g;
            if(bar[p]){ hm[3*p]=0;hm[3*p+1]=0;hm[3*p+2]=120; } }                       // barriers = dim blue
          for(size_t p=0;p<nn;p++){ if(swmap[p]==1){ hm[3*p]=255;hm[3*p+1]=0;hm[3*p+2]=0; }   // backward = red
              else if(swmap[p]==2){ hm[3*p]=255;hm[3*p+1]=180;hm[3*p+2]=0; } }               // jump = orange
          char hf[600]; snprintf(hf,sizeof hf,"%s_switch.ppm",base);
          FILE*hff=fopen(hf,"wb"); if(hff){ fprintf(hff,"P6\n%d %d\n255\n",d,d); fwrite(hm,1,nn*3,hff); fclose(hff); }
          free(hm); printf("wrote %s (switch heatmap: red=backward orange=jump, blue=barrier)\n",hf); }
        free(swmap);
        if(pn>3){ double pm=psum/pn, pv=psum2/pn-pm*pm; printf("METRIC pitch-uniformity  local-pitch mean=%.1f  CoV=%.2f  (n=%ld spacings)\n",pm,sqrt(fmax(0,pv))/fmax(pm,1e-6),pn); }

        // (c) PER-WRAP INTEGRITY: for each integer band, angular coverage (fraction of 360
        // it spans), fragment count (b0, connected pieces), and enclosed holes (b1). A clean
        // wrap is a single arc covering most of 360 with no extra loops; switches break it
        // into many fragments / spurious enclosed holes.
        u8 *bm=malloc(nn); u32 *bl=calloc(nn,sizeof(u32));
        double covs[4096]; int frags[4096],holes[4096]; int nbm=0,nbroken=0;
        for(int w=wlo; w<=whi && nbm<4096; w++){
          long cnt=0; int abins[72]={0};
          for(size_t p=0;p<nn;p++){ int in=(v[p]>=athr && (int)lround(wd[p])==w); bm[p]=(u8)in;
            if(in){cnt++; int yy=(int)(p/d),xx=(int)(p%d); double ang=atan2((double)yy-cyf,(double)xx-cxf); int bi=(int)((ang+M_PI)/(2*M_PI)*72); if(bi<0)bi=0; if(bi>71)bi=71; abins[bi]=1; } }
          if(cnt<150) continue;
          int abocc=0; for(int i=0;i<72;i++)abocc+=abins[i];
          // MORPHOLOGICAL CLOSE (dilate r then erode r): bridge intrinsic delamination air
          // gaps (a wrap is physically broken into arcs by ~24% internal air) so fragments/
          // holes measure FIELD continuity, not the scroll's delamination. A real sheet
          // switch displaces material by ~pitch and survives the close.
          { int R=3; u8*tb2=malloc(nn);
            for(int it2=0;it2<R;it2++){ memcpy(tb2,bm,nn);
              for(int yy=0;yy<d;yy++)for(int xx=0;xx<d;xx++){ size_t q=(size_t)yy*d+xx; if(tb2[q])continue;
                if((xx>0&&tb2[q-1])||(xx<d-1&&tb2[q+1])||(yy>0&&tb2[q-d])||(yy<d-1&&tb2[q+d])) bm[q]=1; } }
            for(int it2=0;it2<R;it2++){ memcpy(tb2,bm,nn);
              for(int yy=0;yy<d;yy++)for(int xx=0;xx<d;xx++){ size_t q=(size_t)yy*d+xx; if(!tb2[q])continue;
                if((xx>0&&!tb2[q-1])||(xx<d-1&&!tb2[q+1])||(yy>0&&!tb2[q-d])||(yy<d-1&&!tb2[q+d])) bm[q]=0; } }
            free(tb2); }
          u32 nc=cc_label(bm,1,d,d,TOPO_CONN26,bl); long*ca2=calloc((size_t)nc+1,sizeof(long));
          for(size_t p=0;p<nn;p++)ca2[bl[p]]++; int frag=0; for(u32 L=1;L<=nc;L++) if(ca2[L]>=20)frag++; free(ca2);
          for(size_t p=0;p<nn;p++) bm[p]=!bm[p];                        // complement for hole-counting
          u32 hc=cc_label(bm,1,d,d,TOPO_CONN6,bl); char*tb=calloc((size_t)hc+1,1);
          for(int x=0;x<d;x++){tb[bl[x]]=1; tb[bl[(size_t)(d-1)*d+x]]=1;} for(int y=0;y<d;y++){tb[bl[(size_t)y*d]]=1; tb[bl[(size_t)y*d+d-1]]=1;}
          long*ha=calloc((size_t)hc+1,sizeof(long)); for(size_t p=0;p<nn;p++)ha[bl[p]]++;
          int hole=0; for(u32 L=1;L<=hc;L++) if(!tb[L]&&ha[L]>=15)hole++; free(tb);free(ha);
          covs[nbm]=abocc/72.0; frags[nbm]=frag; holes[nbm]=hole;
          if(covs[nbm]<0.5) nbroken++; nbm++;   // broken = spans < half the revolution (frag count is delamination-dominated)
        }
        free(bm);free(bl);
        if(nbm>0){ for(int i=0;i<nbm;i++)for(int j=i+1;j<nbm;j++){
            if(covs[j]<covs[i]){double t=covs[i];covs[i]=covs[j];covs[j]=t;}
            if(frags[j]<frags[i]){int t=frags[i];frags[i]=frags[j];frags[j]=t;}
            if(holes[j]<holes[i]){int t=holes[i];holes[i]=holes[j];holes[j]=t;} }
          printf("METRIC per-wrap          bands=%d  median-coverage=%.0f%%  median-fragments=%d  median-holes=%d  broken=%d/%d\n",
                 nbm,100*covs[nbm/2],frags[nbm/2],holes[nbm/2],nbroken,nbm); }
      }
      PHASE("metrics");

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
