/* wind_poisson — INDEPENDENT winding solve for cross-validating sheet_sep3d (review P1).
 *
 * Solves the winding with a COMPLETELY DIFFERENT mechanism: the unwrap library's
 * winding_field_solve (a Poisson/Laplace diffusion warm-started from a contour-warped Archimedean
 * init), vs sheet_sep3d's radial-step-edge graph IRLS. Same archive, same region, same umbilicus,
 * same pitch -> if the two AGREE, the winding is method-independent (robust); if they diverge, the
 * disagreement localizes where one method is wrong. Dumps a _vol.f32 (same layout as sheet_sep3d's:
 * int32 {dz,dy,dx,z0,y0,x0} then dz*dy*dx float32, NaN off-material) so they compare directly.
 *
 * NOTE: winding_field produces a SPIRAL coordinate (r/pitch + theta/2pi); sheet_sep3d's is ~RADIAL
 * (constant along a sheet). To compare per-voxel, convert sheet_sep3d's w -> w + theta/2pi first
 * (scripts/compare_winding.py does this). The WRAP COUNT (span) is directly comparable.
 *
 *   wind_poisson ARCHIVE OUTBASE lod z0 y0 x0 dz dy dx pitch [iters=120] [anchor=0.02]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "io/mca.h"
#include "annotate/umbilicus.h"
#include "unwrap/winding_field.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"

static int air_threshold(const u8 *v, size_t n){
  long h[256]={0}; long tot=0; double allsum=0; long allnz=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x){ allsum+=x; allnz++; } if(x>=1&&x<=254){ h[x]++; tot++; } }
  if(tot<256){ int t=allnz?(int)(0.5*allsum/allnz+0.5):1; if(t<1)t=1; if(t>254)t=254; return t; }
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sumB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){ wB+=h[t]; if(wB==0)continue; double wF=(double)tot-wB; if(wF<=0)break;
    sumB+=(double)t*h[t]; double mB=sumB/wB,mF=(sum-sumB)/wF,btw=wB*wF*(mB-mF)*(mB-mF);
    if(btw>best){best=btw;thr=t;} }
  return thr;
}

// Tiled structure-tensor detection: st_sheet_detect allocates 8 full-volume f32 buffers internally
// (~8.8GB at dz96) -> the memory peak. The tensor is LOCAL (depends only on a ~3*sigma neighborhood),
// so compute it in 3D blocks with a halo and stitch the interiors. Bounds ST memory to one block.
// Result is identical to the whole-volume call except within `halo` of the global boundary (where the
// blurs see a truncated stencil either way). v=u8 CT (no full f32 copy needed).
static void st_tiled(const u8 *v, int dz, int dy, int dx, const st_params *sp,
                     f32 *sheet, f32 *nrm, int bz, int by, int bx) {
  double smax = sp->sigma_tensor > sp->sigma_grad ? sp->sigma_tensor : sp->sigma_grad;
  int H = (int)ceil(3.0 * smax) + 2;
  for (int z0t = 0; z0t < dz; z0t += bz)
    for (int y0t = 0; y0t < dy; y0t += by)
      for (int x0t = 0; x0t < dx; x0t += bx) {
        int z1t = z0t+bz<dz?z0t+bz:dz, y1t=y0t+by<dy?y0t+by:dy, x1t=x0t+bx<dx?x0t+bx:dx;
        int ez0=z0t-H<0?0:z0t-H, ey0=y0t-H<0?0:y0t-H, ex0=x0t-H<0?0:x0t-H;
        int ez1=z1t+H>dz?dz:z1t+H, ey1=y1t+H>dy?dy:y1t+H, ex1=x1t+H>dx?dx:x1t+H;
        int ez=ez1-ez0, ey=ey1-ey0, ex=ex1-ex0;
        size_t en=(size_t)ez*ey*ex;
        f32 *sub=malloc(en*sizeof(f32)), *ssh=malloc(en*sizeof(f32)), *snr=malloc(3*en*sizeof(f32));
        for(int z=0;z<ez;z++)for(int y=0;y<ey;y++)for(int x=0;x<ex;x++)
          sub[((size_t)z*ey+y)*ex+x]=v[((size_t)(z+ez0)*dy+(y+ey0))*dx+(x+ex0)];
        st_sheet_detect(sub,ez,ey,ex,sp,ssh,snr);
        for(int z=z0t;z<z1t;z++)for(int y=y0t;y<y1t;y++)for(int x=x0t;x<x1t;x++){
          size_t gi=((size_t)z*dy+y)*dx+x;
          size_t li=((size_t)(z-ez0)*ey+(y-ey0))*ex+(x-ex0);
          sheet[gi]=ssh[li]; nrm[3*gi]=snr[3*li]; nrm[3*gi+1]=snr[3*li+1]; nrm[3*gi+2]=snr[3*li+2];
        }
        free(sub);free(ssh);free(snr);
      }
}

// ridge voxel (one z-slice) for intrinsic radial-count winding; sorted by angular bin then radius.
typedef struct { float r, th; int li, bin; } rv_t;
static int rv_cmp(const void *a, const void *b){
  const rv_t *p=a,*q=b; if(p->bin!=q->bin) return p->bin<q->bin?-1:1; return p->r<q->r?-1:p->r>q->r?1:0;
}

// TRUE-3D global normal orientation (replaces the 2.5D umbilicus-outward flip).
// The structure-tensor sheet normal is sign-ambiguous (+-n). For the winding gradient n/pitch to
// increment CONSISTENTLY across sheets we need n oriented the same spiral way everywhere -- in 3D,
// following the real sheet geometry (incl. concave fold-backs the center-outward assumption rounds
// off). Method = confidence-ordered region grow (Hoppe-style): visit material voxels in DECREASING
// sheetness (most reliable normals first); flip each to agree (positive dot) with the RESULTANT of
// its already-oriented 6-neighbors. Sheet interiors fix first and propagate outward; adjacent wraps'
// outward normals agree across the (touching, at LOD2) gap, so orientation crosses sheets coherently.
// No umbilicus, no azimuth, no per-slice anything. O(N) via a 256-bin radix on sheetness.
static long orient_normals_3d(const u8 *mask, f32 *nrm, const f32 *sheet, int dz, int dy, int dx){
  size_t nn=(size_t)dz*dy*dx, nynx=(size_t)dy*dx; const int NB=256;
  u8 *st=calloc(nn,1); long flips=0, ncomp=0;   // st: 0 unseen, 1 queued, 2 done
  // Weak outward-from-centroid bias on the orientation vote: pure local propagation has NO global
  // anchor, so over the full cross-section the front develops sign SEAMS (W wanders +-200 instead of
  // 0..~130). The bias only breaks TIES (where neighbor agreement ~0 = the seams); a real concave
  // fold's neighbor resultant (magnitude ~few) overrides the small bias, so folds are preserved.
  double cz=0,cy=0,cx=0,cc=0;
  for(size_t p=0;p<nn;p++) if(mask[p]){ int z=(int)(p/nynx),r=(int)(p%nynx); cz+=z; cy+=r/dx; cx+=r%dx; cc++; }
  if(cc>0){ cz/=cc; cy/=cc; cx/=cc; }
  double bias=0.5; { const char*e=getenv("ORIENT_BIAS"); if(e)bias=atof(e); }
  // SINGLE-FRONT priority region grow: grow ONE connected front (high-sheetness frontier first) so the
  // WHOLE connected component gets ONE consistent global orientation (the per-voxel sheetness-order grow
  // seeded disconnected regions with arbitrary signs -> sign-domains -> globally non-monotone W of +-200
  // wraps). 256 sheetness buckets of frontier voxels = O(N) priority queue. Only voxels ADJACENT to the
  // already-oriented set are ever oriented, so orientation propagates contiguously from a single seed.
  size_t *bcnt=calloc(NB,sizeof(size_t)), *bcap=calloc(NB,sizeof(size_t)), **buk=calloc(NB,sizeof(size_t*));
  #define SBIN(p) ({ int _b=(int)(sheet[p]*255.0f); _b<0?0:(_b>255?255:_b); })
  #define PUSH(p) do{ int _b=SBIN(p); if(bcnt[_b]==bcap[_b]){ bcap[_b]=bcap[_b]?bcap[_b]*2:64; buk[_b]=realloc(buk[_b],bcap[_b]*sizeof(size_t)); } buk[_b][bcnt[_b]++]=(p); st[p]=1; if(_b>top)top=_b; }while(0)
  int top=-1;
  // material voxels radix-sorted by DESCENDING sheetness, ONCE -> O(N) component seeding (a pointer
  // advances through it; avoids an O(nn) rescan per disconnected speckle, which hung at L2).
  size_t *cnt=calloc(NB,sizeof(size_t)), tot=0;
  for(size_t p=0;p<nn;p++) if(mask[p]){ cnt[SBIN(p)]++; tot++; }
  size_t *posb=malloc(NB*sizeof(size_t)), accb=0; for(int b=NB-1;b>=0;b--){ posb[b]=accb; accb+=cnt[b]; }
  size_t *order=malloc((tot?tot:1)*sizeof(size_t));
  for(size_t p=0;p<nn;p++) if(mask[p]) order[posb[SBIN(p)]++]=p;
  free(cnt); free(posb);
  if(tot==0){ free(st);free(bcnt);free(bcap);free(buk);free(order); return 0; }
  size_t op=0; PUSH(order[0]); ncomp++;
  for(;;){
    int b=top; while(b>=0 && bcnt[b]==0) b--; top=b;
    if(b<0){ // frontier empty: advance the sorted pointer to the next unseen voxel = next component seed
      while(op<tot && st[order[op]]!=0) op++;
      if(op>=tot) break; PUSH(order[op]); ncomp++; continue; }
    size_t p=buk[b][--bcnt[b]]; if(st[p]==2) continue;
    int z=(int)(p/nynx), rem=(int)(p%nynx), y=rem/dx, x=rem%dx;
    double nx=nrm[3*p], ny=nrm[3*p+1], nz=nrm[3*p+2], vx=0,vy=0,vz=0;
    size_t nb[6]; int nv=0;
    if(x>0)nb[nv++]=p-1; if(x<dx-1)nb[nv++]=p+1;
    if(y>0)nb[nv++]=p-dx; if(y<dy-1)nb[nv++]=p+dx;
    if(z>0)nb[nv++]=p-nynx; if(z<dz-1)nb[nv++]=p+nynx;
    for(int i=0;i<nv;i++){ size_t q=nb[i]; if(mask[q]&&st[q]==2){ vx+=nrm[3*q]; vy+=nrm[3*q+1]; vz+=nrm[3*q+2]; } }
    if(bias>0){ double ox=x-cx,oy=y-cy,oz=z-cz,ol=sqrt(ox*ox+oy*oy+oz*oz); if(ol>1e-6){ vx+=bias*ox/ol; vy+=bias*oy/ol; vz+=bias*oz/ol; } }
    if(vx*nx+vy*ny+vz*nz<0){ nrm[3*p]=(f32)-nx; nrm[3*p+1]=(f32)-ny; nrm[3*p+2]=(f32)-nz; flips++; }
    st[p]=2;
    for(int i=0;i<nv;i++){ size_t q=nb[i]; if(mask[q]&&st[q]==0) PUSH(q); }
  }
  for(int b=0;b<NB;b++) free(buk[b]);
  free(bcnt);free(bcap);free(buk);free(order);
  fprintf(stderr,"3D orient: single-front grow, %ld components\n",ncomp);
  u8 *done=st;   // alias for the relaxation below (st==2 everywhere material)
  #undef SBIN
  #undef PUSH
  // RELAXATION: the single greedy pass can lock early errors (wrong flip propagates). Sweep: flip any
  // voxel whose normal disagrees (negative resultant dot) with its already-oriented neighbors. Gauss-
  // Seidel, a few passes -> isolated orientation errors self-correct, cutting the field's overcount jumps.
  int relax=0; { const char*e=getenv("ORIENT_RELAX"); if(e)relax=atoi(e); }
  for(int it=0; it<relax; it++){ long ch=0;
    for(size_t p=0;p<nn;p++){ if(!mask[p])continue; int z=(int)(p/nynx),rem=(int)(p%nynx),y=rem/dx,x=rem%dx;
      double nx=nrm[3*p],ny=nrm[3*p+1],nz=nrm[3*p+2],vx=0,vy=0,vz=0;
      if(x>0&&mask[p-1]){vx+=nrm[3*(p-1)];vy+=nrm[3*(p-1)+1];vz+=nrm[3*(p-1)+2];}
      if(x<dx-1&&mask[p+1]){vx+=nrm[3*(p+1)];vy+=nrm[3*(p+1)+1];vz+=nrm[3*(p+1)+2];}
      if(y>0&&mask[p-dx]){vx+=nrm[3*(p-dx)];vy+=nrm[3*(p-dx)+1];vz+=nrm[3*(p-dx)+2];}
      if(y<dy-1&&mask[p+dx]){vx+=nrm[3*(p+dx)];vy+=nrm[3*(p+dx)+1];vz+=nrm[3*(p+dx)+2];}
      if(z>0&&mask[p-nynx]){vx+=nrm[3*(p-nynx)];vy+=nrm[3*(p-nynx)+1];vz+=nrm[3*(p-nynx)+2];}
      if(z<dz-1&&mask[p+nynx]){vx+=nrm[3*(p+nynx)];vy+=nrm[3*(p+nynx)+1];vz+=nrm[3*(p+nynx)+2];}
      if(vx*nx+vy*ny+vz*nz<0){ nrm[3*p]=(f32)-nx; nrm[3*p+1]=(f32)-ny; nrm[3*p+2]=(f32)-nz; ch++; } }
    flips+=ch; if(!ch)break;
  }
  free(done);
  return flips;
}

// TRUE-3D local sheet PITCH from the 3D ridges (replaces the uniform-pitch prior). At each ridge
// voxel march along +-n (unit normal) to the next ridge -> the local sheet spacing in voxels. Then
// multi-source BFS propagates each material voxel the pitch of its nearest ridge. forcing=n/pitch_local
// then increments W by ~1 across EVERY sheet regardless of how the spacing varies (folds vs open) ->
// tightens the increment distribution the uniform pitch leaves spread. Fully 3D (marches in 3D, BFS in 3D).
// RIDGE DEFRAGMENTATION (tangent-plane closing): ridge_nms ridges are porous (gaps ALONG sheets) ->
// the pitch march skips to every-OTHER sheet (measures ~2x) and surfaces fragment. Grow each ridge
// voxel ONLY within its own sheet plane (neighbor offset ~perpendicular to the normal), iterated, so
// coplanar fragments connect ACROSS gaps but the ridge never thickens ACROSS the normal (adjacent
// sheets stay separate). Fully 3D. Returns voxels added.
static long ridge_defrag_tangent(const u8 *mask, u8 *ridge, const f32 *nrm, int dz, int dy, int dx,
                                 int R, double tthr){
  // FILL HOLES ONLY: mark a non-ridge voxel as ridge iff it sits in a GAP of a sheet -- ridge present
  // in OPPOSITE in-plane directions (so a hole gets bridged, but sheet boundaries, ridge on one side
  // only, do NOT grow). Tangent basis t1,t2 from the voxel's own normal; probe +-t1, +-t2 out to R.
  size_t nn=(size_t)dz*dy*dx, nynx=(size_t)dy*dx; u8 *add=calloc(nn,1); long nadd=0; (void)tthr;
  #pragma omp parallel for schedule(dynamic,8) reduction(+:nadd)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*nynx+(size_t)y*dx+x;
    if(ridge[p]||!mask[p])continue;
    double nx=nrm[3*p],ny=nrm[3*p+1],nz=nrm[3*p+2], nl=sqrt(nx*nx+ny*ny+nz*nz);
    if(nl<1e-6)continue; nx/=nl;ny/=nl;nz/=nl;
    double ex= fabs(nx)<0.9?1:0, ey= fabs(nx)<0.9?0:1, ez=0;     // a vector not parallel to n
    double t1x=ny*ez-nz*ey, t1y=nz*ex-nx*ez, t1z=nx*ey-ny*ex, t1l=sqrt(t1x*t1x+t1y*t1y+t1z*t1z);
    if(t1l<1e-6)continue; t1x/=t1l;t1y/=t1l;t1z/=t1l;
    double t2x=ny*t1z-nz*t1y, t2y=nz*t1x-nx*t1z, t2z=nx*t1y-ny*t1x;   // n x t1 (unit)
    int hp1=0,hm1=0,hp2=0,hm2=0;
    for(int r=1;r<=R;r++){
      #define RDG(ax,ay,az) ({ int zz=(int)lround(z+(az)),yy=(int)lround(y+(ay)),xx=(int)lround(x+(ax)); \
        (zz>=0&&zz<dz&&yy>=0&&yy<dy&&xx>=0&&xx<dx&&ridge[(size_t)zz*nynx+(size_t)yy*dx+xx]&&mask[(size_t)zz*nynx+(size_t)yy*dx+xx]); })
      if(!hp1&&RDG(r*t1x,r*t1y,r*t1z))hp1=1; if(!hm1&&RDG(-r*t1x,-r*t1y,-r*t1z))hm1=1;
      if(!hp2&&RDG(r*t2x,r*t2y,r*t2z))hp2=1; if(!hm2&&RDG(-r*t2x,-r*t2y,-r*t2z))hm2=1;
      #undef RDG
    }
    if((hp1&&hm1)||(hp2&&hm2)){ add[p]=1; nadd++; }
  }
  for(size_t p=0;p<nn;p++) if(add[p]) ridge[p]=1;
  free(add); return nadd;
}

// GLOBAL-LOCAL pitch: real sheet pitch is variable but SMOOTH (fairly consistent regionally), so a
// per-voxel ridge-march is the wrong tool (porosity makes it noisy/skip sheets). Instead aggregate to a
// COARSE BLOCK GRID: per block take the MEDIAN ridge-spacing (robust to porosity/outliers), fill empty
// blocks + smooth the coarse grid, then TRILINEAR-interpolate back to a smooth per-voxel pitch. The
// march itself is porosity-tolerant (3x3x3 check at each step bridges gaps in the adjacent ridge).
static void measure_pitch_3d(const u8 *mask, const u8 *ridge, const f32 *nrm, int dz, int dy, int dx,
                             double pitch0, double minp, double maxp, f32 *pitchf){
  size_t nynx=(size_t)dy*dx;
  int BS=48; { const char*e=getenv("VARPITCH_BLOCK"); if(e)BS=atoi(e); }
  int nbz=(dz+BS-1)/BS, nby=(dy+BS-1)/BS, nbx=(dx+BS-1)/BS; size_t nb=(size_t)nbz*nby*nbx;
  const int HB=40; int *hist=calloc(nb*HB,sizeof(int));
  int MAXT=(int)(maxp*1.6)+2; long nseed=0;
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*nynx+(size_t)y*dx+x;
    if(!ridge[p]||!mask[p])continue;
    double nx=nrm[3*p],ny=nrm[3*p+1],nz=nrm[3*p+2]; int g=0;
    for(int s=1;s>=-1&&!g;s-=2){
      int phase=0, gapc=0;                            // phase0: still on own sheet; phase1: in gap, seek next ridge
      for(int t=1;t<=MAXT;t++){ int zz=(int)lround(z+s*t*nz),yy=(int)lround(y+s*t*ny),xx=(int)lround(x+s*t*nx);
        if(zz<1||zz>=dz-1||yy<1||yy>=dy-1||xx<1||xx>=dx-1)break;
        int isr=0;                                    // 3x3x3 tolerance: bridge porous gaps in the ridge
        for(int a=-1;a<=1&&!isr;a++)for(int b=-1;b<=1&&!isr;b++)for(int c=-1;c<=1&&!isr;c++){
          size_t pp=(size_t)(zz+a)*nynx+(size_t)(yy+b)*dx+(xx+c); if(ridge[pp]&&mask[pp])isr=1; }
        if(phase==0){ if(!isr){ if(++gapc>=2) phase=1; } else gapc=0; }   // leave own sheet after 2 clear steps
        else if(isr){ g=t; break; }                   // first ridge AFTER the gap = the adjacent sheet
      } }
    if(g>=(int)minp&&g<=(int)maxp){ size_t bi=((size_t)(z/BS)*nby+(y/BS))*nbx+(x/BS); int gi=g<HB?g:HB-1; hist[bi*HB+gi]++; nseed++; }
  }
  // per-block median pitch (robust); blocks with too few ridge samples = unknown (filled by diffusion)
  float *bp=malloc(nb*sizeof(float)); u8 *bok=calloc(nb,1);
  for(size_t bi=0;bi<nb;bi++){ int tot=0; for(int h=0;h<HB;h++) tot+=hist[bi*HB+h];
    if(tot<10){ bp[bi]=(float)pitch0; continue; }
    int cum=0,med=(int)pitch0; for(int h=0;h<HB;h++){ cum+=hist[bi*HB+h]; if(cum*2>=tot){med=h;break;} }
    bp[bi]=(float)med; bok[bi]=1; }
  free(hist);
  // diffuse into empty blocks + smooth the coarse grid (Jacobi; known blocks pinned-ish via heavy self-weight)
  float *tb=malloc(nb*sizeof(float));
  int bsm=25; { const char*e=getenv("VARPITCH_BSMOOTH"); if(e)bsm=atoi(e); }   // few passes: fill empties, KEEP regional variation
  for(int it=0;it<bsm;it++){
    for(int bz=0;bz<nbz;bz++)for(int by=0;by<nby;by++)for(int bx=0;bx<nbx;bx++){ size_t bi=((size_t)bz*nby+by)*nbx+bx;
      double s=bp[bi]*(bok[bi]?4.0:0.0); double c=bok[bi]?4.0:0.0;
      if(bx>0){s+=bp[bi-1];c++;} if(bx<nbx-1){s+=bp[bi+1];c++;}
      if(by>0){s+=bp[bi-nbx];c++;} if(by<nby-1){s+=bp[bi+nbx];c++;}
      if(bz>0){s+=bp[bi-(size_t)nby*nbx];c++;} if(bz<nbz-1){s+=bp[bi+(size_t)nby*nbx];c++;}
      tb[bi]=(float)(c>0?s/c:pitch0); }
    memcpy(bp,tb,nb*sizeof(float)); }
  free(tb); free(bok);
  // trilinear interpolate coarse block-center grid -> smooth per-voxel pitch (clamped)
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*nynx+(size_t)y*dx+x;
    double fz=((double)z+0.5)/BS-0.5, fy=((double)y+0.5)/BS-0.5, fx=((double)x+0.5)/BS-0.5;
    int bz=(int)floor(fz),by=(int)floor(fy),bx=(int)floor(fx); double tz=fz-bz,ty=fy-by,tx=fx-bx;
    double val=0,wsum=0;
    for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){ int Z=bz+a,Y=by+b,X=bx+c;
      if(Z<0||Z>=nbz||Y<0||Y>=nby||X<0||X>=nbx)continue;
      double w=(a?tz:1-tz)*(b?ty:1-ty)*(c?tx:1-tx); val+=w*bp[((size_t)Z*nby+Y)*nbx+X]; wsum+=w; }
    double pv = wsum>0? val/wsum : pitch0; if(pv<minp)pv=minp; if(pv>maxp)pv=maxp;
    pitchf[p]=(f32)pv; }
  // report the pitch field distribution (so the "global-local" value is visible)
  double psum=0; long pc=0; float pmin=1e9,pmax=-1e9;
  for(size_t p=0;p<(size_t)dz*nynx;p++) if(mask[p]){ psum+=pitchf[p]; pc++; if(pitchf[p]<pmin)pmin=pitchf[p]; if(pitchf[p]>pmax)pmax=pitchf[p]; }
  free(bp);
  fprintf(stderr,"VARPITCH global-local: %ld ridge spacings, block=%d -> pitch %.1f..%.1f vox (mean %.2f)\n",
          nseed,BS,pmin,pmax,pc?psum/pc:0);
}

int main(int argc,char**argv){
  if(argc<11){ fprintf(stderr,"usage: %s ARC OUTBASE lod z0 y0 x0 dz dy dx pitch [iters=120] [anchor=0.02]\n",argv[0]); return 2; }
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]);
  int dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]); double pitch=atof(argv[10]);
  int iters=argc>11?atoi(argv[11]):120; double anchor=argc>12?atof(argv[12]):0.02;
  int usenormals=argc>13?atoi(argv[13]):1;   // 3D-DEFAULT ON: forcing = div(oriented 3D sheet-normal/pitch)

  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  // global umbilicus from coarse LOD5 (same as sheet_sep3d)
  int cl=5; double cs=(double)(1<<cl); int ccz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,ccz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse fail\n");return 1;}
  size_t ccn=(size_t)ccz*ccy*ccx; u8*ccm=malloc(ccn); for(size_t i=0;i<ccn;i++)ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,ccz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;} free(ccm);
  double ls=(double)(1<<lod);
  // subvol-local umbilicus: sample the global axis at NC control z, convert LOD5->subvol-local
  const int NC=9; umbilicus ul; ul.n=NC; ul.z=malloc(NC*sizeof(f32)); ul.y=malloc(NC*sizeof(f32)); ul.x=malloc(NC*sizeof(f32));
  for(int i=0;i<NC;i++){ double zz=(double)i*(dz-1)/(NC-1); double coarse_z=(z0+zz)*ls/cs;
    if(coarse_z<0)coarse_z=0; if(coarse_z>ccz-1)coarse_z=ccz-1; f32 ucy,ucx; umbilicus_center(&umb,(f32)coarse_z,&ucy,&ucx);
    ul.z[i]=(f32)zz; ul.y[i]=(f32)(ucy*cs/ls - y0); ul.x[i]=(f32)(ucx*cs/ls - x0); }
  fprintf(stderr,"subvol z%ld y%ld x%ld + %dx%dx%d; center z0=(%.0f,%.0f) zmid=(%.0f,%.0f) pitch=%.1f\n",
          z0,y0,x0,dz,dy,dx,ul.y[0],ul.x[0],ul.y[NC/2],ul.x[NC/2],pitch);

  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,dz,dy,dx); if(!v){fprintf(stderr,"subvol read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr?1:0; nmat+=mask[p]; }
  fprintf(stderr,"air<%d, %ld/%zu material\n",athr,nmat,nn);

  // P2 NORMAL FORCING: instead of the uniform-pitch radial assumption, drive the Poisson solve to
  // match the ACTUAL sheet geometry. Desired winding gradient g = n_out/pitch (n_out = sheet normal
  // oriented outward from the umbilicus); forcing = div(g) turns the Laplacian into a Poisson solve
  // whose level sets follow the real sheets, capturing variable density a uniform-pitch field misses.
  // ANISO (env, e.g. 0.1): anisotropic along-sheet smoothing. Cross-sheet diffusivity alpha in
  // (0,1]; lower = sharper across sheets. 0/unset disables. Built from the UNORIENTED structure-
  // tensor sheet normal -> no umbilicus-outward orientation (which fails in folds). The structure
  // tensor is computed once and shared with the (optional) normal forcing.
  double aniso_alpha=0.05; { const char*e=getenv("ANISO"); if(e)aniso_alpha=atof(e); }   // 3D-DEFAULT: along-sheet anisotropy
  int force_iter=6; { const char*e=getenv("FORCE_ITER"); if(e)force_iter=atoi(e); if(force_iter<1)force_iter=1; }   // 3D-DEFAULT: propagate + gradW self-orient passes converge rate->1.0
  // PLANARITY weighting (env PLANARITY=1, default on): the structure tensor's sheetness in [0,1]
  // = (l0-l1)/(l0+eps) is a per-voxel confidence in the normal. Where it's low (ambiguous/noisy),
  // down-weight BOTH the anisotropy (-> isotropic) and the forcing (-> let smoothing rule), so bad
  // normals stop injecting wrong gradients. This directly attacks the structure-tensor noise that
  // caps the metric.
  int planarity=0; { const char*e=getenv("PLANARITY"); if(e)planarity=atoi(e); }   // default OFF: down-weighting noisy normals empirically HURT (even noisy normals carry useful avg direction; weakening anisotropy lets the isotropic Laplacian smear back in)
  int tensor6_on=0; { const char*e=getenv("TENSOR6"); if(e)tensor6_on=atoi(e); }
  int wind_anchor=0; { const char*e=getenv("WIND_ANCHOR"); if(e)wind_anchor=atoi(e); }   // Phase 3: ridge-anchored re-solve
  int varpitch=0; { const char*e=getenv("VARPITCH"); if(e)varpitch=atoi(e); }   // 3D local pitch from ridge spacing
  f32*forcing=NULL; f32*aniso=NULL; f32*tensor6=NULL; f32*nrm=NULL; f32*sheet=NULL; u8*ridge=NULL; f32*pitchf=NULL;
  if(usenormals || aniso_alpha>0.0){
    nrm=malloc(3*nn*sizeof(f32)); sheet=malloc(nn*sizeof(f32));
    // 3D-DEFAULT sigma_tensor ~0.5*pitch: empirically smoother normals (less curl) lift the wrap-count
    // rate toward 1.0 and tighten its spread; clamp [1.5,4]. (was pitch*0.12 -> too noisy.)
    st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=(f32)(pitch*0.5); if(sp.sigma_tensor<1.5f)sp.sigma_tensor=1.5f; if(sp.sigma_tensor>4.0f)sp.sigma_tensor=4.0f;
    { const char*e=getenv("ST_SIGMA"); if(e){ sp.sigma_tensor=(f32)atof(e); fprintf(stderr,"ST_SIGMA override -> sigma_tensor=%.2f\n",sp.sigma_tensor); } }
    { const char*e=getenv("ST_SGRAD"); if(e){ sp.sigma_grad=(f32)atof(e); } }
    // tiled ST (env ST_TILE=N block size, default 64; bounds the 8-buffer ST peak to one block)
    int stb=64; { const char*e=getenv("ST_TILE"); if(e)stb=atoi(e); }
    if(stb>0){ st_tiled(v,dz,dy,dx,&sp,sheet,nrm,stb,stb,stb); fprintf(stderr,"tiled ST (block=%d, halo-stitched)\n",stb); }
    else { f32*vf=malloc(nn*sizeof(f32)); for(size_t q=0;q<nn;q++)vf[q]=v[q]; st_sheet_detect(vf,dz,dy,dx,&sp,sheet,nrm); free(vf); }
    // ST_MULTISCALE=f: also detect at sigma_tensor*f and, per voxel, KEEP the normal from whichever
    // scale shows the clearer sheet (higher sheetness). Denoises by SELECTION (vs down-weighting,
    // which hurt): a sheet that's blurred at one scale is often crisp at another.
    { const char*e=getenv("ST_MULTISCALE"); double f= e?atof(e):0.0;
      if(f>0.0){ st_params sp2=sp; sp2.sigma_tensor=(f32)(sp.sigma_tensor*f);
        f32*nrm2=malloc(3*nn*sizeof(f32)),*sh2=malloc(nn*sizeof(f32));
        if(stb>0) st_tiled(v,dz,dy,dx,&sp2,sh2,nrm2,stb,stb,stb);
        else { f32*vf=malloc(nn*sizeof(f32)); for(size_t q=0;q<nn;q++)vf[q]=v[q]; st_sheet_detect(vf,dz,dy,dx,&sp2,sh2,nrm2); free(vf); }
        long swapped=0;
        for(size_t p=0;p<nn;p++) if(sh2[p]>sheet[p]){ sheet[p]=sh2[p]; nrm[3*p]=nrm2[3*p]; nrm[3*p+1]=nrm2[3*p+1]; nrm[3*p+2]=nrm2[3*p+2]; swapped++; }
        free(nrm2);free(sh2);
        fprintf(stderr,"ST_MULTISCALE sigma*%.2f: %ld/%zu voxels took the finer/coarser-scale normal\n",f,swapped,nn); }
    }
    // PHASE 2 (env WIND_RIDGE=base): detect sheet CENTERLINE ridges (NMS along the normal) from the
    // REAL sheetness+normal (before the planarity reset below). These are the actual sheets the Phase-3
    // winding anchoring will pin to. Dump as a u8 mask (_vol-style header) for validation vs GT.
    { const char*rb=getenv("WIND_RIDGE");
      if(rb || wind_anchor || varpitch){ f32*vf2=malloc(nn*sizeof(f32)); for(size_t q=0;q<nn;q++)vf2[q]=v[q];
        ridge=malloc(nn);
        f32 smin=0.25f,imin=70.f,step=1.0f; { const char*e=getenv("RIDGE_SMIN"); if(e)smin=(f32)atof(e); } { const char*e=getenv("RIDGE_IMIN"); if(e)imin=(f32)atof(e); }
        size_t nr=ridge_nms(vf2,sheet,nrm,dz,dy,dx,smin,imin,step,ridge); free(vf2);
        int defrag=0; { const char*e=getenv("DEFRAG"); if(e)defrag=atoi(e); }   // tangent hole-fill (off: L2 ridges MISS sheets, not holey; useful at L1)
        if(defrag>0){ double tthr=0.5; { const char*e=getenv("DEFRAG_THR"); if(e)tthr=atof(e); }
          long added=ridge_defrag_tangent(mask,ridge,nrm,dz,dy,dx,defrag,tthr);
          fprintf(stderr,"ridge defrag: +%ld voxels over %d tangent-close iters (thr=%.2f)\n",added,defrag,tthr); }
        if(rb){ char rf[760]; snprintf(rf,sizeof rf,"%s_ridge.u8",rb); FILE*rfp=fopen(rf,"wb");
          if(rfp){ int hdr[6]={dz,dy,dx,(int)z0,(int)y0,(int)x0}; fwrite(hdr,sizeof(int),6,rfp); fwrite(ridge,1,nn,rfp); fclose(rfp); } }
        fprintf(stderr,"PHASE2 ridge: %zu ridge voxels (smin=%.2f imin=%.0f)\n",nr,smin,imin);
        if(!wind_anchor && !varpitch){ free(ridge); ridge=NULL; }
      } }
    if(!planarity){ for(size_t p=0;p<nn;p++) sheet[p]=1.0f; }   // disable weighting
    if(aniso_alpha>0.0 && tensor6_on==1){
      // ON-THE-FLY full tensor: pass the (already-resident) normal+sheetness; the solver builds D per
      // voxel -> NO 6*N array (saves ~6.6GB at dz96). Set below via wp.normal/sheetness/tensor_alpha.
      fprintf(stderr,"FULL anisotropic TENSOR (off-diagonal, ON-THE-FLY, no array) ON (alpha=%.3f)\n",aniso_alpha);
    } else if(aniso_alpha>0.0 && tensor6_on>=2){
      // FULL diffusion tensor D = I-(1-a)*s*n n^T, 6 components/voxel [Dzz,Dyy,Dxx,Dyz,Dxz,Dxy] (explicit array; for equivalence checks).
      tensor6=malloc(6*nn*sizeof(f32)); double a=aniso_alpha;
      for(size_t p=0;p<nn;p++){
        if(!mask[p]){ tensor6[6*p]=tensor6[6*p+1]=tensor6[6*p+2]=1.0f; tensor6[6*p+3]=tensor6[6*p+4]=tensor6[6*p+5]=0.0f; continue; }
        double nx=nrm[3*p+0],ny=nrm[3*p+1],nz=nrm[3*p+2], nn2=nx*nx+ny*ny+nz*nz;
        double s=sheet[p]; if(s<0)s=0; if(s>1)s=1; double k=(1.0-a)*s;
        if(nn2<1e-6){ tensor6[6*p]=tensor6[6*p+1]=tensor6[6*p+2]=1.0f; tensor6[6*p+3]=tensor6[6*p+4]=tensor6[6*p+5]=0.0f; continue; }
        tensor6[6*p+0]=(f32)(1.0-k*nz*nz);  // Dzz
        tensor6[6*p+1]=(f32)(1.0-k*ny*ny);  // Dyy
        tensor6[6*p+2]=(f32)(1.0-k*nx*nx);  // Dxx
        tensor6[6*p+3]=(f32)(-k*ny*nz);     // Dyz
        tensor6[6*p+4]=(f32)(-k*nx*nz);     // Dxz
        tensor6[6*p+5]=(f32)(-k*nx*ny);     // Dxy
      }
      fprintf(stderr,"FULL anisotropic TENSOR (off-diagonal) ON (alpha=%.3f)\n",a);
    } else if(aniso_alpha>0.0){
      // per-voxel axis diffusivity w_axis = 1-(1-alpha)*s*n_axis^2 (s=sheetness); low sheetness -> ~isotropic.
      aniso=malloc(3*nn*sizeof(f32)); double a=aniso_alpha;
      for(size_t p=0;p<nn;p++){
        if(!mask[p]){ aniso[3*p]=aniso[3*p+1]=aniso[3*p+2]=1.0f; continue; }
        double nx=nrm[3*p+0],ny=nrm[3*p+1],nz=nrm[3*p+2], nn2=nx*nx+ny*ny+nz*nz;
        double s=sheet[p]; if(s<0)s=0; if(s>1)s=1;
        if(nn2<1e-6){ aniso[3*p]=aniso[3*p+1]=aniso[3*p+2]=1.0f; continue; }   // isotropic where no sheet
        aniso[3*p+0]=(f32)(1.0-(1.0-a)*s*nz*nz);  // wz
        aniso[3*p+1]=(f32)(1.0-(1.0-a)*s*ny*ny);  // wy
        aniso[3*p+2]=(f32)(1.0-(1.0-a)*s*nx*nx);  // wx
      }
      fprintf(stderr,"ANISOTROPIC along-sheet smoothing ON (alpha=%.3f, planarity-weighted=%d)\n",a,planarity);
    }
  }
  // ---- TRUE-3D init: NO contour-warp (per-slice), NO umbilicus-polar theta (2.5D azimuth). ----
  // Orient the 3D structure-tensor normals globally (geometry-following), then warm-start the solve
  // from a coarse 3D distance-from-material-centroid / pitch (a monotone outward guess; the LINEAR
  // solve is init-independent per pass, the forcing+oriented-normals determine the result). One
  // Dirichlet pin at the centroid fixes the global constant (no anchor needed).
  f32*wind=malloc(nn*sizeof(f32));
  if(usenormals && nrm){
    long fl=orient_normals_3d(mask,nrm,sheet,dz,dy,dx);
    fprintf(stderr,"3D normal orientation: %ld flips over %ld material\n",fl,nmat);
  }
  if(varpitch && ridge && nrm){
    double minp=pitch*0.4, maxp=pitch*2.5; { const char*e=getenv("VARPITCH_MIN"); if(e)minp=atof(e); } { const char*e=getenv("VARPITCH_MAX"); if(e)maxp=atof(e); }
    pitchf=malloc(nn*sizeof(f32)); measure_pitch_3d(mask,ridge,nrm,dz,dy,dx,pitch,minp,maxp,pitchf);
  }
  double cz=0,cy=0,cx=0; { double cc=0;
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x; if(mask[p]){cz+=z;cy+=y;cx+=x;cc++;} }
    if(cc>0){cz/=cc;cy/=cc;cx/=cc;} }
  f32*seedv3=malloc(nn*sizeof(f32)); u8*seedm3=calloc(nn,1);
  size_t pc=0; double bestd=1e30;
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x;
    double d=sqrt((z-cz)*(z-cz)+(y-cy)*(y-cy)+(x-cx)*(x-cx)); wind[p]= mask[p]?(f32)(d/pitch):0.0f;
    if(mask[p]&&d<bestd){bestd=d;pc=p;} }
  seedm3[pc]=1; seedv3[pc]=wind[pc];   // pin the most-central material voxel to its init value
  wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)pitch; wp.iters=iters;
  // SOR: the diagonal/isotropic operator is stable up to omega~1 and converges ~2-3x faster than the
  // old under-relaxed 0.3. The off-diagonal cross terms (lagged) need heavy damping -> stay at ~0.4.
  wp.omega = (tensor6_on==1 || tensor6) ? 0.4f : 1.0f;
  { const char*e=getenv("OMEGA"); if(e){ wp.omega=(f32)atof(e); } }
  fprintf(stderr,"omega=%.2f\n",wp.omega);
  wp.cross_relax=1.0f; { const char*e=getenv("CROSS"); if(e){ wp.cross_relax=(f32)atof(e); fprintf(stderr,"CROSS (off-diag relax) -> %.3f\n",wp.cross_relax); } }
  wp.warm_start=1; wp.anchor_lambda=0.0f; wp.aniso=aniso; wp.tensor6=tensor6;   // 3D: Dirichlet pin (not anchor) fixes the constant; (void)anchor
  if(aniso_alpha>0.0 && tensor6_on==1){ wp.normal=nrm; wp.sheetness=sheet; wp.tensor_alpha=(f32)aniso_alpha; }   // on-the-fly D (no array)

  // FORCE_ITER passes: pass 0 uses the TRUE-3D globally-oriented normals (orient_normals_3d, geometry-
  // following, no umbilicus); later passes re-orient to agree with the CURRENT winding gradient grad(W)
  // -> self-consistent refinement. Each pass rebuilds forcing=div(oriented n/pitch) and re-solves.
  int npass = usenormals ? force_iter : 1;
  double gradbias=0.4; { const char*e=getenv("GRADBIAS"); if(e)gradbias=atof(e); }   // outward tie-break in gradW re-orient
  for(int pass=0; pass<npass; pass++){
    if(usenormals){
      f32*gx=calloc(nn,sizeof(f32)),*gy=calloc(nn,sizeof(f32)),*gz=calloc(nn,sizeof(f32));
      for(int z=0;z<dz;z++){
        for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x; if(!mask[p])continue;
          double nx=nrm[3*p+0],ny=nrm[3*p+1],nz=nrm[3*p+2];
          if(pass==0){ /* already globally oriented in 3D; use as-is */ }
          else { // orient along the current field gradient (central diff, masked) + outward tie-breaker
            double gxf = (x+1<dx&&mask[p+1]&&x>0&&mask[p-1]) ? 0.5*(wind[p+1]-wind[p-1]) : 0.0;
            double gyf = (y+1<dy&&mask[p+dx]&&y>0&&mask[p-dx]) ? 0.5*(wind[p+dx]-wind[p-dx]) : 0.0;
            double gzf = (z+1<dz&&mask[p+(size_t)dy*dx]&&z>0&&mask[p-(size_t)dy*dx]) ? 0.5*(wind[p+(size_t)dy*dx]-wind[p-(size_t)dy*dx]) : 0.0;
            // outward tie-breaker (scaled to |gradW| so it only decides ambiguous/noisy voxels, where a
            // raw gradW flip would create global sign-domains; strong fold gradients still win).
            double gm=sqrt(gxf*gxf+gyf*gyf+gzf*gzf), ox=x-cx,oy=y-cy,oz=z-cz,ol=sqrt(ox*ox+oy*oy+oz*oz);
            double odot = ol>1e-6 ? (nx*ox+ny*oy+nz*oz)/ol : 0.0;
            if(nx*gxf+ny*gyf+nz*gzf + gradbias*gm*odot < 0){nx=-nx;ny=-ny;nz=-nz;}
          }
          double pl = pitchf ? (double)pitchf[p] : pitch;   // 3D local pitch (VARPITCH) or uniform
          gx[p]=(f32)(nx/pl); gy[p]=(f32)(ny/pl); gz[p]=(f32)(nz/pl); } }   // full magnitude (sheetness gates aniso only)
      if(!forcing) forcing=calloc(nn,sizeof(f32));
      for(int z=1;z<dz-1;z++)for(int y=1;y<dy-1;y++)for(int x=1;x<dx-1;x++){ size_t p=(size_t)z*dy*dx+(size_t)y*dx+x; if(!mask[p]){forcing[p]=0;continue;}
        double div=0.5*((double)gx[p+1]-gx[p-1]) + 0.5*((double)gy[p+dx]-gy[p-dx]) + 0.5*((double)gz[p+(size_t)dy*dx]-gz[p-(size_t)dy*dx]);
        forcing[p]=(f32)div; }
      free(gx);free(gy);free(gz);
      fprintf(stderr,"forcing pass %d/%d (%s-oriented)\n",pass+1,npass,pass==0?"3D-propagated":"gradW");
    }
    wp.forcing=forcing;
    if(winding_field_solve(mask,dz,dy,dx,&ul,&wp,seedv3,seedm3,wind)){fprintf(stderr,"solve fail\n");return 1;}
  }
  free(seedv3); free(seedm3);

  // PHASE 3: RIDGE-ANCHORED re-solve (env WIND_ANCHOR=1). A sheet is one wrap, so the RADIAL winding
  // Wr = W - theta/2pi is ~constant along it. Snap Wr to its nearest integer AT each detected ridge
  // and pin seed = round(Wr)+theta/2pi as a Dirichlet condition, then re-solve. This forces the field
  // through the REAL sheets at integer iso-surfaces -> deforms exactly as the scroll does, no pitch
  // prior. The diagonal anisotropic tensor interpolates between ridges ALONG the sheets.
  if(wind_anchor>=2 && ridge){
    // INTRINSIC radial-count winding (built FROM the ridges, NOT the smooth field): per z-slice, bin
    // ridge voxels by azimuth, sort by radius, cluster into ring-crossings by radius gaps -> each ring
    // gets its radial RANK. seed W = rank + theta/2pi (auto-continuous across the theta=+-pi cut: rank
    // jumps +1 exactly as theta/2pi drops 1). Binning+clustering is robust to ridge porosity (a gap at
    // one angle is covered by neighbours in the bin). This removes the pitch prior entirely.
    const double TWO_PI=6.283185307179586;
    f32*seedv=malloc(nn*sizeof(f32)); u8*seedm=calloc(nn,1);
    int NB=360; double gapthr=pitch*0.45; if(gapthr<3)gapthr=3;
    { const char*e=getenv("ANCHOR_GAP"); if(e)gapthr=atof(e); }
    size_t cap=(size_t)dy*dx; rv_t*rv=malloc(cap*sizeof(rv_t));
    long ns=0; double rankmax=0;
    for(int z=0;z<dz;z++){ f32 cyz,cxz; umbilicus_center(&ul,(f32)z,&cyz,&cxz);
      size_t zoff=(size_t)z*dy*dx; int m=0;
      for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=zoff+(size_t)y*dx+x;
        if(!ridge[p]||!mask[p])continue;
        double ry=y-cyz,rx=x-cxz,r=sqrt(ry*ry+rx*rx),th=atan2(ry,rx);
        rv[m].r=(float)r; rv[m].th=(float)(th/TWO_PI); rv[m].li=y*dx+x;
        int b=(int)((th+M_PI)/TWO_PI*NB); if(b<0)b=0; if(b>=NB)b=NB-1; rv[m].bin=b; m++; }
      qsort(rv,m,sizeof(rv_t),rv_cmp);
      int i=0;
      while(i<m){ int b=rv[i].bin; int j=i; double rank=0; double lastr=-1e9;
        while(j<m && rv[j].bin==b){
          if(lastr>-1e8 && rv[j].r-lastr>gapthr) rank+=1.0;   // radius gap = next ring outward
          lastr=rv[j].r;
          double W=rank + rv[j].th;
          size_t p2=zoff+(size_t)rv[j].li; seedv[p2]=(f32)W; seedm[p2]=1; ns++;
          if(rank>rankmax)rankmax=rank; j++; }
        i=j; }
    }
    free(rv);
    fprintf(stderr,"PHASE3 radial-count: %ld ridge seeds, max rank %.0f (gap=%.1f); re-solving\n",ns,rankmax,gapthr);
    wp.forcing=NULL; wp.warm_start=1;
    winding_field_solve(mask,dz,dy,dx,&ul,&wp,seedv,seedm,wind);
    free(seedv); free(seedm);
  } else if(wind_anchor && ridge){
    const double TWO_PI=6.283185307179586;
    f32*seedv=malloc(nn*sizeof(f32)); u8*seedm=calloc(nn,1);
    // [WIND_ANCHOR=1, reference/failed] Per-z 2D connected components, snap median Wr of smooth field.
    int *lab=malloc((size_t)dy*dx*sizeof(int));
    int *stk=malloc((size_t)dy*dx*sizeof(int));
    int *mbuf2=malloc((size_t)dy*dx*sizeof(int));   // component pixel list (collected during DFS)
    float *mbuf=malloc((size_t)dy*dx*sizeof(float));
    long ns=0, ncomp=0, drop=0;
    for(int z=0;z<dz;z++){ f32 cyz,cxz; umbilicus_center(&ul,(f32)z,&cyz,&cxz);
      size_t zoff=(size_t)z*dy*dx;
      for(size_t i=0;i<(size_t)dy*dx;i++) lab[i]=0;
      int cur=0;
      for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ int li=y*dx+x; size_t p=zoff+li;
        if(!ridge[p]||!mask[p]||lab[li])continue;
        cur++; int sp=0; stk[sp++]=li; lab[li]=cur; int mc=0;
        int *comp=mbuf2; int csz=0;   // DFS collects the component's pixel list
        while(sp){ int q=stk[--sp]; comp[csz++]=q; int qy=q/dx,qx=q%dx; size_t pq=zoff+q;
          double th=atan2((double)qy-cyz,(double)qx-cxz)/TWO_PI; mbuf[mc++]=(float)((double)wind[pq]-th);
          for(int dyy=-1;dyy<=1;dyy++)for(int dxx=-1;dxx<=1;dxx++){ if(!dyy&&!dxx)continue;
            int ny2=qy+dyy,nx2=qx+dxx; if(ny2<0||nx2<0||ny2>=dy||nx2>=dx)continue; int nl=ny2*dx+nx2;
            if(lab[nl]||!ridge[zoff+nl]||!mask[zoff+nl])continue; lab[nl]=cur; stk[sp++]=nl; } }
        if(mc<4){ drop++; continue; }   // tiny speckle component -> skip (don't seed)
        for(int a=1;a<mc;a++){ float vv=mbuf[a]; int b=a-1; while(b>=0&&mbuf[b]>vv){mbuf[b+1]=mbuf[b];b--;} mbuf[b+1]=vv; }
        double Wri=floor((double)mbuf[mc/2]+0.5); ncomp++;
        for(int c=0;c<csz;c++){ int li2=comp[c]; int y2=li2/dx,x2=li2%dx;
          double th=atan2((double)y2-cyz,(double)x2-cxz)/TWO_PI; size_t p2=zoff+li2;
          seedv[p2]=(f32)(Wri+th); seedm[p2]=1; ns++; }
      }
    }
    free(lab);free(stk);free(mbuf);free(mbuf2);
    fprintf(stderr,"PHASE3 anchor: %ld comps (%ld speckle dropped), %ld seeds; re-solving\n",ncomp,drop,ns);
    wp.forcing=NULL; wp.warm_start=1;
    winding_field_solve(mask,dz,dy,dx,&ul,&wp,seedv,seedm,wind);
    free(seedv); free(seedm);
  }
  if(nrm)free(nrm); if(sheet)free(sheet); if(ridge)free(ridge);

  double wmin=1e30,wmax=-1e30; for(size_t p=0;p<nn;p++) if(mask[p]){ if(wind[p]<wmin)wmin=wind[p]; if(wind[p]>wmax)wmax=wind[p]; }
  fprintf(stderr,"POISSON winding: %.1f..%.1f (%.0f wraps) over %dx%dx%d\n",wmin,wmax,wmax-wmin,dz,dy,dx);

  char fn[700]; snprintf(fn,sizeof fn,"%s_vol.f32",base); FILE*f=fopen(fn,"wb");
  if(f){ int hdr[6]={dz,dy,dx,(int)z0,(int)y0,(int)x0}; fwrite(hdr,sizeof(int),6,f);
    float*wo=malloc(nn*sizeof(float)); for(size_t p=0;p<nn;p++) wo[p]=mask[p]?wind[p]:NAN;
    fwrite(wo,sizeof(float),nn,f); free(wo); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  umbilicus_free(&umb); free(ul.z);free(ul.y);free(ul.x); free(v);free(mask);free(wind); if(forcing)free(forcing); if(aniso)free(aniso); if(tensor6)free(tensor6); if(pitchf)free(pitchf); mca_close(arc);
  return 0;
}
