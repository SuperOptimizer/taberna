/* wind_tv — Phase 3a of docs/touching-sheets-plan.md: the GLOBAL ordered-label entry the research flags
 * as "lower-friction" (continuous-max-flow / convex relaxation, NOT a from-scratch Boykov-Kolmogorov
 * maxflow). The two lightweight families are both exhausted and converge here: the downstream local
 * relabelings (3b) and the upstream sheet-normal field forcing (2a) BOTH lose to plain floor(W), because
 * the dominant defect — floor(W) flipping ALONG a sheet where the winding W crosses an integer over the
 * sheet's length — is GLOBAL (both halves of the sheet are locally self-consistent, so no local filter
 * has a minority to flip).
 *
 * FORMULATION. Regularize the winding field W into u by the weighted-TV-L2 (ROF) energy
 *      min_u   lambda/2 * ||u - W||^2  +  integral g(x) |grad u|        (over MATERIAL voxels only)
 * with the edge weight g = sheetness. Where g is HIGH (a sheet interior) the TV term forces u FLAT, so
 * the whole sheet collapses to ONE level => round/floor(u) is constant along the sheet and the flip is
 * gone — a GLOBAL effect a local mode-filter structurally cannot produce. Where g is LOW (gaps / sheet
 * edges) u is free to jump by ~1 to the next wrap, so adjacent wraps stay separated; the data term
 * lambda*||u-W|| anchors u to the winding so wraps don't collapse together. Solved by Chambolle-Pock
 * primal-dual (matrix-free, tileable, convex — verifiable, unlike a hand-rolled graph maxflow).
 *
 * SCOPE of this build: fixes defect 1 (along-sheet floor flip), the pervasive one. Defect 2 (LEAKED
 * touches, where W itself failed to count the turn) needs the additional monotone-OUTWARD ordering
 * constraint (Ishikawa) — W can't separate what it didn't resolve, and neither can a data term tied to
 * it. That ordering term is the next increment on top of this solver.
 *
 * PHASE 3b ORDERING (order_mu>0, default 0=OFF): TESTED-NEGATIVE on L1. A soft one-sided MIN-SLOPE
 * penalty du/dn >= delta=1/pitch along the smoothed-grad-W outward normal (n = normalize(blur(grad W)),
 * which stays defined THROUGH a leak by borrowing direction from its edges). Goal: force the +1 increment
 * W missed across a fused touch. RESULT: it cannot. Broad threshold -> fires everywhere TV softened the
 * slope, globally re-steepening u back toward floor(W) (flip RISES 1.58->1.79%, RMS|u-W| falls); severe
 * (leak-only) threshold -> near no-op (flip +0.02, the TV-merged wrap NEVER recovered, labels stay 30).
 * MECHANISM: splitting a GAPLESS fused touch is a discrete cut-PLACEMENT decision ("a +1 must cross here,
 * at the lowest-sheetness location") that a soft continuous penalty structurally cannot make — it smears
 * the increment or skips it, never pins a sharp boundary where no local gap exists. This empirically
 * confirms the research: defect 2 needs the HARD ordered-label min-cut (Ishikawa/LOGISMOS) or the spiral
 * fit, not a soft term. The ordering scaffolding (smoothed outward-normal field, masked trilinear/blur) is
 * kept default-OFF as reusable infra for that hard solver (it already builds the per-voxel column direction).
 *
 * MEASURED (canonical touch region L2 z3936 y718 x1580, 32x1024x1024, prior = wind_poisson field):
 * along-sheet floor-flip fraction (2D structure-tensor sheetness + along-tangent label change) drops
 * monotonically with the data weight while all 31 wrap labels survive (no collapse) — floor(W) 3.35% ->
 * TV lambda 0.50 2.68% (-20%) -> 0.25 2.41% (-28%) -> 0.12 2.14% (-36%); distinct-labels-per-radial-ray
 * 19.63 -> 19.22 (lambda 0.25), a ~2% change much of which is spurious-flip removal (a flip splits one
 * sheet into 2 labels, inflating the count). lambda~0.3 default = effective but conservative on merging.
 * This is the FIRST relabeling that BEATS plain floor(W) (the 3b/2a lightweight attempts all tied/lost).
 *
 *   wind_tv ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [lambda=0.3] [niter=250] [zc=mid] [gmin=0.03] [pitch=0:auto] [order_mu=0]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"

static int air_threshold(const u8*v,size_t n){ long h[256]={0},t=0; double s=0; long z=0;
  for(size_t i=0;i<n;i++){int x=v[i]; if(x){s+=x;z++;} if(x>=1&&x<=254){h[x]++;t++;}}
  if(t<256)return z?(int)(0.5*s/z+0.5):1; double sum=0; for(int i=1;i<=254;i++)sum+=(double)i*h[i];
  double sB=0,wB=0,b=-1; int thr=1; for(int k=1;k<=254;k++){wB+=h[k]; if(!wB)continue; double wF=t-wB; if(wF<=0)break;
    sB+=(double)k*h[k]; double mB=sB/wB,mF=(sum-sB)/wF,bb=wB*wF*(mB-mF)*(mB-mF); if(bb>b){b=bb;thr=k;}} return thr; }

static double cw_trilin(const float*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0,acc=0,wt=0;
  int zs[2]={z0,z1},ys[2]={y0,y1},xs[2]={x0,x1}; double zw[2]={1-dz,dz},yw[2]={1-dy,dy},xw[2]={1-dx,dx};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){
    double vv=w[((size_t)zs[a]*ny+ys[b])*nx+xs[c]]; if(!isfinite(vv))continue;
    double gg=zw[a]*yw[b]*xw[c]; acc+=gg*vv; wt+=gg; }
  return wt>1e-9? acc/wt : NAN;
}

// masked trilinear sample of a field over the working dz x dy x dx grid; material-weighted, NAN if no material
static double tri_m(const f32*f,const u8*mask,int dz,int dy,int dx,double z,double y,double x){
  if(z<0||y<0||x<0||z>dz-1||y>dy-1||x>dx-1) return NAN;
  int z0=(int)z,y0=(int)y,x0=(int)x,z1=z0<dz-1?z0+1:z0,y1=y0<dy-1?y0+1:y0,x1=x0<dx-1?x0+1:x0;
  double fz=z-z0,fy=y-y0,fx=x-x0,acc=0,wt=0;
  int zs[2]={z0,z1},ys[2]={y0,y1},xs[2]={x0,x1}; double zw[2]={1-fz,fz},yw[2]={1-fy,fy},xw[2]={1-fx,fx};
  for(int a=0;a<2;a++)for(int b=0;b<2;b++)for(int c=0;c<2;c++){ size_t q=((size_t)zs[a]*dy+ys[b])*dx+xs[c];
    if(!mask[q])continue; double g=zw[a]*yw[b]*xw[c]; acc+=g*f[q]; wt+=g; }
  return wt>1e-6?acc/wt:NAN;
}
// separable masked box blur (radius r) along one axis via running window over material voxels only
static void blur_axis(f32*a,const u8*mask,int dz,int dy,int dx,int r,int axis){
  size_t plane=(size_t)dy*dx; int n=axis==0?dz:axis==1?dy:dx; size_t st=axis==0?plane:axis==1?(size_t)dx:1;
  int outer1=axis==0?dy:dz, outer2=axis==0?dx:axis==1?dx:dy;
  #pragma omp parallel
  { double*sv=malloc((n+1)*sizeof(double)); double*sc=malloc((n+1)*sizeof(double)); f32*line=malloc(n*sizeof(f32));
    #pragma omp for collapse(2) schedule(static)
    for(int o1=0;o1<outer1;o1++)for(int o2=0;o2<outer2;o2++){
      size_t base; if(axis==0) base=(size_t)o1*dx+o2; else if(axis==1) base=(size_t)o1*plane+o2; else base=(size_t)o1*plane+(size_t)o2*dx;
      sv[0]=sc[0]=0;
      for(int i=0;i<n;i++){ size_t q=base+(size_t)i*st; int m=mask[q]; sv[i+1]=sv[i]+(m?a[q]:0.0); sc[i+1]=sc[i]+(m?1.0:0.0); }
      for(int i=0;i<n;i++){ int lo=i-r<0?0:i-r, hi=i+r+1>n?n:i+r+1; double c=sc[hi]-sc[lo];
        line[i]=(f32)(c>0?(sv[hi]-sv[lo])/c:a[base+(size_t)i*st]); }
      for(int i=0;i<n;i++) a[base+(size_t)i*st]=line[i]; }
    free(sv);free(sc);free(line); }
}

int main(int argc,char**argv){
  if(argc<12){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [lambda=0.3] [niter=250] [zc=mid] [gmin=0.03] [pitch=0:auto] [order_mu=0]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  const char*pvf=argv[10]; int plod=atoi(argv[11]);
  double lambda=argc>12?atof(argv[12]):0.3; int niter=argc>13?atoi(argv[13]):250;
  int zc=argc>14?atoi(argv[14]):dz/2; double gmin=argc>15?atof(argv[15]):0.03;
  double pitch=argc>16?atof(argv[16]):0.0; double order_mu=argc>17?atof(argv[17]):0.0;  // Phase 3b ordering

  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr?1:0; nmat+=mask[p]; }
  fprintf(stderr,"air<%d, %ld/%zu material\n",athr,nmat,nn);

  // sheetness (small sigma; Phase 0a: lone sheets peak <=0.7) -> edge weight g = gmin + sheetness
  f32*vf=malloc(nn*sizeof(f32)),*sh=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params spp=st_default_params(); spp.sigma_grad=1.0f; spp.sigma_tensor=0.8f;
  st_sheet_detect(vf,dz,dy,dx,&spp,sh,NULL); free(vf);
  { double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p]; if(smx<1e-6)smx=1;
    for(size_t p=0;p<nn;p++) sh[p]=(f32)(gmin+sh[p]/smx); }   // g in [gmin, gmin+1]

  // prior winding W per material voxel (NaN elsewhere); warm-start u = W
  FILE*pf=fopen(pvf,"rb"); if(!pf){fprintf(stderr,"priorvol open fail\n");return 1;}
  int hd[6]; if(fread(hd,sizeof(int),6,pf)!=6){fprintf(stderr,"hdr fail\n");return 1;}
  int cnz=hd[0],cny=hd[1],cnx=hd[2]; long cz0=hd[3],cy0=hd[4],cx0=hd[5];
  size_t cn=(size_t)cnz*cny*cnx; float*cw=malloc(cn*sizeof(float));
  if(fread(cw,sizeof(float),cn,pf)!=cn){fprintf(stderr,"data fail\n");return 1;} fclose(pf);
  double scl=ldexp(1.0,lod-plod);
  f32*W=malloc(nn*sizeof(f32)),*u=malloc(nn*sizeof(f32));
  long nlab=0;
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
    W[p]=NAN; if(mask[p]){ double wv=cw_trilin(cw,cnz,cny,cnx,(z0+z)*scl-cz0,(y0+y)*scl-cy0,(x0+x)*scl-cx0);
      if(isfinite(wv)){ W[p]=(f32)wv; nlab++; } else mask[p]=0; } u[p]=mask[p]?W[p]:0; }
  free(cw);
  fprintf(stderr,"%ld winding-anchored material voxels; CP weighted-TV lambda=%.3f niter=%d\n",nlab,lambda,niter);
  size_t plane=(size_t)dy*dx;

  // ---- Phase 3b ordering field: outward across-sheet direction = normalize(smooth(grad W)). Raw grad W
  // is flat AT a leak, so smoothing over ~pitch/4 borrows the outward direction from the leak's EDGES ->
  // defined THROUGH the touch. The CP loop then imposes a one-sided MIN-SLOPE  du/dn >= delta=1/pitch
  // along it: a flat leak (slope~0 < delta) is forced to ramp ~1 wrap (splitting a fused touch W missed);
  // where slope already >= delta the penalty is slack. order_mu=0 -> exact Phase 3a (no ordering). ----
  f32 *onx=NULL,*ony=NULL,*onz=NULL,*otmp=NULL; double delta=0, ostep=0;
  if(order_mu>0){
    onx=calloc(nn,sizeof(f32)); ony=calloc(nn,sizeof(f32)); onz=calloc(nn,sizeof(f32)); otmp=malloc(nn*sizeof(f32));
    double gsum=0; long gn=0;
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x; if(!mask[p])continue;
      double gx=0,gy=0,gz=0;
      if(x>0&&x<dx-1&&mask[p-1]&&mask[p+1])           gx=0.5*(W[p+1]-W[p-1]);
      if(y>0&&y<dy-1&&mask[p-dx]&&mask[p+dx])         gy=0.5*(W[p+dx]-W[p-dx]);
      if(z>0&&z<dz-1&&mask[p-plane]&&mask[p+plane])   gz=0.5*(W[p+plane]-W[p-plane]);
      onx[p]=(f32)gx; ony[p]=(f32)gy; onz[p]=(f32)gz;
      double m=sqrt(gx*gx+gy*gy+gz*gz); if(m>1e-4){ gsum+=m; gn++; } }
    if(pitch<=0){ pitch = gn? (double)gn/gsum : 40.0; fprintf(stderr,"[order] pitch auto = %.1f (1/mean|gradW|)\n",pitch); }
    int br=(int)(pitch/4); if(br<1)br=1;
    for(int ax=0;ax<3;ax++){ blur_axis(onx,mask,dz,dy,dx,br,ax); blur_axis(ony,mask,dz,dy,dx,br,ax); blur_axis(onz,mask,dz,dy,dx,br,ax); }
    for(size_t p=0;p<nn;p++){ if(!mask[p]){onx[p]=ony[p]=onz[p]=0;continue;}
      double m=sqrt((double)onx[p]*onx[p]+(double)ony[p]*ony[p]+(double)onz[p]*onz[p]);
      if(m>1e-6){ onx[p]/=m; ony[p]/=m; onz[p]/=m; } else { onx[p]=ony[p]=onz[p]=0; } }
    delta=1.0/pitch; ostep=pitch*0.5;
    fprintf(stderr,"[order] ON: mu=%.3f delta=%.4f/vox step=%.1f blur_r=%d\n",order_mu,delta,ostep,br);
  }

  // ---- Chambolle-Pock: min_u  lambda/2 |u-W|^2 + ||g grad u||_{2,1}  over material ----
  // forward-difference gradient (px,py,pz per voxel; edge inactive => stays 0), backward-diff divergence.
  f32*px=calloc(nn,sizeof(f32)),*py=calloc(nn,sizeof(f32)),*pz=calloc(nn,sizeof(f32)),*ub=malloc(nn*sizeof(f32));
  memcpy(ub,u,nn*sizeof(f32));
  const double L2=12.0; double sigma=1.0/sqrt(L2), tau=1.0/sqrt(L2);   // sigma*tau*L2 = 1
  for(int it=0;it<niter;it++){
    // dual ascent on p = proj_{|p|<=g}( p + sigma*grad(ub) )
    #pragma omp parallel for schedule(static)
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
      if(!mask[p]) continue;
      double gx=0,gy=0,gz=0;
      if(x+1<dx && mask[p+1])     gx=ub[p+1]-ub[p];
      if(y+1<dy && mask[p+dx])    gy=ub[p+dx]-ub[p];
      if(z+1<dz && mask[p+plane]) gz=ub[p+plane]-ub[p];
      double nx2=px[p]+sigma*gx, ny2=py[p]+sigma*gy, nz2=pz[p]+sigma*gz;
      double nrm=sqrt(nx2*nx2+ny2*ny2+nz2*nz2), g=sh[p];
      double sc = nrm>g ? g/nrm : 1.0;
      px[p]=(f32)(nx2*sc); py[p]=(f32)(ny2*sc); pz[p]=(f32)(nz2*sc); }
    // primal descent: u_new = (u + tau*(div p + lambda*W)) / (1 + tau*lambda); ub = 2 u_new - u
    double tl=tau*lambda;
    #pragma omp parallel for schedule(static)
    for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
      if(!mask[p]) continue;
      double div=px[p]+py[p]+pz[p];                         // -(backward neighbor) terms below
      if(x>0 && mask[p-1])     div-=px[p-1];
      if(y>0 && mask[p-dx])    div-=py[p-dx];
      if(z>0 && mask[p-plane]) div-=pz[p-plane];
      double un=(u[p]+tau*(div+lambda*W[p]))/(1.0+tl);
      if(order_mu>0){ u[p]=(f32)un; }                 // ub refreshed after the ordering pass (theta=0)
      else { ub[p]=(f32)(2.0*un-u[p]); u[p]=(f32)un; }
    }
    // ---- Phase 3b ordering pass: one-sided min-slope GD along the smoothed outward normal. Race-free
    // GATHER into otmp (reads u at p +/- step*n via masked trilinear). r>0 means the outward slope fell
    // short of delta -> push u[p] toward satisfying it (down vs the too-flat outward side, up vs inward).
    if(order_mu>0){
      double want=delta*ostep;
      #pragma omp parallel for schedule(static)
      for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
        if(!mask[p]){ otmp[p]=u[p]; continue; }
        double nxp=onx[p],nyp=ony[p],nzp=onz[p];
        if(nxp==0&&nyp==0&&nzp==0){ otmp[p]=u[p]; continue; }
        double F=0;
        double uo=tri_m(u,mask,dz,dy,dx, z+ostep*nzp, y+ostep*nyp, x+ostep*nxp);
        double ui=tri_m(u,mask,dz,dy,dx, z-ostep*nzp, y-ostep*nyp, x-ostep*nxp);
        double sev=want*0.35;   // LEAK-SPECIFIC: fire only on a SEVERE slope deficit (true leak/reversal),
                                // not the mild slope softening TV applies everywhere (which re-steepens globally)
        if(isfinite(uo)){ double slope=uo-u[p]; if(slope<sev) F-=(sev-slope); }   // outward step nearly flat -> leak
        if(isfinite(ui)){ double slope=u[p]-ui; if(slope<sev) F+=(sev-slope); }
        otmp[p]=(f32)(u[p]+order_mu*F); }
      memcpy(u,otmp,nn*sizeof(f32)); memcpy(ub,u,nn*sizeof(f32));   // theta=0: no over-relaxation with ordering
    }
  }
  free(px);free(py);free(pz);free(ub);
  if(order_mu>0){ free(onx);free(ony);free(onz);free(otmp); }

  // residual + label stats
  double res=0; long n2=0; for(size_t p=0;p<nn;p++) if(mask[p]){ double d=u[p]-W[p]; res+=d*d; n2++; }
  fprintf(stderr,"CP done: RMS|u-W|=%.4f wraps\n",sqrt(res/(n2?n2:1)));

  s32*lab=malloc(nn*sizeof(s32)); long wmin=1<<30,wmax=-(1<<30);
  for(size_t p=0;p<nn;p++){ lab[p]=-1; if(mask[p]){ int Lq=(int)floor(u[p]); lab[p]=Lq; if(Lq<wmin)wmin=Lq; if(Lq>wmax)wmax=Lq; } }
  fprintf(stderr,"floor(u): wrap bands %ld..%ld (%ld wraps)\n",wmin,wmax,wmax-wmin+1);

  // outputs: regularized field _tv.f32 (vol header), label _lab.i32, mid-z render _label.ppm
  char fn[700];
  snprintf(fn,sizeof fn,"%s_tv.f32",outp); FILE*tf=fopen(fn,"wb");
  if(tf){ int h2[6]={dz,dy,dx,z0,y0,x0}; fwrite(h2,sizeof(int),6,tf);
    float*wo=malloc(nn*sizeof(float)); for(size_t p=0;p<nn;p++) wo[p]=mask[p]?u[p]:NAN;
    fwrite(wo,sizeof(float),nn,tf); free(wo); fclose(tf); fprintf(stderr,"wrote %s\n",fn); }
  snprintf(fn,sizeof fn,"%s_lab.i32",outp); FILE*lf=fopen(fn,"wb");
  if(lf){ int h2[6]={dz,dy,dx,z0,y0,x0}; fwrite(h2,sizeof(int),6,lf); fwrite(lab,sizeof(s32),nn,lf); fclose(lf);
    fprintf(stderr,"wrote %s (%dx%dx%d i32 wrap labels)\n",fn,dz,dy,dx); }
  static const u8 pal[6][3]={{230,60,60},{240,180,40},{60,200,80},{50,180,220},{90,90,230},{220,90,210}};
  u8*rgb=calloc(plane*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
    if(lab[p]>=0){ int q=((lab[p]%6)+6)%6; rgb[o]=pal[q][0]; rgb[o+1]=pal[q][1]; rgb[o+2]=pal[q][2]; } }
  snprintf(fn,sizeof fn,"%s_label.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,plane*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb);
  free(sh);free(W);free(u);free(lab);free(mask);free(v); mca_close(arc); return 0;
}
