/* wind_label — Phase 3 of docs/touching-sheets-plan.md: turn a winding field into a clean integer
 * wrap-index INSTANCE-LABEL volume. Phase 1 established that floor(winding) is already concentric-per-
 * wrap (the field separates wraps; the supervoxel/MWS detour only fragmented). This makes that the
 * product: label[v] = floor(W[v]), DESPECKLED with an in-plane mode filter to remove the band-boundary
 * flicker where W crosses an integer. Writes a label volume (i32, same _vol-style header) + a colored
 * mid-z render (ordered palette so concentric wraps read as continuous cycling rings).
 *
 * NOT done here (the harder Phase 3, see IMPLEMENTATION LOG): forcing touches the field LEAKED to split.
 * The ordered-label idea with straight RADIAL rays as columns (LOGISMOS-style) FAILS on this data: the
 * scroll is DEFORMED, so a straight ray cuts ACROSS the deformation and crosses wraps out of winding
 * order (W sampled along a ray is non-monotone even with the exact umbilicus). The columns must follow
 * the winding GRADIENT (streamlines of grad W, perpendicular to sheets), or use a global continuous-
 * max-flow ordered-label solve on the voxel field. That is the next build; this tool ships the robust
 * instance map the field already supports.
 *
 *   wind_label ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [zc=mid] [modr=1]
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

int main(int argc,char**argv){
  if(argc<12){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx priorvol priorlod [zc=mid] [modr=1]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  const char*pvf=argv[10]; int plod=atoi(argv[11]);
  int zc=argc>12?atoi(argv[12]):dz/2; int modr=argc>13?atoi(argv[13]):1;
  double sthr=argc>14?atof(argv[14]):0.15;   // sheetness threshold for connected-component "patches"
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);

  // sheetness up front (small sigma; Phase 0a: lone sheets peak <=0.7) -> normalized [0,1]; used to GATE
  // the label smoothing (propagate a label only WITHIN a sheet, never across a gap) and for the CC diag.
  f32*vf=malloc(nn*sizeof(f32)),*sh=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params spp=st_default_params(); spp.sigma_grad=1.0f; spp.sigma_tensor=0.8f;
  st_sheet_detect(vf,dz,dy,dx,&spp,sh,NULL); free(vf);
  { double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p]; if(smx<1e-6)smx=1;
    for(size_t p=0;p<nn;p++) sh[p]=(f32)(sh[p]/smx); }

  FILE*pf=fopen(pvf,"rb"); if(!pf){fprintf(stderr,"priorvol open fail\n");return 1;}
  int hd[6]; if(fread(hd,sizeof(int),6,pf)!=6){fprintf(stderr,"hdr fail\n");return 1;}
  int cnz=hd[0],cny=hd[1],cnx=hd[2]; long cz0=hd[3],cy0=hd[4],cx0=hd[5];
  size_t cn=(size_t)cnz*cny*cnx; float*cw=malloc(cn*sizeof(float));
  if(fread(cw,sizeof(float),cn,pf)!=cn){fprintf(stderr,"data fail\n");return 1;} fclose(pf);
  double scl=ldexp(1.0,lod-plod);

  // raw label = floor(W) per material voxel (-1 = air / no winding); also keep raw W on slice zc
  s32*lab=malloc(nn*sizeof(s32)); float*Wsl=malloc((size_t)dy*dx*sizeof(float));
  long wmin=1<<30,wmax=-(1<<30),nlab=0;
  for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
    lab[p]=-1; double wv=NAN;
    if(v[p]>=athr){ wv=cw_trilin(cw,cnz,cny,cnx,(z0+z)*scl-cz0,(y0+y)*scl-cy0,(x0+x)*scl-cx0);
      if(isfinite(wv)){ int L=(int)floor(wv); lab[p]=L; nlab++; if(L<wmin)wmin=L; if(L>wmax)wmax=L; } }
    if(z==zc) Wsl[(size_t)y*dx+x]=(float)wv; }
  free(cw);
  fprintf(stderr,"floor(W): %ld labelled voxels, wrap bands %ld..%ld (%ld wraps)\n",nlab,wmin,wmax,wmax-wmin+1);

  // SHEETNESS-GATED in-plane mode filter: each sheet voxel takes the MODE label of its high-sheetness
  // (within-sheet) neighbors -> a label propagates ALONG a sheet (fixing the along-sheet floor(W) flip)
  // but does NOT cross a gap (low-sheetness neighbors are excluded). At a touch (no gap) majority wins,
  // so the wrap boundary shifts rather than collapsing (unlike pure CC). modr iterations, radius 1.
  long changed=0;
  if(modr>0){ s32*tmp=malloc(nn*sizeof(s32));
    for(int it=0;it<modr;it++){ memcpy(tmp,lab,nn*sizeof(s32));
      #pragma omp parallel for schedule(static) reduction(+:changed)
      for(int z=0;z<dz;z++)for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)z*dy+y)*dx+x;
        if(tmp[p]<0||sh[p]<sthr)continue; int cnt[9]={0}; int val[9]; int nv=0;   // local label histogram
        for(int ddy=-1;ddy<=1;ddy++)for(int ddx=-1;ddx<=1;ddx++){ int yy=y+ddy,xx=x+ddx;
          if(yy<0||yy>=dy||xx<0||xx>=dx)continue; size_t q=((size_t)z*dy+yy)*dx+xx;
          if(sh[q]<sthr)continue;                                  // GATE: only propagate within a sheet
          s32 L=tmp[q]; if(L<0)continue;
          int f=-1; for(int u=0;u<nv;u++) if(val[u]==L){f=u;break;} if(f<0){val[nv]=L;cnt[nv]=1;nv++;} else cnt[f]++; }
        int bi=-1,bc=0; for(int u=0;u<nv;u++) if(cnt[u]>bc){bc=cnt[u];bi=u;}
        if(bi>=0&&val[bi]!=tmp[p]){ lab[p]=val[bi]; changed++; } } }
    free(tmp); fprintf(stderr,"sheetness-gated mode-filter (x%d, sthr=%.2f): %ld relabeled\n",modr,sthr,changed); }

  // write label volume (i32) with the _vol-style int header
  char fn[700]; snprintf(fn,sizeof fn,"%s_lab.i32",outp); FILE*lf=fopen(fn,"wb");
  if(lf){ int h2[6]={dz,dy,dx,z0,y0,x0}; fwrite(h2,sizeof(int),6,lf); fwrite(lab,sizeof(s32),nn,lf); fclose(lf);
    fprintf(stderr,"wrote %s (%dx%dx%d i32 wrap labels)\n",fn,dz,dy,dx); }

  // colored mid-z render, ordered cycling palette
  static const u8 pal[6][3]={{230,60,60},{240,180,40},{60,200,80},{50,180,220},{90,90,230},{220,90,210}};
  u8*rgb=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
    if(lab[p]>=0){ int q=((lab[p]%6)+6)%6; rgb[o]=pal[q][0]; rgb[o+1]=pal[q][1]; rgb[o+2]=pal[q][2]; } }
  snprintf(fn,sizeof fn,"%s_label.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb);

  // ---- Phase 3b ATTEMPT (TESTED-NEGATIVE on densely-touching regions): connected-component sheet
  // labeling anchored to winding (Thaumato "patch + winding number"). A sheet = a connected component of
  // high-sheetness voxels; label = round(mean winding). HOPE: label constant per sheet + touching wraps
  // (connectivity-separable) split into 2 patches. RESULT: in a densely-touching region the high-
  // sheetness voxels form ONE connected mass spanning ALL wraps (touches bridge everything) -> one giant
  // patch with a meaningless mean-W; raising sthr (0.15->0.6) only erodes sheets into noise fragments
  // without breaking the bridges (793->1395 tiny patches around the same mega-blob). So pure connectivity
  // CANNOT separate wraps here. The robust labeler stays floor(W) above; the genuine touch fix needs the
  // GLOBAL ordered-label optimization (continuous-max-flow / graph-cut: integer L minimizing |L-W| +
  // smoothness, boundaries pinned to low-sheetness gaps) -- the heavy Phase 3b build. _cc.ppm kept as a
  // diagnostic of the bridging. Mid-z slice (per-slice 2D CC).
  u8*sm=calloc((size_t)dy*dx,1);                                     // sheet mask on slice zc
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x;
    if(v[p]>=athr && sh[p]>=sthr) sm[(size_t)y*dx+x]=1; }
  free(sh);
  s32*cc=malloc((size_t)dy*dx*sizeof(s32)); for(size_t i=0;i<(size_t)dy*dx;i++) cc[i]=-1;
  int *stk=malloc((size_t)dy*dx*sizeof(int)); int ncc=0; long merged2=0;
  // BFS flood fill, 8-connectivity; per-component accumulate mean W -> wrap label
  s32*cclab=malloc((size_t)dy*dx*sizeof(s32)); int ccn=0; int cccap=dy*dx;
  for(int y0i=0;y0i<dy;y0i++)for(int x0i=0;x0i<dx;x0i++){ size_t s0=(size_t)y0i*dx+x0i;
    if(!sm[s0]||cc[s0]>=0)continue; int top=0; stk[top++]=(int)s0; cc[s0]=ncc; double wsum=0; long wc=0;
    while(top){ int q=stk[--top]; int qy=q/dx,qx=q%dx; double wv=Wsl[q]; if(isfinite(wv)){wsum+=wv;wc++;}
      for(int ddy=-1;ddy<=1;ddy++)for(int ddx=-1;ddx<=1;ddx++){ if(!ddy&&!ddx)continue; int ny=qy+ddy,nx=qx+ddx;
        if(ny<0||ny>=dy||nx<0||nx>=dx)continue; size_t nq=(size_t)ny*dx+nx;
        if(sm[nq]&&cc[nq]<0){ cc[nq]=ncc; stk[top++]=(int)nq; } } }
    (void)cccap; cclab[ncc]= wc? (s32)lround(wsum/wc) : -1; ncc++; }
  fprintf(stderr,"Phase 3b: %d sheet patches (sthr=%.2f) on slice z=%d\n",ncc,sthr,zc);
  free(stk);
  // render: color each patch by its winding-anchored wrap label (ordered palette)
  u8*rgb2=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3,s=(size_t)y*dx+x;
    int g=v[p]/4; rgb2[o]=rgb2[o+1]=rgb2[o+2]=(u8)g;
    if(cc[s]>=0&&cclab[cc[s]]>=0){ int L=cclab[cc[s]],q=((L%6)+6)%6; rgb2[o]=pal[q][0]; rgb2[o+1]=pal[q][1]; rgb2[o+2]=pal[q][2]; } }
  snprintf(fn,sizeof fn,"%s_cc.ppm",outp); FILE*f2=fopen(fn,"wb");
  if(f2){ fprintf(f2,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb2,1,(size_t)dy*dx*3,f2); fclose(f2); fprintf(stderr,"wrote %s\n",fn); }
  (void)merged2; free(rgb2);free(cc);free(cclab);free(sm);free(Wsl);
  free(lab);free(v); mca_close(arc); return 0;
}
