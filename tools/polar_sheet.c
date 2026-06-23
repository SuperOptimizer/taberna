/* polar_sheet — polar-unroll a z-slice's raw CT AND structure-tensor sheetness around the umbilicus,
 * for testing angular optical-flow sheet-tracking (does OF follow the spiral / bridge touches, and does
 * the sheetness preprocessing help vs raw CT?). Output: OUT_ct.pgm and OUT_sh.pgm, both (NA cols = angle
 * 0..2pi, NR rows = radius 0..Rmax). Sheets are near-horizontal bands; sweeping columns = going around.
 *
 *   polar_sheet ARC OUT lod z0 y0 x0 dz dy dx [zc=mid] [NA=2048]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"
#include "annotate/umbilicus.h"
#include "segmentation/sheet_tensor.h"

static int air_threshold(const u8 *v, size_t n){
  long h[256]={0},tot=0; double allsum=0; long allnz=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x){allsum+=x;allnz++;} if(x>=1&&x<=254){h[x]++;tot++;} }
  if(tot<256){ int t=allnz?(int)(0.5*allsum/allnz+0.5):1; return t<1?1:(t>254?254:t); }
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sumB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){ wB+=h[t]; if(wB==0)continue; double wF=tot-wB; if(wF<=0)break;
    sumB+=(double)t*h[t]; double mB=sumB/wB,mF=(sum-sumB)/wF,btw=wB*wF*(mB-mF)*(mB-mF); if(btw>best){best=btw;thr=t;} }
  return thr;
}
static double bilin(const f32*im,int dy,int dx,double y,double x){
  if(y<0||x<0||y>dy-1||x>dx-1)return 0; int y0=(int)y,x0=(int)x,y1=y0<dy-1?y0+1:y0,x1=x0<dx-1?x0+1:x0;
  double fy=y-y0,fx=x-x0;
  return (1-fy)*((1-fx)*im[y0*dx+x0]+fx*im[y0*dx+x1])+fy*((1-fx)*im[y1*dx+x0]+fx*im[y1*dx+x1]);
}
static void wppgm(const char*fn,const u8*img,int w,int h){ FILE*f=fopen(fn,"wb"); if(f){fprintf(f,"P5\n%d %d\n255\n",w,h);fwrite(img,1,(size_t)w*h,f);fclose(f);fprintf(stderr,"wrote %s (%dx%d)\n",fn,w,h);} }

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [zc=mid] [NA=2048]\n",argv[0]); return 2; }
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]);
  int dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int zc=argc>10?atoi(argv[10]):dz/2; int NA=argc>11?atoi(argv[11]):2048;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  int cl=5; double cs=(double)(1<<cl); int ccz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,ccz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse fail\n");return 1;}
  size_t ccn=(size_t)ccz*ccy*ccx; u8*ccm=malloc(ccn); for(size_t i=0;i<ccn;i++)ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,ccz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;} free(ccm);
  double ls=(double)(1<<lod); double coarse_z=(z0+zc)*ls/cs; if(coarse_z<0)coarse_z=0; if(coarse_z>ccz-1)coarse_z=ccz-1;
  f32 ucy,ucx; umbilicus_center(&umb,(f32)coarse_z,&ucy,&ucx);
  double cy=ucy*cs/ls-y0, cx=ucx*cs/ls-x0;
  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  fprintf(stderr,"zc=%d center(%.0f,%.0f) air<%d\n",zc,cy,cx,athr);
  // sheetness on the whole subvol (3D tensor), then take slice zc
  f32*sh=malloc(nn*sizeof(f32)),*nrm=malloc(3*nn*sizeof(f32)),*vf=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=1.5f;
  st_sheet_detect(vf,dz,dy,dx,&sp,sh,nrm); free(vf); free(nrm);
  for(size_t p=0;p<nn;p++) if(v[p]<athr) sh[p]=0;
  double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p];
  for(size_t p=0;p<nn;p++) sh[p]=(f32)(sh[p]/(smx>1e-6?smx:1)*255.0);
  // radius extent
  double Rmx=0; const u8*vz=v+(size_t)zc*dy*dx; const f32*shz=sh+(size_t)zc*dy*dx;
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++) if(vz[y*dx+x]>=athr){ double r=hypot(y-cy,x-cx); if(r>Rmx)Rmx=r; }
  int NR=(int)Rmx; if(NR<8){fprintf(stderr,"degenerate radius\n");return 1;}
  fprintf(stderr,"polar %d angles x %d radii (Rmax=%.0f)\n",NA,NR,Rmx);
  u8*ctp=malloc((size_t)NA*NR),*shp=malloc((size_t)NA*NR);
  for(int a=0;a<NA;a++){ double th=2*M_PI*a/NA, sy=sin(th), sx=cos(th);
    for(int r=0;r<NR;r++){ double yy=cy+r*sy, xx=cx+r*sx; size_t o=(size_t)r*NA+a;
      double y=yy,x=xx; int ok=!(y<0||x<0||y>dy-1||x>dx-1);
      double ctval=0,shval=0;
      if(ok){ int iy=(int)y,ix=(int)x,iy1=iy<dy-1?iy+1:iy,ix1=ix<dx-1?ix+1:ix; double fy=y-iy,fx=x-ix;
        ctval=(1-fy)*((1-fx)*vz[iy*dx+ix]+fx*vz[iy*dx+ix1])+fy*((1-fx)*vz[iy1*dx+ix]+fx*vz[iy1*dx+ix1]);
        shval=bilin(shz,dy,dx,y,x); }
      ctp[o]=(u8)(ctval<0?0:ctval>255?255:ctval); shp[o]=(u8)(shval<0?0:shval>255?255:shval); } }
  char fn[700]; snprintf(fn,sizeof fn,"%s_ct.pgm",base); wppgm(fn,ctp,NA,NR);
  snprintf(fn,sizeof fn,"%s_sh.pgm",base); wppgm(fn,shp,NA,NR);
  free(v);free(sh);free(ctp);free(shp); umbilicus_free(&umb); mca_close(arc); return 0;
}
