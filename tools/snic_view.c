/* snic_view — run 3D SNIC supervoxels on the sheetness field of a real mca region and color a slice
 * by supervoxel, to check whether supervoxels cleanly HUG the sheets (the prerequisite for clean
 * supervoxel-level signed affinity -- the designed fix for voxel-resolution over-fragmentation).
 *
 *   snic_view ARC OUT lod z0 y0 x0 dz dy dx [d_seed=4] [compactness=20] [zc=mid]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/snic.h"

static int air_threshold(const u8*v,size_t n){ long h[256]={0},t=0; double s=0; long z=0;
  for(size_t i=0;i<n;i++){int x=v[i]; if(x){s+=x;z++;} if(x>=1&&x<=254){h[x]++;t++;}}
  if(t<256)return z?(int)(0.5*s/z+0.5):1; double sum=0; for(int i=1;i<=254;i++)sum+=(double)i*h[i];
  double sB=0,wB=0,b=-1; int thr=1; for(int k=1;k<=254;k++){wB+=h[k]; if(!wB)continue; double wF=t-wB; if(wF<=0)break;
    sB+=(double)k*h[k]; double mB=sB/wB,mF=(sum-sB)/wF,bb=wB*wF*(mB-mF)*(mB-mF); if(bb>b){b=bb;thr=k;}} return thr; }

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [d_seed=4] [compact=20] [zc=mid]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int dseed=argc>10?atoi(argv[10]):4; float compact=argc>11?atof(argv[11]):20.0f; int zc=argc>12?atoi(argv[12]):dz/2;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  f32*sh=malloc(nn*sizeof(f32)),*nrm=malloc(3*nn*sizeof(f32)),*vf=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=1.5f;
  st_sheet_detect(vf,dz,dy,dx,&sp,sh,nrm); free(vf); free(nrm);
  double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p];
  for(size_t p=0;p<nn;p++) sh[p]=(f32)(sh[p]/(smx>1e-6?smx:1)*255.0);   // SNIC wants a scalar field

  int nsp=snic_superpixel_count(dz,dy,dx,dseed);
  u32*labels=calloc(nn,sizeof(u32)); Superpixel*sps=calloc((size_t)nsp+1,sizeof(Superpixel));
  int rc=snic(sh,dz,dy,dx,dseed,compact,80.0f,160.0f,labels,sps);
  fprintf(stderr,"SNIC: d_seed=%d -> %d supervoxels (rc=%d) over %dx%dx%d\n",dseed,nsp,rc,dz,dy,dx);

  // color zc slice by supervoxel; dim CT background
  u8*rgb=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;
    if(v[p]>=athr){ u32 L=labels[p],h=L*2654435761u;
      rgb[o]=(u8)(70+(h&150)); rgb[o+1]=(u8)(70+((h>>8)&150)); rgb[o+2]=(u8)(70+((h>>16)&150)); } }
  char fn[700]; snprintf(fn,sizeof fn,"%s_snic.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb);free(labels);free(sps);free(sh);free(v); mca_close(arc); return 0;
}
