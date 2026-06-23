/* affinity_seg — run the SIGNED-AFFINITY + Mutex-Watershed wrap segmentation (the touching-sheets
 * machinery: affinity.c + partition.c) on a REAL mca region, and color a mid-z slice by segment so we
 * can SEE whether physically-touching wraps land in SEPARATE labels (the thing sheet_sep3d's radial
 * valley-counting merges). This is the decisive test of whether the classical signed repulsion
 * (across-sheet contact through high-sheetness material = "crossing into the touching next wrap")
 * actually separates touches on real CT, not just synthetic shells (run_e1).
 *
 *   affinity_seg ARC OUT lod z0 y0 x0 dz dy dx [krepel=1.0] [kattract=1.0] [zc=mid]
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/affinity.h"
#include "segmentation/partition.h"

static int air_threshold(const u8 *v, size_t n){
  long h[256]={0},tot=0; double s=0; long nz=0;
  for(size_t i=0;i<n;i++){ int x=v[i]; if(x){s+=x;nz++;} if(x>=1&&x<=254){h[x]++;tot++;} }
  if(tot<256){ int t=nz?(int)(0.5*s/nz+0.5):1; return t<1?1:t; }
  double sum=0; for(int i=1;i<=254;i++) sum+=(double)i*h[i];
  double sB=0,wB=0,best=-1; int thr=1;
  for(int t=1;t<=254;t++){ wB+=h[t]; if(!wB)continue; double wF=tot-wB; if(wF<=0)break;
    sB+=(double)t*h[t]; double mB=sB/wB,mF=(sum-sB)/wF,b=wB*wF*(mB-mF)*(mB-mF); if(b>best){best=b;thr=t;} }
  return thr;
}
int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [krepel=1] [kattract=1] [zc=mid]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  double krepel=argc>10?atof(argv[10]):1.0, kattract=argc>11?atof(argv[11]):1.0;
  int zc=argc>12?atoi(argv[12]):dz/2;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  size_t nn=(size_t)dz*dy*dx; int athr=air_threshold(v,nn);
  u8*mask=malloc(nn); long nmat=0; for(size_t p=0;p<nn;p++){ mask[p]=v[p]>=athr; nmat+=mask[p]; }
  fprintf(stderr,"region %dx%dx%d air<%d, %ld material\n",dz,dy,dx,athr,nmat);

  f32*sh=malloc(nn*sizeof(f32)),*nrm=malloc(3*nn*sizeof(f32)),*vf=malloc(nn*sizeof(f32));
  for(size_t p=0;p<nn;p++) vf[p]=v[p];
  st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=1.5f;
  st_sheet_detect(vf,dz,dy,dx,&sp,sh,nrm); free(vf);
  double smx=0; for(size_t p=0;p<nn;p++) if(sh[p]>smx)smx=sh[p];
  for(size_t p=0;p<nn;p++) sh[p]=(f32)(sh[p]/(smx>1e-6?smx:1));   // sheetness in [0,1] for affinity

  affinity_params ap={ .k_attract=(f32)kattract, .k_repel=(f32)krepel };
  sgraph g; s32*node_of=NULL; u32*voxel_of=NULL;
  if(affinity_build_voxel(mask,sh,nrm,dz,dy,dx,&ap,&g,&node_of,&voxel_of)){ fprintf(stderr,"affinity fail\n"); return 1; }
  fprintf(stderr,"signed graph: %d nodes, %d edges (krepel=%.1f)\n",g.nnodes,g.nedges,krepel);
  u32*labels=malloc((size_t)g.nnodes*sizeof(u32));
  int ncl=mws_partition(&g,labels);
  fprintf(stderr,"Mutex-Watershed: %d segments over %d material nodes\n",ncl,g.nnodes);

  // color the zc slice by segment label (random color per label); save PPM
  u8*rgb=calloc((size_t)dy*dx*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x; size_t o=((size_t)y*dx+x)*3;
    int g0=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g0;       // dim CT background
    if(mask[p] && node_of[p]>=0){ u32 L=labels[node_of[p]]; u32 h=L*2654435761u;
      rgb[o]=(u8)(80+(h&127)); rgb[o+1]=(u8)(80+((h>>8)&127)); rgb[o+2]=(u8)(80+((h>>16)&127)); } }
  char fn[700]; snprintf(fn,sizeof fn,"%s_seg.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,(size_t)dy*dx*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb); sgraph_free(&g); free(node_of);free(voxel_of);free(labels);free(sh);free(nrm);free(mask);free(v); mca_close(arc);
  return 0;
}
