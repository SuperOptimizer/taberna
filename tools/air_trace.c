/* air_trace.c — trace the AIR (void) to seed sheet surfaces. The air IS the structure:
 *   - the EXTERIOR air (connected to the volume border) -> its boundary with material is
 *     the OUTER scroll surface, which gives the macro deformation away from a circle;
 *   - each ENCLOSED air pocket (a gap where wraps separate) is bracketed by two reliable
 *     sheet surfaces: the RECTO of the outer wrap and the VERSO of the inner wrap.
 * These are unambiguous (clear void on one side) -> ideal seeds for tracing/winding.
 *
 * Outputs a mid-z slice visualization OUT_air.tif (gray levels: exterior air=20,
 * enclosed air=40, material=70, outer surface=130, verso seed=180, recto seed=255) and
 * a seed-class volume OUT_seed.{raw via tiff} if requested. 3D: components are labeled
 * over the whole slab; recto/verso split is by radius about the per-z umbilicus.
 *
 * Usage: air_trace ARCHIVE.mca OUTBASE clod zc [minpocket=10]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "io/tiff_vol.h"
#include "postproc/morph.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"

int main(int argc,char**argv){
  if(argc<5){fprintf(stderr,"usage: %s ARCHIVE OUTBASE clod zc [minpocket=10]\n",argv[0]);return 2;}
  const char*path=argv[1],*base=argv[2]; int clod=atoi(argv[3]); long zc=atol(argv[4]);
  int minpk=argc>5?atoi(argv[5]):10;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  double s=(double)(1<<clod); int cz=(int)(fz/s),cy=(int)(fy/s),cx=(int)(fx/s);
  size_t cn=(size_t)cz*cy*cx, nynx=(size_t)cy*cx;
  u8*cv=mca_read(arc,clod,0,0,0,cz,cy,cx); if(!cv){fprintf(stderr,"read fail\n");return 1;}
  // RAW threshold throughout — NO connected-component / majority cleaning at LOD5: it
  // merges/removes real sheet fragments and opens the thin inter-wrap gaps (which is
  // exactly the air structure we trace). The umbilicus estimator uses outer-boundary
  // symmetry, robust to internal speckle, so it doesn't need a cleaned mask.
  u8*cm=malloc(cn);
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<cn;i++) cm[i]=cv[i]>=80;
  umbilicus umb; if(umbilicus_estimate(cm,cz,cy,cx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;}

  // air = RAW sub-threshold void. Label PER SLICE (each z-ring encloses its own gap
  // network; in 3D the gaps leak to the exterior along z).
  u8*air=malloc(cn);
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<cn;i++) air[i]=!cm[i];
  u8*seed=calloc(cn,1); long nout=0,nrec=0,nver=0,npk=0;
  #pragma omp parallel for schedule(static) reduction(+:nout,nrec,nver,npk)
  for(int z=0;z<cz;z++){
    f32 cyf,cxf; umbilicus_center(&umb,(f32)z,&cyf,&cxf);
    const u8*as=air+(size_t)z*nynx; const u8*ms=cm+(size_t)z*nynx; u8*ss=seed+(size_t)z*nynx;
    u32*lbl=calloc(nynx,sizeof(u32));
    u32 nc=cc_label(as,1,cy,cx,TOPO_CONN6,lbl);     // 2D air components in this slice
    size_t*area=calloc((size_t)nc+1,sizeof(size_t)); char*isext=calloc((size_t)nc+1,1);
    double*prs=calloc((size_t)nc+1,sizeof(double));
    for(int y=0;y<cy;y++)for(int x=0;x<cx;x++){ size_t p=(size_t)y*cx+x; u32 L=lbl[p]; if(!L)continue;
      area[L]++; if(y==0||y==cy-1||x==0||x==cx-1)isext[L]=1;
      double dy=y-cyf,dx=x-cxf; prs[L]+=sqrt(dy*dy+dx*dx); }
    for(u32 L=1;L<=nc;L++) if(!isext[L]&&area[L]>=(size_t)minpk) npk++;
    for(int y=0;y<cy;y++)for(int x=0;x<cx;x++){ size_t p=(size_t)y*cx+x; if(!ms[p])continue;
      double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx);
      int nb[4][2]={{y,x-1},{y,x+1},{y-1,x},{y+1,x}}; int cls=0;
      for(int k=0;k<4;k++){int yy=nb[k][0],xx=nb[k][1]; if(yy<0||yy>=cy||xx<0||xx>=cx)continue;
        u32 L=lbl[(size_t)yy*cx+xx]; if(!L)continue;
        if(isext[L]){ if(cls<1)cls=1; }
        else if(area[L]>=(size_t)minpk){ double pr=prs[L]/area[L]; cls=(r>pr)?3:2; break; } }
      if(cls){ ss[p]=(u8)cls; if(cls==1)nout++; else if(cls==2)nver++; else nrec++; } }
    free(lbl);free(area);free(isext);free(prs);
  }
  printf("pockets(>=%d)=%ld  outer-surf=%ld  verso=%ld  recto=%ld\n",minpk,npk,nout,nver,nrec);

  // visualize the slice at zc
  int z=(int)(zc/s); if(z<0)z=0; if(z>=cz)z=cz-1; size_t off=(size_t)z*nynx;
  u8*vis=malloc(nynx);
  for(size_t p=0;p<nynx;p++){ size_t i=off+p; u8 v;
    if(seed[i]==3)v=255; else if(seed[i]==2)v=180; else if(seed[i]==1)v=130;
    else if(cm[i])v=70; else v=40;
    vis[p]=v; }
  char fn[600]; snprintf(fn,sizeof fn,"%s_air.tif",base); tiff_save_u8(fn,vis,1,cy,cx);
  printf("wrote %s (recto=255 verso=180 outer=130 material=70 pocket-air=40 ext-air=20)\n",fn);
  return 0;
}
