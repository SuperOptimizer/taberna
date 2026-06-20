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
#include "io/mca.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"

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

  // read region (single slice)
  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,1,d,d); if(!v){fprintf(stderr,"region read fail\n");return 1;}
  size_t nn=(size_t)d*d;

  // (1) recto faces: material with air OUTWARD along the radius
  u8*recto=calloc(nn,1);
  #pragma omp parallel for schedule(static)
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!v[p])continue;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1)continue; double uy=dy/r,ux=dx/r;
    for(int k=1;k<=3;k++){ int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k);
      if(oy>=0&&oy<d&&ox>=0&&ox<d&&v[(size_t)oy*d+ox]==0){ recto[p]=1; break; } } }

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

  // (3) radial-adjacency graph: from each recto voxel walk outward to the next recto
  // voxel of a DIFFERENT segment -> edge[a][b]++ ("b is one wrap outward of a").
  long *edge=calloc((size_t)K*K,sizeof(long));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!recto[p])continue;
    int La=cid[sl[p]]; if(La<0)continue;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1)continue; double uy=dy/r,ux=dx/r;
    for(int k=2;k<40;k++){ int ny_=(int)lround(y+uy*k),nx_=(int)lround(x+ux*k);
      if(ny_<0||ny_>=d||nx_<0||nx_>=d)break; size_t q=(size_t)ny_*d+nx_;
      if(recto[q]){ int Lb=cid[sl[q]]; if(Lb>=0&&Lb!=La){ edge[(size_t)La*K+Lb]++; break; } } } }

  // (4) fit winding: relax wind[b]=wind[a]+1 over edges (weighted). Anchor innermost to 0.
  double *wind=calloc(K,sizeof(double));
  int inner=0; for(int i=1;i<K;i++) if(meanr[i]<meanr[inner]) inner=i;
  for(int it=0;it<2000;it++){
    double *nw=calloc(K,sizeof(double)); double *den=calloc(K,sizeof(double));
    for(int a=0;a<K;a++)for(int b=0;b<K;b++){ long w=edge[(size_t)a*K+b]; if(w<3)continue;
      nw[a]+=w*(wind[b]-1); den[a]+=w;          // a wants to be one less than b
      nw[b]+=w*(wind[a]+1); den[b]+=w; }         // b wants to be one more than a
    nw[inner]+=5*0.0; den[inner]+=5;             // anchor
    for(int i=0;i<K;i++) if(den[i]>0) wind[i]=0.5*wind[i]+0.5*(nw[i]/den[i]);
    free(nw);free(den);
  }
  // residual + range
  double sres=0; long ne=0; double wmin=1e30,wmax=-1e30;
  for(int a=0;a<K;a++)for(int b=0;b<K;b++){ long w=edge[(size_t)a*K+b]; if(w<3)continue;
    sres+=fabs(wind[b]-wind[a]-1); ne++; }
  for(int i=0;i<K;i++){ if(wind[i]<wmin)wmin=wind[i]; if(wind[i]>wmax)wmax=wind[i]; }
  printf("segments=%d edges=%ld  winding %.1f..%.1f (%.0f wraps)  edge-residual=%.3f\n",
         K,ne,wmin,wmax,wmax-wmin,ne?sres/ne:0);

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
  return 0;
}
