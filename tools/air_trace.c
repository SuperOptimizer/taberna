/* air_trace.c — trace the AIR (void) to label recto/verso sheet surfaces, at ANY LOD.
 * The .mca masks air as 0 (a voxel is 0 IFF air, nonzero IFF papyrus), so NO threshold
 * and NO cleaning — air = (v==0). Per z-slice (3D labels leak gaps to the exterior
 * along z): label 2D air components with cc_label; the component touching the slice
 * border is the EXTERIOR -> material bordering it is the outer scroll surface. Each
 * ENCLOSED pocket is a real delamination gap, bracketed by two reliable surfaces: the
 * RECTO (outer-wrap side, larger radius from the umbilicus) and VERSO (inner-wrap side).
 *
 * Works on a REGION at native LOD0: a quick coarse read estimates the global umbilicus
 * (scaled into region coords) so recto/verso orientation is consistent on a local patch.
 *
 * Usage: air_trace ARCHIVE.mca OUTBASE lod z0 y0 x0 d [nz=1] [minpocket=8]
 *   (region [z0,z0+nz) x [y0,y0+d) x [x0,x0+d) at `lod`; outputs OUT_air.tif of the
 *    middle slice: recto=255 verso=180 outer=130 material=70 enclosed-air=40 ext-air=20)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "io/tiff_vol.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"

int main(int argc,char**argv){
  if(argc<8){fprintf(stderr,"usage: %s ARCHIVE OUTBASE lod z0 y0 x0 d [nz=1] [minpocket=8]\n",argv[0]);return 2;}
  const char*path=argv[1],*base=argv[2]; int lod=atoi(argv[3]);
  long z0=atol(argv[4]),y0=atol(argv[5]),x0=atol(argv[6]); int d=atoi(argv[7]);
  int nz=argc>8?atoi(argv[8]):1; int minpk=argc>9?atoi(argv[9]):8;
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);

  // global umbilicus from a quick coarse (LOD5) read, scaled into this region's coords
  int cl=5; double cs=(double)(1<<cl); int cz=(int)(fz/cs),ccy=(int)(fy/cs),ccx=(int)(fx/cs);
  u8*coarse=mca_read(arc,cl,0,0,0,cz,ccy,ccx); if(!coarse){fprintf(stderr,"coarse read fail\n");return 1;}
  size_t ccn=(size_t)cz*ccy*ccx; u8*ccm=malloc(ccn);
  for(size_t i=0;i<ccn;i++) ccm[i]=coarse[i]!=0; free(coarse);
  umbilicus umb; if(umbilicus_estimate(ccm,cz,ccy,ccx,9,&umb)){fprintf(stderr,"umb fail\n");return 1;}
  free(ccm);
  // umbilicus center at THIS lod (use mid-z control point), in region-local coords
  double ls=(double)(1<<lod); f32 ucy,ucx; umbilicus_center(&umb,(f32)(cz/2),&ucy,&ucx);
  double cyf=ucy*cs/ls - y0, cxf=ucx*cs/ls - x0;   // far away for a local LOD0 patch — fine
  fprintf(stderr,"umbilicus in region coords: (y=%.0f, x=%.0f)\n",cyf,cxf);

  // read the region at `lod`
  u8*v=mca_read(arc,lod,(int)z0,(int)y0,(int)x0,nz,d,d); if(!v){fprintf(stderr,"region read fail\n");return 1;}
  size_t nynx=(size_t)d*d;
  u8*seed=calloc((size_t)nz*nynx,1); long nrec=0,nver=0;
  (void)minpk;

  // Local labeling (no enclosure: inter-sheet gaps are open channels in a patch). For
  // each papyrus voxel, step a few voxels OUTWARD and INWARD along the radius from the
  // umbilicus: if air lies outward -> this face borders the OUTER gap = RECTO (green);
  // if air lies inward -> VERSO (blue). The two alternate through the wraps.
  #pragma omp parallel for schedule(static) reduction(+:nrec,nver)
  for(int z=0;z<nz;z++){
    const u8*ms=v+(size_t)z*nynx; u8*ss=seed+(size_t)z*nynx;
    for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!ms[p])continue;
      double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1)continue;
      double uy=dy/r,ux=dx/r;                       // outward unit (away from umbilicus)
      int outair=0,inair=0;
      for(int k=1;k<=3;k++){
        int oy=(int)lround(y+uy*k),ox=(int)lround(x+ux*k);
        int iy=(int)lround(y-uy*k),ix=(int)lround(x-ux*k);
        if(oy>=0&&oy<d&&ox>=0&&ox<d&&ms[(size_t)oy*d+ox]==0) outair=1;
        if(iy>=0&&iy<d&&ix>=0&&ix<d&&ms[(size_t)iy*d+ix]==0) inair=1;
      }
      if(outair){ ss[p]=3; nrec++; }                // air outward -> recto (green)
      else if(inair){ ss[p]=2; nver++; }            // air inward  -> verso (blue)
    }
  }
  printf("region %dx%dx%d @lod%d (%ld,%ld,%ld): recto(green)=%ld verso(blue)=%ld\n",
         nz,d,d,lod,z0,y0,x0,nrec,nver);

  int zc=nz/2; size_t off=(size_t)zc*nynx;
  char fn[600];
  // (1) recto/verso overlay on the slice data
  u8*rgb=malloc(nynx*3);
  for(size_t p=0;p<nynx;p++){ size_t i=off+p;
    int g=v[i]*2; if(g>255)g=255; u8 R=g,G=g,B=g;
    if(seed[i]==3){ R=40; G=255; B=40; } else if(seed[i]==2){ R=50; G=120; B=255; }
    rgb[3*p+0]=R; rgb[3*p+1]=G; rgb[3*p+2]=B; }
  snprintf(fn,sizeof fn,"%s.ppm",base);
  FILE*f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nynx*3,f); fclose(f);
    printf("wrote %s (recto=GREEN verso=BLUE on slice data)\n",fn); }

  // (2) TRACE: connect each recto rim into a continuous sheet surface (connected
  // component) and color each traced sheet distinctly -> individual sheets pulled out
  // of the air structure.
  u8*recto=malloc(nynx);
  for(size_t p=0;p<nynx;p++) recto[p]=(seed[off+p]==3);
  u32*sl=calloc(nynx,sizeof(u32)); u32 nsheet=cc_label(recto,1,d,d,TOPO_CONN26,sl);
  // drop tiny fragments
  size_t*sa=calloc((size_t)nsheet+1,sizeof(size_t));
  for(size_t p=0;p<nynx;p++) sa[sl[p]]++;
  long nbig=0; for(u32 L=1;L<=nsheet;L++) if(sa[L]>=20) nbig++;
  // vivid distinct color per traced sheet on a dim-data background; dilate rims a bit
  for(size_t p=0;p<nynx;p++){ size_t i=off+p; int g=v[i]/4; rgb[3*p+0]=rgb[3*p+1]=rgb[3*p+2]=(u8)g; }
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; u32 L=sl[p]; if(!L||sa[L]<20)continue;
    // bright HSV-ish color from the label
    double h=(L*47%360)/60.0; int hi=(int)h; double fr=h-hi;
    int vv=255,pp=40,qq=(int)(255*(1-0.85*fr))+40*0,tt=(int)(255*(0.15+0.85*fr));
    if(qq>255)qq=255; if(tt>255)tt=255;
    u8 R,G,B; switch(hi%6){case 0:R=vv;G=tt;B=pp;break;case 1:R=qq;G=vv;B=pp;break;
      case 2:R=pp;G=vv;B=tt;break;case 3:R=pp;G=qq;B=vv;break;case 4:R=tt;G=pp;B=vv;break;
      default:R=vv;G=pp;B=qq;}
    for(int dy=-1;dy<=1;dy++)for(int dx=-1;dx<=1;dx++){int yy=y+dy,xx=x+dx;
      if(yy<0||yy>=d||xx<0||xx>=d)continue; size_t q=(size_t)yy*d+xx;
      rgb[3*q+0]=R;rgb[3*q+1]=G;rgb[3*q+2]=B; } }
  snprintf(fn,sizeof fn,"%s_sheets.ppm",base);
  f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nynx*3,f); fclose(f); }
  printf("traced %ld sheet surfaces (>=20px) -> %s\n",nbig,fn);

  // (3) WINDING from the traced sheets: relative winding = number of recto surfaces
  // crossed going outward from the umbilicus. Propagate radius-by-radius: each material
  // voxel inherits its inward neighbour's winding +1 when it ENTERS a recto rim. This
  // follows the actual (deformed, variable-pitch) sheets — no circular assumption.
  f32*wind=malloc(nynx*sizeof(f32)); for(size_t p=0;p<nynx;p++) wind[p]=-1.f;
  int rmaxp=(int)ceil(sqrt((d+fabs(cyf))*(d+fabs(cyf))+(d+fabs(cxf))*(d+fabs(cxf))))+1;
  u32*cnt=calloc(rmaxp+1,sizeof(u32)); // bucket sizes
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ if(!v[(size_t)y*d+x])continue;
    double dy=y-cyf,dx=x-cxf; int r=(int)sqrt(dy*dy+dx*dx); if(r<=rmaxp)cnt[r]++; }
  u32*boff=calloc(rmaxp+2,sizeof(u32)); for(int r=0;r<=rmaxp;r++) boff[r+1]=boff[r]+cnt[r];
  u32 nmat=boff[rmaxp+1]; size_t*order=malloc((size_t)nmat*sizeof(size_t));
  u32*cur=calloc(rmaxp+1,sizeof(u32));
  for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x; if(!v[p])continue;
    double dy=y-cyf,dx=x-cxf; int r=(int)sqrt(dy*dy+dx*dx); if(r>rmaxp)continue;
    order[boff[r]+cur[r]++]=p; }
  for(u32 k=0;k<nmat;k++){ size_t p=order[k]; int y=p/d,x=p%d;
    double dy=y-cyf,dx=x-cxf,r=sqrt(dy*dy+dx*dx); if(r<1){wind[p]=0;continue;}
    double uy=dy/r,ux=dx/r;
    int iy=(int)lround(y-uy),ix=(int)lround(x-ux); size_t ip=(size_t)iy*d+ix;
    if(iy>=0&&iy<d&&ix>=0&&ix<d&&v[ip]&&wind[ip]>=0){
      wind[p]=wind[ip];                                  // same sheet (material inward)
    } else {
      // inward neighbour is air/edge -> we just entered a new sheet across a gap.
      // search further inward (across the gap) for the last sheet's winding, +1.
      double best=-1;
      for(int st=2;st<=60;st++){ int jy=(int)lround(y-uy*st),jx=(int)lround(x-ux*st);
        if(jy<0||jy>=d||jx<0||jx>=d) break; size_t jp=(size_t)jy*d+jx;
        if(v[jp]&&wind[jp]>=0){ best=wind[jp]; break; } }
      wind[p]=(best>=0)?best+1.f:0.f;
    }
  }
  // smooth out the ray off-by-one stripes: average winding over material 4-neighbours
  // (won't cross the air gaps, so wrap separation is preserved) a few passes.
  f32*tmp=malloc(nynx*sizeof(f32));
  for(int it=0;it<25;it++){
    #pragma omp parallel for schedule(static)
    for(int y=0;y<d;y++)for(int x=0;x<d;x++){ size_t p=(size_t)y*d+x;
      if(!v[p]||wind[p]<0){tmp[p]=wind[p];continue;}
      double s=wind[p];int c=1;
      int nb[4][2]={{y,x-1},{y,x+1},{y-1,x},{y+1,x}};
      for(int k=0;k<4;k++){int yy=nb[k][0],xx=nb[k][1]; if(yy<0||yy>=d||xx<0||xx>=d)continue;
        size_t q=(size_t)yy*d+xx; if(v[q]&&wind[q]>=0){s+=wind[q];c++;} }
      tmp[p]=(f32)(s/c); }
    f32*sw=wind;wind=tmp;tmp=sw;
  }
  free(tmp);
  // visualize: each wrap = one hue cycle (frac of winding), on dim data
  for(size_t p=0;p<nynx;p++){ size_t i=off+p; int g=v[i]/4;
    u8 R=g,G=g,B=g;
    if(v[i]&&wind[p]>=0){ double f6=fmod(wind[p],6.0); int hi=(int)f6; double fr=f6-hi;
      int vv=255,pp=30,qq=(int)(255*(1-fr)),tt=(int)(255*fr);
      switch(hi){case 0:R=vv;G=tt;B=pp;break;case 1:R=qq;G=vv;B=pp;break;case 2:R=pp;G=vv;B=tt;break;
        case 3:R=pp;G=qq;B=vv;break;case 4:R=tt;G=pp;B=vv;break;default:R=vv;G=pp;B=qq;} }
    rgb[3*p+0]=R;rgb[3*p+1]=G;rgb[3*p+2]=B; }
  snprintf(fn,sizeof fn,"%s_wind.ppm",base);
  f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nynx*3,f); fclose(f); }
  float wmx=0; for(size_t p=0;p<nynx;p++) if(wind[p]>wmx)wmx=wind[p];
  printf("winding from traced sheets: %.0f wraps across patch -> %s\n",wmx,fn);
  return 0;
}
