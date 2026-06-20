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

  // RGB overlay on the slice data: recto=green, verso=blue (outer=faint red), PPM image
  int zc=nz/2; size_t off=(size_t)zc*nynx;
  u8*rgb=malloc(nynx*3);
  for(size_t p=0;p<nynx;p++){ size_t i=off+p;
    int g=v[i]*2; if(g>255)g=255;                 // brighten the (faint) papyrus
    u8 R=g,G=g,B=g;
    if(seed[i]==3){ R=40; G=255; B=40; }           // recto = green
    else if(seed[i]==2){ R=50; G=120; B=255; }     // verso = blue
    else if(seed[i]==1){ R=(u8)(120+g/3); G=g/2; B=g/2; }  // outer = faint red
    rgb[3*p+0]=R; rgb[3*p+1]=G; rgb[3*p+2]=B; }
  char fn[600]; snprintf(fn,sizeof fn,"%s.ppm",base);
  FILE*f=fopen(fn,"wb"); if(f){ fprintf(f,"P6\n%d %d\n255\n",d,d); fwrite(rgb,1,nynx*3,f); fclose(f);
    printf("wrote %s (recto=GREEN verso=BLUE outer=red, on slice data)\n",fn); }
  return 0;
}
