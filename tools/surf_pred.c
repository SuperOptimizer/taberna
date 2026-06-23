/* surf_pred — generate a SMOOTH surface prediction volume (taberna's classical analog of VC3D's
 * *_surf.zarr): a [0,1] sheet-probability field that is coherent ALONG sheets and sharp ACROSS gaps,
 * so a tracer/unroll follows it without sheet-switching.
 *
 * Pipeline:  CT -> coherence-enhancing diffusion (CED: diffuse along the sheet plane, not across,
 *            connecting porous/broken sheet signal while keeping adjacent wraps apart)
 *               -> structure-tensor sheetness at SMALL sigma (Phase 0a: lone sheets peak at sigma<=0.7;
 *                  sigma~1.5 fuses two wraps) -> normalize [0,1].
 * Writes the smooth prediction (_surf.f32, _vol-style header), the RAW (no-CED) sheetness (_raw.f32)
 * for comparison, and mid-z grayscale renders of both.
 *
 *   surf_pred ARC OUT lod z0 y0 x0 dz dy dx [ced_iters=12] [sigma_tensor=0.8] [c_norm=0.05] [zc=mid] [aircut_aggr=1.5]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ced.h"
#include "segmentation/ridge.h"
#include "fysics.h"   // fy_valley_depth, fy_box_smooth — air-cut (no longer done in the export preprocess)

/* fysics air-cut: box-smooth the region, find the air<->material histogram VALLEY, cut there. The export
 * pipeline no longer zeros air, so we do it here; cutting at the valley (vs a low Otsu split) removes the
 * faint low-density WISPS that aren't dense enough to be papyrus but still trip the structure tensor.
 * `aggr`: 0 -> physics floor (conservative), 1 -> valley (removes faint material), >1 -> past valley.
 * Returns the u8 cut and writes mask[i]=1 where the SMOOTHED density is >= cut. */
static int aircut(const u8*v,int dz,int dy,int dx,double aggr,u8*mask){
  size_t n=(size_t)dz*dy*dx;
  f32*in=malloc(n*sizeof(f32)),*sm=malloc(n*sizeof(f32)),*tmp=malloc(n*sizeof(f32));
  for(size_t i=0;i<n;i++) in[i]=v[i]/255.0f;
  fy_box_smooth(in,sm,tmp,dz,dy,dx,5);              // scratch smooth (tightens the histogram modes)
  long hist[256]={0}; for(size_t i=0;i<n;i++){ int b=(int)(sm[i]*255.0f+0.5f); if(b<0)b=0; if(b>255)b=255; hist[b]++; }
  int dark=0,light=0,valley=0; double d=fy_valley_depth(hist,&dark,&light,&valley);
  int phys=39, cut;
  if(d>=0 && valley>phys){ if(aggr<0)aggr=0; if(aggr>2)aggr=2;
    cut = aggr<=1.0 ? (int)(phys+aggr*(valley-phys)+0.5) : (int)(valley+(aggr-1.0)*(light>valley?light-valley:0)+0.5); }
  else cut = (valley>phys)?valley:phys;
  long nm=0; for(size_t i=0;i<n;i++){ mask[i]=(sm[i]*255.0f>=cut); nm+=mask[i]; }
  fprintf(stderr,"aircut: dark=%d valley=%d light=%d -> cut=%d (aggr=%.2f); material %ld/%zu (%.1f%%)\n",
          dark,valley,light,cut,aggr,nm,n,100.0*nm/n);
  free(in);free(sm);free(tmp); return cut;
}

static void norm01(f32*s,size_t n){ double mx=0; for(size_t i=0;i<n;i++) if(s[i]>mx)mx=s[i]; if(mx<1e-6)mx=1;
  for(size_t i=0;i<n;i++) s[i]=(f32)(s[i]/mx); }
static void write_vol(const char*outp,const char*suf,const f32*s,const u8*mask,int dz,int dy,int dx,int z0,int y0,int x0){
  char fn[700]; snprintf(fn,sizeof fn,"%s%s",outp,suf); FILE*f=fopen(fn,"wb"); if(!f)return;
  int hd[6]={dz,dy,dx,z0,y0,x0}; fwrite(hd,sizeof(int),6,f); size_t n=(size_t)dz*dy*dx;
  f32*o=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) o[i]=(!mask||mask[i])?s[i]:0.0f;
  fwrite(o,sizeof(f32),n,f); free(o); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
static void write_render(const char*outp,const char*suf,const f32*s,const u8*mask,int dy,int dx,int zc,size_t plane){
  char fn[700]; snprintf(fn,sizeof fn,"%s%s",outp,suf); FILE*f=fopen(fn,"wb"); if(!f)return;
  u8*g=malloc(plane); for(size_t i=0;i<plane;i++){ size_t p=(size_t)zc*plane+i;
    double v=(mask&&!mask[p])?0.0:s[p]*255.0; g[i]=v<0?0:v>255?255:(u8)v; }   // air -> 0 (honor the cut)
  fprintf(f,"P5\n%d %d\n255\n",dx,dy); fwrite(g,1,plane,f); free(g); fclose(f); fprintf(stderr,"wrote %s\n",fn); }

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [ced_iters=12] [sigma_tensor=0.8] [c_norm=0.05] [zc=mid] [aircut_aggr=1.5]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  int ced_iters=argc>10?atoi(argv[10]):12; float sigt=argc>11?(float)atof(argv[11]):0.8f;
  float c_norm=argc>12?(float)atof(argv[12]):0.05f; int zc=argc>13?atoi(argv[13]):dz/2;
  double aggr=argc>14?atof(argv[14]):1.5;   // air-cut aggr: 1=valley, 1.5 default (drops faint wisps),
                                             // 2=material peak (max wisp removal, sheets stay separated)
  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;}
  mca_close(arc);
  size_t nn=(size_t)dz*dy*dx, plane=(size_t)dy*dx;
  // AIR-CUT (fysics valley) — replaces the export-preprocess step we no longer run; drops faint wisps.
  u8*mask=malloc(nn); int cut=aircut(v,dz,dy,dx,aggr,mask);
  long nmat=0; for(size_t p=0;p<nn;p++) nmat+=mask[p];
  fprintf(stderr,"CED iters=%d sigma_tensor=%.2f c_norm=%.3f\n",ced_iters,sigt,c_norm);

  // CT zeroed in air so the structure tensor never fires on wisp edges
  f32*vf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vf[p]=mask[p]?v[p]:0.0f;
  // RAW sheetness (no CED) for comparison
  f32*raw=malloc(nn*sizeof(f32));
  st_params sp=st_default_params(); sp.sigma_grad=1.0f; sp.sigma_tensor=sigt;
  st_sheet_detect(vf,dz,dy,dx,&sp,raw,NULL); norm01(raw,nn);
  write_vol(outp,"_raw.f32",raw,mask,dz,dy,dx,z0,y0,x0);
  write_render(outp,"_raw.ppm",raw,mask,dy,dx,zc,plane);

  // CED-enhanced: diffuse the CT along sheets, then detect sheetness on the connected field
  fprintf(stderr,"coherence-enhancing diffusion...\n");
  ced_params cp=ced_default_params(); cp.iters=ced_iters; cp.c_norm=c_norm; cp.st.sigma_tensor=sigt;
  ced_diffuse(vf,dz,dy,dx,&cp);
  f32*surf=malloc(nn*sizeof(f32)),*nrm=malloc(3*nn*sizeof(f32));
  st_sheet_detect(vf,dz,dy,dx,&sp,surf,nrm); norm01(surf,nn);
  write_vol(outp,"_surf.f32",surf,mask,dz,dy,dx,z0,y0,x0);
  write_render(outp,"_surf.ppm",surf,mask,dy,dx,zc,plane);

  // THINNING: collapse the thick CED'd sheetness to its ~1-voxel medial centerline (ridge NMS along the
  // normal) = taberna's customThinning analog -> the clean, trace-friendly smooth-surface prediction.
  u8*ridge=malloc(nn); size_t nr=ridge_nms(vf,surf,nrm,dz,dy,dx,0.3f,(f32)cut,1.0f,ridge);
  fprintf(stderr,"ridge centerline: %zu voxels (%.2f%% of material)\n",nr,100.0*nr/nmat);
  f32*rf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) rf[p]=ridge[p]?1.0f:0.0f;
  write_vol(outp,"_ridge.f32",rf,mask,dz,dy,dx,z0,y0,x0);
  write_render(outp,"_ridge.ppm",rf,mask,dy,dx,zc,plane);
  free(nrm);free(ridge);free(rf);

  // quick quality proxies: coverage of confident sheet (>0.5) and mean over material
  long rc=0,sc=0; double rm=0,sm=0; for(size_t p=0;p<nn;p++){ if(!mask[p])continue;
    if(raw[p]>0.5)rc++; if(surf[p]>0.5)sc++; rm+=raw[p]; sm+=surf[p]; }
  fprintf(stderr,"confident-sheet(>0.5) frac: raw %.3f  CED %.3f ; mean sheetness: raw %.3f CED %.3f\n",
          (double)rc/nmat,(double)sc/nmat,rm/nmat,sm/nmat);
  free(v);free(vf);free(raw);free(surf);free(mask); return 0;
}
