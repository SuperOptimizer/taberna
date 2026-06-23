/* trace_mca — run taberna's advancing-front sheet tracer (the grower) on an .mca region and emit
 * PER-SHEET labels. This is the #2 (coherence: the front bridges porosity by snapping over gaps) +
 * #3 (grower: normal-consistency gate stops wrap-jumps) step — turns the fragmented per-slice ridge
 * into coherent, separated, traced sheets. Air-cut (fysics valley) first so it never seeds on wisps.
 *
 *   trace_mca ARC OUT lod z0 y0 x0 dz dy dx [aggr=1.5] [seed=0.4] [cont=0.2] [normcos=0.6] [minsize=400] [zc]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "segmentation/trace.h"
#include "fysics.h"

static int aircut(const u8*v,int dz,int dy,int dx,double aggr,u8*mask){
  size_t n=(size_t)dz*dy*dx;
  f32*in=malloc(n*sizeof(f32)),*sm=malloc(n*sizeof(f32)),*tmp=malloc(n*sizeof(f32));
  for(size_t i=0;i<n;i++) in[i]=v[i]/255.0f;
  fy_box_smooth(in,sm,tmp,dz,dy,dx,5);
  long hist[256]={0}; for(size_t i=0;i<n;i++){ int b=(int)(sm[i]*255.0f+0.5f); if(b<0)b=0; if(b>255)b=255; hist[b]++; }
  int dark=0,light=0,valley=0; double d=fy_valley_depth(hist,&dark,&light,&valley);
  int phys=39, cut;
  if(d>=0 && valley>phys){ if(aggr<0)aggr=0; if(aggr>2)aggr=2;
    cut = aggr<=1.0 ? (int)(phys+aggr*(valley-phys)+0.5) : (int)(valley+(aggr-1.0)*(light>valley?light-valley:0)+0.5); }
  else cut=(valley>phys)?valley:phys;
  long nm=0; for(size_t i=0;i<n;i++){ mask[i]=(sm[i]*255.0f>=cut); nm+=mask[i]; }
  fprintf(stderr,"aircut: valley=%d cut=%d (aggr=%.2f); material %ld/%zu (%.1f%%)\n",valley,cut,aggr,nm,n,100.0*nm/n);
  free(in);free(sm);free(tmp); return cut;
}

int main(int argc,char**argv){
  if(argc<10){ fprintf(stderr,"usage: %s ARC OUT lod z0 y0 x0 dz dy dx [aggr=1.5] [seed=0.4] [cont=0.2] [normcos=0.6] [minsize=400] [zc]\n",argv[0]); return 2; }
  const char*path=argv[1],*outp=argv[2]; int lod=atoi(argv[3]);
  int z0=atoi(argv[4]),y0=atoi(argv[5]),x0=atoi(argv[6]),dz=atoi(argv[7]),dy=atoi(argv[8]),dx=atoi(argv[9]);
  double aggr=argc>10?atof(argv[10]):1.5; float seed=argc>11?(float)atof(argv[11]):0.4f;
  float cont=argc>12?(float)atof(argv[12]):0.2f; float ncos=argc>13?(float)atof(argv[13]):0.6f;
  int minsize=argc>14?atoi(argv[14]):400; int zc=argc>15?atoi(argv[15]):dz/2;

  mca_handle*arc=mca_open(path); if(!arc){fprintf(stderr,"open fail\n");return 1;}
  u8*v=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!v){fprintf(stderr,"read fail\n");return 1;} mca_close(arc);
  size_t nn=(size_t)dz*dy*dx, plane=(size_t)dy*dx;
  u8*mask=malloc(nn); int cut=aircut(v,dz,dy,dx,aggr,mask);
  // pass the ORIGINAL CT (zeroing air creates artificial edges that wreck the structure-tensor normals
  // the front follows); the tracer's i_min=cut gate already keeps it out of air/wisps.
  f32*vf=malloc(nn*sizeof(f32)); for(size_t p=0;p<nn;p++) vf[p]=v[p];

  trace_params tp=trace_default_params();
  tp.st.sigma_grad=0.5f; tp.st.sigma_tensor=2.0f;   // LARGE sigma for a stable across-sheet normal (the
  { char*sg=getenv("TR_SIGT"); if(sg)tp.st.sigma_tensor=(f32)atof(sg); }  // front follows the normal)
  tp.seed_thresh=seed; tp.cont_thresh=cont; tp.i_min=(f32)cut; tp.normal_cos=ncos;
  { char*es=getenv("TR_STEP"); if(es)tp.step=(f32)atof(es); char*sn=getenv("TR_SNAP"); if(sn)tp.snap_radius=(f32)atof(sn); }
  fprintf(stderr,"step=%.2f snap=%.2f\n",tp.step,tp.snap_radius);
  if(getenv("TR_UNION")){ u8*u=malloc(nn); int n=sheet_trace(vf,dz,dy,dx,&tp,u);
    long cov=0; for(size_t p=0;p<nn;p++)cov+=u[p]; fprintf(stderr,"UNION sheet_trace: %d fronts, coverage %ld (%.1f%% of material)\n",n,cov,100.0*cov/nn); free(u); }
  s32*lab=malloc(nn*sizeof(s32));
  int ns=sheet_trace_lab(vf,dz,dy,dx,&tp,lab);
  fprintf(stderr,"traced %d fronts (seed=%.2f cont=%.2f normcos=%.2f i_min=%d)\n",ns,seed,cont,ncos,cut);

  // size each label; keep only >= minsize (drop noise fronts)
  long*sz=calloc((size_t)ns+1,sizeof(long));
  for(size_t p=0;p<nn;p++) if(lab[p]>0) sz[lab[p]]++;
  int kept=0; for(int i=1;i<=ns;i++) if(sz[i]>=minsize) kept++;
  fprintf(stderr,"kept %d/%d fronts with >=%d voxels\n",kept,ns,minsize);

  // write label volume + colored mid-z render (distinct color per kept sheet, small ones dropped)
  char fn[700]; snprintf(fn,sizeof fn,"%s_trace.i32",outp); FILE*lf=fopen(fn,"wb");
  if(lf){ int h2[6]={dz,dy,dx,z0,y0,x0}; fwrite(h2,sizeof(int),6,lf); fwrite(lab,sizeof(s32),nn,lf); fclose(lf); }
  u8*rgb=calloc(plane*3,1);
  for(int y=0;y<dy;y++)for(int x=0;x<dx;x++){ size_t p=((size_t)zc*dy+y)*dx+x,o=((size_t)y*dx+x)*3;
    int g=v[p]/4; rgb[o]=rgb[o+1]=rgb[o+2]=(u8)g;                 // CT background (dim)
    int L=lab[p]; if(L>0 && sz[L]>=minsize){ unsigned h=(unsigned)L*2654435761u;
      rgb[o]=80+(h&127); rgb[o+1]=80+((h>>8)&127); rgb[o+2]=80+((h>>16)&127); } }
  snprintf(fn,sizeof fn,"%s_trace.ppm",outp); FILE*f=fopen(fn,"wb");
  if(f){ fprintf(f,"P6\n%d %d\n255\n",dx,dy); fwrite(rgb,1,plane*3,f); fclose(f); fprintf(stderr,"wrote %s\n",fn); }
  free(rgb);free(sz);free(lab);free(vf);free(mask);free(v); return 0;
}
