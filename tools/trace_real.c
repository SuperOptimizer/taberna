/* trace_real.c — run detection/tracing on a REAL .mca scroll region (no ground
 * truth) and report no-reference quality: topology (b0 sheets, b1 tunnels, b2
 * cavities), surface fraction, and coherence (largest-component share). Good
 * sheets = few large coherent components, low b1/b2.
 *
 * Usage: trace_real ARCHIVE.mca [d=384] [lod=0] [method=ridge|trace] [z0 y0 x0]
 *        (omit z0/y0/x0 to auto-find a material-rich region)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io/mca.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "segmentation/trace.h"
#include "postproc/morph.h"
#include "eval/topo.h"

static void report(const char *tag, const u8 *m, int nz, int ny, int nx) {
  size_t n=(size_t)nz*ny*nx, fg=0;
  for (size_t i=0;i<n;i++) fg+=m[i]!=0;
  u32 *lab=malloc(n*sizeof(u32));
  u32 nc=cc_label(m,nz,ny,nx,TOPO_CONN26,lab);
  size_t *sz=calloc(nc+1,sizeof(size_t)); for(size_t i=0;i<n;i++) sz[lab[i]]++;
  size_t big=0; for(u32 c=1;c<=nc;c++) if(sz[c]>big) big=sz[c];
  topo_betti b=betti_numbers_6(m,nz,ny,nx);
  printf("%-10s surf=%zu (%.1f%%)  comps=%u  largest=%.0f%%  b0=%ld b1=%ld b2=%ld\n",
         tag, fg, 100.0*fg/n, nc, fg?100.0*big/fg:0.0, b.b0,b.b1,b.b2);
  free(lab);free(sz);
}

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: %s ARCHIVE.mca [d=384] [lod=0] [ridge|trace] [z0 y0 x0]\n",argv[0]);return 2;}
  const char*path=argv[1];
  int d=argc>2?atoi(argv[2]):384, lod=argc>3?atoi(argv[3]):0;
  const char*method=argc>4?argv[4]:"ridge";
  int nz,ny,nx,nlods; float q;
  if(mca_dims(path,&nz,&ny,&nx,&q,&nlods)!=0){fprintf(stderr,"open failed\n");return 1;}
  printf("archive %s: dims z=%d y=%d x=%d  quality=%.3f  nlods=%d\n",path,nz,ny,nx,q,nlods);
  int z0,y0,x0;
  if(argc>=8){ z0=atoi(argv[5]);y0=atoi(argv[6]);x0=atoi(argv[7]); }
  else { fprintf(stderr,"finding material region (frac>=0.4)...\n");
    if(mca_find_region(path,lod,d,0.4f,12345ull,&z0,&y0,&x0)!=0){fprintf(stderr,"no region\n");return 1;} }
  printf("region @ (z=%d y=%d x=%d) size %d^3 lod %d\n",z0,y0,x0,d,lod);

  u8*img=mca_load_region(path,lod,z0,y0,x0,d,d,d);
  if(!img){fprintf(stderr,"read failed\n");return 1;}
  size_t n=(size_t)d*d*d, mat=0; for(size_t i=0;i<n;i++) mat+=img[i]>=80;
  printf("loaded; material(>=80)=%zu (%.1f%%)\n",mat,100.0*mat/n);

  f32*v=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) v[i]=(f32)img[i]; free(img);
  f32*sh=malloc(n*sizeof(f32)),*nm=malloc(3*n*sizeof(f32));
  fprintf(stderr,"structure tensor...\n");
  st_params sp={0.5f,1.0f,1.0f}; st_sheet_detect(v,d,d,d,&sp,sh,nm);
  u8*q1=malloc(n);

  if(!strcmp(method,"trace")){
    trace_params tp=trace_default_params(); tp.seed_thresh=0.4f; tp.snap_radius=1.5f;
    fprintf(stderr,"tracing...\n");
    int ns=sheet_trace(v,d,d,d,&tp,q1);
    fprintf(stderr,"%d fronts\n",ns);
    remove_small_components(q1,d,d,d,TOPO_CONN26,200);
    report("trace+dust",q1,d,d,d);
  } else {
    ridge_nms(v,sh,nm,d,d,d,0.3f,80.f,1.0f,q1);
    remove_small_components(q1,d,d,d,TOPO_CONN26,100);
    report("ridge",q1,d,d,d);
    dilate_ball(q1,q1,d,d,d,1); majority_filter(q1,q1,d,d,d,2); fill_holes(q1,d,d,d);
    report("ridge+clean",q1,d,d,d);
  }
  free(v);free(sh);free(nm);free(q1);
  return 0;
}
