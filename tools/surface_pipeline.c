/* surface_pipeline.c — the full taberna surface pipeline, everything wired in:
 *   CT -> CED (connect field) -> advancing-front tracer (watertight sheets)
 *      -> PH tunnel surgery (cubical_features localizes + caps residual tunnels)
 *      -> dust removal -> exact metric score (+ write pred TIFF).
 *
 * Usage: surface_pipeline IMAGE LABEL OUT.tif
 *        [ced_iters=10] [seed=0.4] [snap=1.5] [surg_radius=2] [surg_iters=3] [min=200]
 */
#include <stdio.h>
#include <stdlib.h>
#include "io/tiff_vol.h"
#include "segmentation/ced.h"
#include "segmentation/trace.h"
#include "postproc/morph.h"
#include "postproc/topo_surgery.h"
#include "eval/score.h"
#include "eval/topo.h"

int main(int argc,char**argv){
  if(argc<4){fprintf(stderr,"usage: %s IMAGE LABEL OUT.tif [ced][seed][snap][surg_r][surg_it][min]\n",argv[0]);return 2;}
  int ced_it=argc>4?atoi(argv[4]):10;
  float seed=argc>5?atof(argv[5]):0.4f, snap=argc>6?atof(argv[6]):1.5f;
  int sr=argc>7?atoi(argv[7]):2, sit=argc>8?atoi(argv[8]):3, minsz=argc>9?atoi(argv[9]):200;
  int nz,ny,nx,lz,ly,lx;
  u8*img=tiff_load_u8(argv[1],&nz,&ny,&nx),*lab=tiff_load_u8(argv[2],&lz,&ly,&lx);
  if(!img||!lab||lz!=nz){fprintf(stderr,"load/dim err\n");return 1;}
  size_t n=(size_t)nz*ny*nx;
  f32*volf=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) volf[i]=(f32)img[i]; free(img);

  if(ced_it>0){ fprintf(stderr,"CED (%d)...\n",ced_it); ced_params cp=ced_default_params(); cp.iters=ced_it; ced_diffuse(volf,nz,ny,nx,&cp); }

  fprintf(stderr,"tracing...\n");
  trace_params tp=trace_default_params(); tp.seed_thresh=seed; tp.snap_radius=snap; tp.min_size=minsz;
  u8*q=malloc(n); int ns=sheet_trace(volf,nz,ny,nx,&tp,q); free(volf);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,minsz);
  comp_score s0=competition_score(q,lab,nz,ny,nx,2,1,2);
  fprintf(stderr,"after trace: Score %.4f b1=%ld (%d fronts)\n",s0.score,s0.pred_b1,ns);

  // cheap morphological pre-reduction so PH surgery faces a manageable residual
  dilate_ball(q,q,nz,ny,nx,1);
  majority_filter(q,q,nz,ny,nx,3);
  comp_score s1=competition_score(q,lab,nz,ny,nx,2,1,2);
  fprintf(stderr,"after morph: Score %.4f b1=%ld\n",s1.score,s1.pred_b1);

  fprintf(stderr,"PH tunnel surgery...\n");
  int tw=fill_tunnels(q,nz,ny,nx, 16, 4, sr, sit);
  fill_holes(q,nz,ny,nx);
  remove_small_components(q,nz,ny,nx,TOPO_CONN26,minsz);
  fprintf(stderr,"surgery: %d tunnel windows\n",tw);

  tiff_save_u8(argv[3],q,nz,ny,nx);
  comp_score s=competition_score(q,lab,nz,ny,nx,2,1,2);
  printf("=== PIPELINE (CED+trace+PH-surgery) ===\n");
  printf("Score %.4f  Topo %.4f  SurfD %.4f  VOI %.4f  b0=%ld b1=%ld\n",
         s.score,s.topo_score,s.surface_dice,s.voi_score,s.pred_b0,s.pred_b1);
  free(q);free(lab);return 0;
}
