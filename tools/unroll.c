/* unroll.c — the full taberna scroll-unrolling pipeline, end to end on real data:
 *   .mca region -> sheet surface (structure-tensor ridge) -> umbilicus ->
 *   winding field (Jacobi) -> Archimedean spiral fit (r=a+b*theta) ->
 *   deformation field + fold check -> FLATTEN the sheet to a 2D (winding, z)
 *   image (the unrolled, readable papyrus) -> no-reference quality metrics.
 *
 * Usage: unroll ARCHIVE.mca OUT_flat.tif [d=512] [z0 y0 x0] [lod=0]
 *        (umbilicus = scroll axis through the volume center; override-free)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "io/mca.h"
#include "io/tiff_vol.h"
#include "segmentation/sheet_tensor.h"
#include "segmentation/ridge.h"
#include "postproc/morph.h"
#include "eval/topo.h"
#include "annotate/umbilicus.h"
#include "unwrap/winding_field.h"
#include "unwrap/spiral_fit.h"
#include "unwrap/wmetrics.h"
#include "unwrap/deform.h"

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr,"usage: %s ARCHIVE.mca OUT_flat.tif [d=512] [z0 y0 x0] [lod=0]\n",argv[0]); return 2; }
  const char *path=argv[1], *outp=argv[2];
  int d = argc>3?atoi(argv[3]):512;
  int fz,fy,fx,nlods; float qual;
  if (mca_dims(path,&fz,&fy,&fx,&qual,&nlods)!=0){fprintf(stderr,"open failed\n");return 1;}
  int lod = argc>7?atoi(argv[7]):0;
  int z0,y0,x0;
  if (argc>=7){ z0=atoi(argv[4]);y0=atoi(argv[5]);x0=atoi(argv[6]); }
  else if (mca_find_region(path,lod,d,0.5f,7ull,&z0,&y0,&x0)!=0){fprintf(stderr,"no region\n");return 1;}
  printf("archive %dx%dx%d; region @ (%d,%d,%d) size %d^3\n",fz,fy,fx,z0,y0,x0,d);

  u8 *img = mca_load_region(path,lod,z0,y0,x0,d,d,d);
  if(!img){fprintf(stderr,"read failed\n");return 1;}
  size_t n=(size_t)d*d*d;

  // 1) surface detection
  f32 *v=malloc(n*sizeof(f32)); for(size_t i=0;i<n;i++) v[i]=(f32)img[i];
  f32 *sh=malloc(n*sizeof(f32)),*nm=malloc(3*n*sizeof(f32));
  fprintf(stderr,"detect...\n"); st_params sp={0.5f,1.0f,1.0f}; st_sheet_detect(v,d,d,d,&sp,sh,nm);
  u8 *mask=malloc(n);
  ridge_nms(v,sh,nm,d,d,d,0.3f,80.f,1.0f,mask);
  remove_small_components(mask,d,d,d,TOPO_CONN26,200);
  size_t surf=0; for(size_t i=0;i<n;i++) surf+=mask[i];
  printf("surface voxels: %zu (%.1f%%)\n",surf,100.0*surf/n);
  free(sh);free(nm);

  // 2) umbilicus = straight scroll axis (volume center, region-relative)
  umbilicus umb; umb.n=2; umb.z=malloc(2*sizeof(f32)); umb.y=malloc(2*sizeof(f32)); umb.x=malloc(2*sizeof(f32));
  f32 cy=(f32)(fy*0.5 - y0), cx=(f32)(fx*0.5 - x0);
  umb.z[0]=0; umb.z[1]=(f32)(d-1); umb.y[0]=umb.y[1]=cy; umb.x[0]=umb.x[1]=cx;
  printf("umbilicus axis at region (y=%.0f, x=%.0f)\n",cy,cx);

  // 3) winding field
  fprintf(stderr,"winding field...\n");
  f32 *wind=malloc(n*sizeof(f32));
  wfield_params wp=winding_default_params();
  if (winding_field_solve(mask,d,d,d,&umb,&wp,NULL,NULL,wind)!=0){fprintf(stderr,"winding failed\n");return 1;}
  float wmin=1e30f,wmax=-1e30f; for(size_t i=0;i<n;i++) if(mask[i]){ if(wind[i]<wmin)wmin=wind[i]; if(wind[i]>wmax)wmax=wind[i]; }
  printf("winding range: %.2f .. %.2f  (~%.1f wraps)\n",wmin,wmax,(wmax-wmin));

  // 4) Archimedean spiral fit + quality metrics
  spiral_params spf=spiral_fit_from_field(mask,wind,d,d,d,&umb,4);
  printf("spiral fit: r = %.2f + %.3f*theta   RMS=%.2f vox  (n=%d)\n",spf.a,spf.b,spf.rms,spf.nsamples);
  double mvf=winding_mvf(mask,wind,d,d,d,&umb,64,d/2);
  printf("winding MVF (monotonicity violation): %.4f\n",mvf);
  f32 *disp=malloc(3*n*sizeof(f32));
  if (deform_build(mask,wind,d,d,d,&umb,spf,disp)==0){
    double fold=jacobian_fold_fraction(disp,d,d,d);
    printf("deformation fold fraction (det J<=0): %.4f\n",fold);
  }
  free(disp);

  // 5) FLATTEN: map each surface voxel to (winding, z) -> unrolled papyrus image
  int W = (int)(wmax-wmin)+1; if(W<1)W=1; if(W>20000)W=20000;
  int H = d;
  double *acc=calloc((size_t)W*H,sizeof(double)); int *cnt=calloc((size_t)W*H,sizeof(int));
  for (int z=0;z<d;z++) for(int y=0;y<d;y++) for(int x=0;x<d;x++){
    size_t i=((size_t)z*d+y)*d+x; if(!mask[i]) continue;
    int u=(int)(wind[i]-wmin); if(u<0)u=0; if(u>=W)u=W-1;
    acc[(size_t)z*W+u]+=img[i]; cnt[(size_t)z*W+u]++;
  }
  u8 *flat=calloc((size_t)W*H,1);
  for (size_t i=0;i<(size_t)W*H;i++) if(cnt[i]) flat[i]=(u8)(acc[i]/cnt[i]);
  tiff_save_u8(outp, flat, 1, H, W);   // single-slice TIFF (z=1, H rows, W cols)
  printf("wrote unrolled papyrus: %s  (%d x %d)\n",outp,W,H);

  free(img);free(v);free(mask);free(wind);free(acc);free(cnt);free(flat);
  free(umb.z);free(umb.y);free(umb.x);
  return 0;
}
