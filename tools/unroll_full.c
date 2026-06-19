/* unroll_full.c — tile + stitch the whole scroll into ONE unrolled image.
 *
 * Coarse global winding (whole scroll at `clod`) gives a single coordinate system;
 * every LOD0 tile then accumulates its native-resolution intensity into the SAME
 * global unrolled image at (global_winding, z). Tiles stitch automatically BECAUSE
 * they all sample the one global field — re-solving winding per tile would let each
 * tile's relaxation drift relative to its neighbors and open seams, so we sample,
 * never re-solve. Accuracy comes from solving the global field at the finest LOD
 * that fits in RAM (LOD4 ~ 12 GB), not from per-tile refinement.
 *
 * Usage: unroll_full ARCHIVE.mca OUT.tif clod  zc0 zh  yc0 xc0 ext  tile ppw
 *   (z-slab [zc0,zc0+zh) at LOD0; cross-section box [yc0,yc0+ext)x[xc0,xc0+ext)
 *    tiled in `tile`-voxel blocks; ppw = output columns per wrap)
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
#include "unwrap/winding_field.h"

static f32 trilin(const f32*w,int nz,int ny,int nx,double z,double y,double x){
  if(z<0)z=0;if(z>nz-1)z=nz-1;if(y<0)y=0;if(y>ny-1)y=ny-1;if(x<0)x=0;if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x,z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  double dz=z-z0,dy=y-y0,dx=x-x0;
  #define G(a,b,c) w[((size_t)(a)*ny+(b))*nx+(c)]
  double c00=G(z0,y0,x0)*(1-dx)+G(z0,y0,x1)*dx,c01=G(z0,y1,x0)*(1-dx)+G(z0,y1,x1)*dx;
  double c10=G(z1,y0,x0)*(1-dx)+G(z1,y0,x1)*dx,c11=G(z1,y1,x0)*(1-dx)+G(z1,y1,x1)*dx;
  #undef G
  return (f32)((c00*(1-dy)+c01*dy)*(1-dz)+(c10*(1-dy)+c11*dy)*dz);
}

int main(int argc,char**argv){
  if(argc<11){fprintf(stderr,"usage: %s ARCHIVE OUT clod zc0 zh yc0 xc0 ext tile ppw\n",argv[0]);return 2;}
  const char*path=argv[1],*outp=argv[2]; int clod=atoi(argv[3]);
  long zc0=atol(argv[4]); int zh=atoi(argv[5]); long yc0=atol(argv[6]),xc0=atol(argv[7]); int ext=atoi(argv[8]);
  int tile=atoi(argv[9]),ppw=atoi(argv[10]);
  mca_handle *arc=mca_open(path); if(!arc){fprintf(stderr,"open failed\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  double s=(double)(1<<clod); double inv_s=1.0/s;  // s is a power of 2 -> inv_s exact
  int cz=(int)(fz/s),cy=(int)(fy/s),cx=(int)(fx/s);

  // coarse global winding (whole scroll at clod)
  fprintf(stderr,"coarse winding LOD%d (%dx%dx%d, %.1f GB vox)...\n",clod,cz,cy,cx,(double)cz*cy*cx/1e9);
  u8 *cv=mca_read(arc,clod,0,0,0,cz,cy,cx); if(!cv){fprintf(stderr,"read fail\n");return 1;}
  size_t cn=(size_t)cz*cy*cx; u8 *cm=malloc(cn);
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<cn;i++) cm[i]=cv[i]>=80;
  majority_filter(cm,cm,cz,cy,cx,1); remove_small_components(cm,cz,cy,cx,TOPO_CONN26,50);
  free(cv);

  // auto-estimate the scroll axis (no annotation) instead of assuming volume center
  umbilicus umb;
  if(umbilicus_estimate(cm,cz,cy,cx,9,&umb)){fprintf(stderr,"umbilicus estimate failed\n");return 1;}
  fprintf(stderr,"umbilicus (LOD%d): ",clod);
  for(int i=0;i<umb.n;i++) fprintf(stderr,"(%.0f,%.0f,%.0f) ",umb.z[i],umb.y[i],umb.x[i]);
  fprintf(stderr,"\n");

  // pitch is physical (~96 LOD0 voxels/wrap here); scale to this LOD's voxels
  f32 *cw=malloc(cn*sizeof(f32)); wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)(96.0/s);
  winding_field_solve(cm,cz,cy,cx,&umb,&wp,NULL,NULL,cw); free(cm);

  // winding range over the cross-section box (sample coarse field)
  float wmin=1e30f,wmax=-1e30f;
  for(int y=0;y<ext;y+=16)for(int x=0;x<ext;x+=16){
    f32 w=trilin(cw,cz,cy,cx,(zc0+zh/2)/s,(yc0+y)/s,(xc0+x)/s);
    if(w<wmin)wmin=w; if(w>wmax)wmax=w; }
  int W=(int)((wmax-wmin)*ppw)+1; if(W<1)W=1; if(W>200000)W=200000; int H=zh;
  printf("global unroll image %dx%d (winding %.1f..%.1f)\n",W,H,wmin,wmax);
  double *acc=calloc((size_t)W*H,sizeof(double)); int *cnt=calloc((size_t)W*H,sizeof(int));
  if(!acc||!cnt){fprintf(stderr,"oom for output\n");return 1;}

  // tile the cross-section box; each tile reads LOD0 + places intensity by SAMPLING
  // the one global winding field (seamless: no per-tile re-solve).
  int ntile=0;
  for(int ty=0; ty<ext; ty+=tile) for(int tx=0; tx<ext; tx+=tile){
    int th=tile<ext-ty?tile:ext-ty, tw=tile<ext-tx?tile:ext-tx;
    fprintf(stderr,"tile (%d,%d) %dx%dx%d...\n",ty,tx,zh,th,tw);
    u8 *img=mca_read(arc,0,(int)zc0,(int)(yc0+ty),(int)(xc0+tx),zh,th,tw);
    if(!img) continue;
    // Output row is p=z*W+u and the parallel-for partitions z across threads, so
    // each thread owns disjoint rows -> no two threads share a p. No atomics needed
    // (profiled as the top self-time hotspot when they were here). Keep z the
    // parallelized (outermost, non-collapsed) loop or this invariant breaks.
    #pragma omp parallel for schedule(static)
    for(int z=0;z<zh;z++)for(int y=0;y<th;y++)for(int x=0;x<tw;x++){
      size_t i=((size_t)z*th+y)*tw+x; if(img[i]<80) continue;
      f32 w=trilin(cw,cz,cy,cx,(zc0+z)*inv_s,(yc0+ty+y)*inv_s,(xc0+tx+x)*inv_s);
      int u=(int)((w-wmin)*ppw); if(u<0||u>=W) continue;
      size_t p=(size_t)z*W+u;
      acc[p]+=img[i];
      cnt[p]++;
    }
    free(img); ntile++;
  }
  u8 *flat=calloc((size_t)W*H,1); size_t fill=0;
  for(size_t i=0;i<(size_t)W*H;i++){ if(cnt[i]){flat[i]=(u8)(acc[i]/cnt[i]);fill++;} }
  tiff_save_u8(outp,flat,1,H,W);
  printf("stitched %d tiles -> %s (%dx%d), %.0f%% filled\n",ntile,outp,W,H,100.0*fill/((size_t)W*H));
  free(cw);free(acc);free(cnt);free(flat);umbilicus_free(&umb);mca_close(arc);
  return 0;
}
