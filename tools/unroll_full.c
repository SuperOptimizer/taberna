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

/* Data-driven radial pitch = median radial spacing between consecutive sheet PEAKS
 * (intensity local maxima) along rays from the umbilicus. The wraps are mostly in
 * contact (the .mca masks only true vacuum as 0), so sheets are intensity maxima, not
 * air-bounded runs — and per-ray (not angle-averaged: averaging smears the spiral).
 * Light radial smoothing + a prominence test suppress LOD-graininess. No thresholding.
 * This IS dr_per_winding (one sheet per wrap). */
static double measure_pitch(const u8*v,int nz,int ny,int nx,const umbilicus*umb,double fallback){
  int z=nz/2; f32 cyf,cxf; umbilicus_center(umb,(f32)z,&cyf,&cxf);
  size_t nynx=(size_t)ny*nx; int maxr=(int)(0.9*(ny<nx?ny:nx)*0.5); if(maxr<12) return fallback;
  int *hist=calloc(maxr+2,sizeof(int)); long total=0;
  double *p=malloc(maxr*sizeof(double)),*sm=malloc(maxr*sizeof(double));
  for(int a=0;a<360;a++){ double ca=cos(a*M_PI/180),sa=sin(a*M_PI/180); int rmax=0;
    for(int r=0;r<maxr;r++){ int yy=(int)(cyf+r*sa),xx=(int)(cxf+r*ca);
      if(yy<0||yy>=ny||xx<0||xx>=nx)break; p[r]=v[(size_t)z*nynx+(size_t)yy*nx+xx]; rmax=r; }
    if(rmax<12) continue;
    for(int r=0;r<=rmax;r++){ double s=0;int c=0; for(int d=-2;d<=2;d++){int rr=r+d; if(rr>=0&&rr<=rmax){s+=p[rr];c++;}} sm[r]=s/c; }
    int last=-1;
    for(int r=3;r<=rmax-3;r++){
      if(sm[r]>sm[r-1]&&sm[r]>=sm[r+1]&&sm[r]>sm[r-3]+2&&sm[r]>sm[r+3]+2){   // prominent peak
        if(last>0){ int g=r-last; if(g>=2&&g<maxr){hist[g]++;total++;} } last=r; } }
  }
  free(p);free(sm);
  if(total<50){ free(hist); return fallback; }
  long acc=0; int med=(int)fallback; for(int g=2;g<maxr;g++){ acc+=hist[g]; if(acc*2>=total){med=g;break;} }
  free(hist); return med>=2?med:fallback;
}

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
  int citers=argc>11?atoi(argv[11]):-1;   // coarse winding iters (-1 = solver default)
  int warp=argc>12?atoi(argv[12]):0;       // 1 = contour-warp init (deform to scroll shape)
  mca_handle *arc=mca_open(path); if(!arc){fprintf(stderr,"open failed\n");return 1;}
  int fz,fy,fx,nl; float ql; mca_handle_dims(arc,&fz,&fy,&fx,&ql,&nl);
  double s=(double)(1<<clod); double inv_s=1.0/s;  // s is a power of 2 -> inv_s exact
  int cz=(int)(fz/s),cy=(int)(fy/s),cx=(int)(fx/s);

  // coarse global winding (whole scroll at clod)
  fprintf(stderr,"coarse winding LOD%d (%dx%dx%d, %.1f GB vox)...\n",clod,cz,cy,cx,(double)cz*cy*cx/1e9);
  u8 *cv=mca_read(arc,clod,0,0,0,cz,cy,cx); if(!cv){fprintf(stderr,"read fail\n");return 1;}
  size_t cn=(size_t)cz*cy*cx; u8 *cm=malloc(cn);
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<cn;i++) cm[i]=cv[i]!=0;   // .mca masks air as 0: material = nonzero.
  // No CC/majority cleaning at LOD5: at 32x downsampling any nonzero voxel is real papyrus.

  // auto-estimate the scroll axis (no annotation) instead of assuming volume center
  umbilicus umb;
  if(umbilicus_estimate(cm,cz,cy,cx,9,&umb)){fprintf(stderr,"umbilicus estimate failed\n");return 1;}
  fprintf(stderr,"umbilicus (LOD%d): ",clod);
  for(int i=0;i<umb.n;i++) fprintf(stderr,"(%.0f,%.0f,%.0f) ",umb.z[i],umb.y[i],umb.x[i]);
  fprintf(stderr,"\n");

  // pitch is data-driven: measured median sheet spacing at this LOD (= 1 wrap)
  double pitch=measure_pitch(cv,cz,cy,cx,&umb,96.0/s); free(cv);
  fprintf(stderr,"measured pitch %.2f vox (LOD%d) -> ~%.0f LOD0 vox/wrap\n",pitch,clod,pitch*s);
  f32 *cw=malloc(cn*sizeof(f32)); wfield_params wp=winding_default_params(); wp.dr_per_winding=(f32)pitch;
  if(citers>=0) wp.iters=citers;
  if(warp){ fprintf(stderr,"contour-warp init (deform to scroll shape)...\n");
    winding_contour_warp(cm,cz,cy,cx,&umb,pitch,cw); wp.warm_start=1;
    if(citers<0) wp.iters=5;   // the warp init is already deformed; heavy relaxation
  }                            // smooths it back toward circular -> only denoise lightly
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
      size_t i=((size_t)z*th+y)*tw+x; if(img[i]==0) continue;   // 0 = archive-masked air
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
  // No-GT quality metric: VERTICAL COHERENCE — a correct unroll maps each sheet to a
  // vertical band (constant winding across z), so adjacent rows correlate highly;
  // winding crossings (the outer-wrap artifact) lower it. Mean Pearson over rows.
  double csum=0; int cnz=0;
  for(int z=0;z+1<H;z++){ const u8*r0=flat+(size_t)z*W,*r1=flat+(size_t)(z+1)*W;
    double sx=0,sy=0,sxx=0,syy=0,sxy=0; int nn=0;
    for(int x=0;x<W;x++){ if(!r0[x]||!r1[x])continue; double a=r0[x],b=r1[x]; sx+=a;sy+=b;sxx+=a*a;syy+=b*b;sxy+=a*b;nn++; }
    if(nn<50)continue; double vx=sxx-sx*sx/nn,vy=syy-sy*sy/nn; if(vx<=1e-6||vy<=1e-6)continue;
    csum+=(sxy-sx*sy/nn)/sqrt(vx*vy); cnz++; }
  double coh=cnz?csum/cnz:0;
  printf("stitched %d tiles -> %s (%dx%d), %.0f%% filled, vert-coherence=%.4f\n",
         ntile,outp,W,H,100.0*fill/((size_t)W*H),coh);
  free(cw);free(acc);free(cnt);free(flat);umbilicus_free(&umb);mca_close(arc);
  return 0;
}
