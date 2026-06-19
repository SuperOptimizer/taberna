/* topo_surgery.c — see topo_surgery.h. PH-localized tunnel filling. */
#include "postproc/topo_surgery.h"
#include "eval/topo.h"       // betti_numbers_6 (fast screen)
#include "topo/cubical.h"    // cubical_features (#2 cycle extraction)

#include <stdlib.h>
#include <string.h>

int fill_tunnels(u8 *mask, int nz, int ny, int nx,
                 int win, int overlap, int radius, int iters) {
  size_t nynx=(size_t)ny*nx;
  int stride = win-overlap; if(stride<1) stride=1;
  size_t wmax=(size_t)win*win*win;
  u8 *sub=malloc(wmax); f32 *inv=malloc(wmax*sizeof(f32)); u8 *loop=malloc(wmax);
  int total=0;
  for (int it=0; it<iters; it++){
    int filled_this=0;
    for (int z0=0; z0<nz; z0+=stride)
    for (int y0=0; y0<ny; y0+=stride)
    for (int x0=0; x0<nx; x0+=stride){
      int dz=win<nz-z0?win:nz-z0, dy=win<ny-y0?win:ny-y0, dx=win<nx-x0?win:nx-x0;
      if (dz<4||dy<4||dx<4) continue;
      size_t m=(size_t)dz*dy*dx, fg=0;
      for (int z=0;z<dz;z++) for(int y=0;y<dy;y++) for(int x=0;x<dx;x++){
        u8 v=mask[(size_t)(z0+z)*nynx+(size_t)(y0+y)*nx+(x0+x)];
        sub[((size_t)z*dy+y)*dx+x]=v; fg+=v; }
      if (fg<4 || fg==m) continue;
      if (fg*2 > m) continue;   // dense interior window: skip (generic reduction is
                                // O(cells^2) on near-solid blocks; tunnels live in
                                // sparse sheet regions anyway)
      topo_betti b=betti_numbers_6(sub,dz,dy,dx);     // fast screen
      if (b.b1<=0) continue;
      if (b.b1>40) continue;   // too tangled for per-window PH; leave to morphology
      // localize tunnels: PH on the inverted window
      for (size_t i=0;i<m;i++) inv[i]= sub[i]?0.0f:1.0f;
      int nf; pers_feat *ft=cubical_features(inv,dz,dy,dx,&nf);
      memset(loop,0,m);
      size_t Nw=m;
      for (int i=0;i<nf;i++){
        if (ft[i].dim!=1 || ft[i].death>=TOPO_INF) continue;
        for (int t=0;t<ft[i].ncells;t++){
          long gid=ft[i].cells[t]; long v=gid%(long)Nw; int axis=(gid-(long)Nw)/(long)Nw;
          long step=(axis==0?1: axis==1?dx : dx*dy), v2=v+step;
          loop[v]=1; if(v2>=0&&v2<(long)Nw) loop[v2]=1;
        }
      }
      for (int i=0;i<nf;i++) free(ft[i].cells);
      free(ft);
      // dilate-fill the loop voxels (radius) into the FULL mask -> caps the tunnel
      for (int z=0;z<dz;z++) for(int y=0;y<dy;y++) for(int x=0;x<dx;x++){
        if (!loop[((size_t)z*dy+y)*dx+x]) continue;
        for (int rz=-radius;rz<=radius;rz++) for(int ry=-radius;ry<=radius;ry++) for(int rx=-radius;rx<=radius;rx++){
          if (rz*rz+ry*ry+rx*rx>radius*radius) continue;
          int gz=z0+z+rz, gy=y0+y+ry, gx=x0+x+rx;
          if (gz<0||gz>=nz||gy<0||gy>=ny||gx<0||gx>=nx) continue;
          mask[(size_t)gz*nynx+(size_t)gy*nx+gx]=1;
        }
      }
      total++; filled_this++;
    }
    if (!filled_this) break;   // converged
  }
  free(sub);free(inv);free(loop);
  return total;
}
