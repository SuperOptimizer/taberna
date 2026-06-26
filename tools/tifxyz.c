/* tifxyz — read/render VC3D "tifxyz" QuadSurfaces (the interop format for surface predictions).
 *
 * A VC3D segment is a directory of single-page float32 TIFFs: x.tif / y.tif / z.tif are (u,v) grids
 * whose pixels hold the VOLUME coordinate (in ZYX volume space) of each flattened surface point;
 * value -1 = invalid/unfilled. (Optional mask.tif; we treat x<0 as invalid.) Sampling a CT archive at
 * each (z,y,x) renders the FLATTENED papyrus — i.e. plug a VC3D surface into our pipeline.
 *
 * COORDINATE NOTE: the paris4_segments meta says _volume_shape=[18946,8174,8174] = the Paris4 L2 grid,
 * so segment coords map 1:1 to a --level 2 .mca (scale=1). Use scale=2 against an L1 .mca, 4 against L0.
 *
 *   tifxyz info   SEGDIR
 *   tifxyz render SEGDIR ARC OUT.pgm [lod=0] [scale=1.0] [uvtile=48]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tiff.h"
#include "io/mca.h"

typedef struct { int H, W; float *X, *Y, *Z; } TifXYZ;   // (u,v) grids of volume coords; X<0 = invalid

static int load_plane(const char *segdir, const char *name, tiff_volume *v){
  char p[1024]; snprintf(p, sizeof p, "%s/%s.tif", segdir, name);
  if (tiff_read_volume(p, v) != TIFF_OK) { fprintf(stderr, "read %s failed\n", p); return -1; }
  if (v->type != TIFF_F32 || v->channels != 1 || v->depth != 1) {
    fprintf(stderr, "%s: expected single-page f32 1-channel (got type=%d ch=%d depth=%u)\n",
            p, v->type, v->channels, v->depth); tiff_volume_free(v); return -1; }
  return 0;
}
static int load_tifxyz(const char *segdir, TifXYZ *s){
  tiff_volume vx, vy, vz;
  if (load_plane(segdir,"x",&vx)) return -1;
  if (load_plane(segdir,"y",&vy)) { tiff_volume_free(&vx); return -1; }
  if (load_plane(segdir,"z",&vz)) { tiff_volume_free(&vx); tiff_volume_free(&vy); return -1; }
  if (vx.width!=vy.width||vx.width!=vz.width||vx.height!=vy.height||vx.height!=vz.height){
    fprintf(stderr,"x/y/z dims differ\n"); return -1; }
  s->W=(int)vx.width; s->H=(int)vx.height;
  s->X=(float*)vx.data; s->Y=(float*)vy.data; s->Z=(float*)vz.data;   // take ownership
  return 0;
}
static void free_tifxyz(TifXYZ *s){ free(s->X); free(s->Y); free(s->Z); }

int main(int argc, char **argv){
  if (argc < 3){ fprintf(stderr,"usage: %s info SEGDIR | %s render SEGDIR ARC OUT.pgm [lod=0] [scale=1.0] [uvtile=48]\n",argv[0],argv[0]); return 2; }
  const char *mode=argv[1], *segdir=argv[2];
  TifXYZ s; if (load_tifxyz(segdir,&s)) return 1;
  size_t N=(size_t)s.H*s.W; long nv=0;
  double zmn=1e30,zmx=-1e30,ymn=1e30,ymx=-1e30,xmn=1e30,xmx=-1e30;
  for(size_t i=0;i<N;i++){ if(s.X[i]<0||!isfinite(s.X[i]))continue; nv++;
    double X=s.X[i],Y=s.Y[i],Z=s.Z[i];
    if(Z<zmn)zmn=Z; if(Z>zmx)zmx=Z; if(Y<ymn)ymn=Y; if(Y>ymx)ymx=Y; if(X<xmn)xmn=X; if(X>xmx)xmx=X; }
  fprintf(stderr,"tifxyz %s: (u,v)=%dx%d, valid %ld/%zu (%.1f%%)\n  bbox z[%.0f,%.0f] y[%.0f,%.0f] x[%.0f,%.0f]\n",
          segdir,s.H,s.W,nv,N,100.0*nv/N,zmn,zmx,ymn,ymx,xmn,xmx);
  if (!strcmp(mode,"info")){ free_tifxyz(&s); return 0; }

  if (strcmp(mode,"render")||argc<5){ fprintf(stderr,"render needs: SEGDIR ARC OUT.pgm\n"); free_tifxyz(&s); return 2; }
  const char *arcpath=argv[3], *outp=argv[4];
  int lod=argc>5?atoi(argv[5]):0; double scale=argc>6?atof(argv[6]):1.0; int uv=argc>7?atoi(argv[7]):48;
  double lods=1.0/(double)(1<<lod);   // seg->lod0 via `scale`, then origin-subtract, then >>lod to the read grid
  mca_handle *arc=mca_open(arcpath); if(!arc){fprintf(stderr,"open arc fail\n"); free_tifxyz(&s); return 1; }
  int az,ay,ax,nl; float ql; mca_handle_dims(arc,&az,&ay,&ax,&ql,&nl);
  int lz=az>>lod, ly=ay>>lod, lx=ax>>lod;
  // VC3D segment coords live in the FULL source-volume (world) frame; our archive is trimmed to a
  // ROI, so archive(0,0,0) = world(roi.origin). Read that origin from the archive's metadata carveout
  // and subtract it (after `scale` maps seg -> this archive's LOD0 world). Without this the surface
  // samples air offset by the trim. TIFXYZ_NO_ROI=1 disables (coords already archive-local).
  int roz=0,roy=0,rox=0;
  if(getenv("TIFXYZ_NO_ROI")) fprintf(stderr,"ROI origin: disabled (TIFXYZ_NO_ROI)\n");
  else if(mca_roi_origin(arc,&roz,&roy,&rox)==0) fprintf(stderr,"ROI origin (subtracted): z%d y%d x%d\n",roz,roy,rox);
  else fprintf(stderr,"ROI origin: none in archive metadata -> assuming 0 (untrimmed)\n");
  fprintf(stderr,"render vs %s lod%d (%dx%dx%d), coord scale=%.3f, uvtile=%d\n",arcpath,lod,lz,ly,lx,scale,uv);

  u8 *out=calloc(N,1); long sampled=0,oob=0;
  // Tile the (u,v) grid; each tile's valid coords span a bounded 3D bbox -> one region read, sample in-memory.
  for(int v0=0; v0<s.H; v0+=uv) for(int u0=0; u0<s.W; u0+=uv){
    int v1=v0+uv<s.H?v0+uv:s.H, u1=u0+uv<s.W?u0+uv:s.W;
    double bz0=1e30,bz1=-1e30,by0=1e30,by1=-1e30,bx0=1e30,bx1=-1e30; int any=0;
    for(int v=v0;v<v1;v++)for(int u=u0;u<u1;u++){ size_t i=(size_t)v*s.W+u; if(s.X[i]<0||!isfinite(s.X[i]))continue;
      double Z=(s.Z[i]*scale-roz)*lods,Y=(s.Y[i]*scale-roy)*lods,X=(s.X[i]*scale-rox)*lods; any=1;
      if(Z<bz0)bz0=Z; if(Z>bz1)bz1=Z; if(Y<by0)by0=Y; if(Y>by1)by1=Y; if(X<bx0)bx0=X; if(X>bx1)bx1=X; }
    if(!any)continue;
    int z0=(int)floor(bz0), y0=(int)floor(by0), x0=(int)floor(bx0);
    int dz=(int)floor(bz1)-z0+2, dy=(int)floor(by1)-y0+2, dx=(int)floor(bx1)-x0+2;
    if(z0<0)z0=0; if(y0<0)y0=0; if(x0<0)x0=0;
    if(z0+dz>lz)dz=lz-z0; if(y0+dy>ly)dy=ly-y0; if(x0+dx>lx)dx=lx-x0;
    if(dz<=0||dy<=0||dx<=0)continue;
    if((double)dz*dy*dx > 600e6){ fprintf(stderr,"  tile (%d,%d) bbox %dx%dx%d too big -> skip (raise this cap or shrink uvtile)\n",v0,u0,dz,dy,dx); continue; }
    u8 *reg=mca_read(arc,lod,z0,y0,x0,dz,dy,dx); if(!reg)continue;
    for(int v=v0;v<v1;v++)for(int u=u0;u<u1;u++){ size_t i=(size_t)v*s.W+u; if(s.X[i]<0||!isfinite(s.X[i]))continue;
      int Z=(int)lround((s.Z[i]*scale-roz)*lods)-z0, Y=(int)lround((s.Y[i]*scale-roy)*lods)-y0, X=(int)lround((s.X[i]*scale-rox)*lods)-x0;
      if(Z<0||Y<0||X<0||Z>=dz||Y>=dy||X>=dx){ oob++; continue; }
      out[i]=reg[((size_t)Z*dy+Y)*dx+X]; sampled++; }
    free(reg);
  }
  mca_close(arc);
  FILE *f=fopen(outp,"wb"); if(f){ fprintf(f,"P5\n%d %d\n255\n",s.W,s.H); fwrite(out,1,N,f); fclose(f);
    fprintf(stderr,"wrote %s (%dx%d); sampled %ld, oob %ld\n",outp,s.W,s.H,sampled,oob); }
  free(out); free_tifxyz(&s); return 0;
}
