/* trace.c — see trace.h. Advancing-front sheet tracer. */
#include "segmentation/trace.h"
#include "segmentation/ridge.h"   // (shares the trilinear idea; own copy below)

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))

trace_params trace_default_params(void) {
  trace_params p;
  p.st = st_default_params(); p.st.sigma_grad=0.5f; p.st.sigma_tensor=1.0f;
  p.seed_thresh=0.40f; p.cont_thresh=0.20f; p.step=2.0f; p.snap_radius=1.0f; p.i_min=90.0f;
  p.normal_cos=0.6f; p.min_size=200;
  return p;
}

static f32 trilin(const f32 *v, int nz,int ny,int nx, f32 z,f32 y,f32 x){
  size_t nynx=(size_t)ny*nx;
  if(z<0)z=0; if(z>nz-1)z=nz-1; if(y<0)y=0; if(y>ny-1)y=ny-1; if(x<0)x=0; if(x>nx-1)x=nx-1;
  int z0=(int)z,y0=(int)y,x0=(int)x, z1=z0<nz-1?z0+1:z0,y1=y0<ny-1?y0+1:y0,x1=x0<nx-1?x0+1:x0;
  f32 dz=z-z0,dy=y-y0,dx=x-x0;
  f32 c00=v[IDX(z0,y0,x0)]*(1-dx)+v[IDX(z0,y0,x1)]*dx, c01=v[IDX(z0,y1,x0)]*(1-dx)+v[IDX(z0,y1,x1)]*dx;
  f32 c10=v[IDX(z1,y0,x0)]*(1-dx)+v[IDX(z1,y0,x1)]*dx, c11=v[IDX(z1,y1,x0)]*(1-dx)+v[IDX(z1,y1,x1)]*dx;
  return (c00*(1-dy)+c01*dy)*(1-dz)+(c10*(1-dy)+c11*dy)*dz;
}

/* normal at float pos (interpolate the 3 components, renormalize) */
static void normal_at(const f32 *normal, int nz,int ny,int nx, f32 z,f32 y,f32 x, f32 n[3]){
  // sample at nearest voxel (normal is an axis field; interpolation can cancel signs)
  int zi=(int)lroundf(z), yi=(int)lroundf(y), xi=(int)lroundf(x);
  if(zi<0)zi=0; if(zi>=nz)zi=nz-1; if(yi<0)yi=0; if(yi>=ny)yi=ny-1; if(xi<0)xi=0; if(xi>=nx)xi=nx-1;
  size_t nynx=(size_t)ny*nx, i=IDX(zi,yi,xi);
  n[0]=normal[3*i+2]; n[1]=normal[3*i+1]; n[2]=normal[3*i+0]; // store as (z,y,x); field is (x,y,z)
  f32 L=sqrtf(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]); if(L>1e-6f){n[0]/=L;n[1]/=L;n[2]/=L;}
}

static void set_vox(u8 *o,int nz,int ny,int nx,int z,int y,int x){ size_t nynx=(size_t)ny*nx;
  if(z>=0&&z<nz&&y>=0&&y<ny&&x>=0&&x<nx) o[IDX(z,y,x)]=1; }
static void draw_line(u8 *o,int nz,int ny,int nx, f32 z0,f32 y0,f32 x0,f32 z1,f32 y1,f32 x1){
  f32 dz=z1-z0,dy=y1-y0,dx=x1-x0; f32 L=fabsf(dz); if(fabsf(dy)>L)L=fabsf(dy); if(fabsf(dx)>L)L=fabsf(dx);
  int m=(int)L+1; for(int i=0;i<=m;i++){f32 t=m?(f32)i/m:0; set_vox(o,nz,ny,nx,(int)lroundf(z0+t*dz),(int)lroundf(y0+t*dy),(int)lroundf(x0+t*dx));}
}

/* snap a position onto the sheet ridge: search along +/- normal for the CT max */
static int snap(const f32 *ct,const f32 *sheet, int nz,int ny,int nx,
                f32 R, f32 i_min, f32 s_min, f32 n[3], f32 *z,f32 *y,f32 *x){
  f32 best=-1; f32 bz=*z,by=*y,bx=*x;
  for(f32 d=-R; d<=R+1e-3f; d+=0.5f){
    f32 zz=*z+d*n[0], yy=*y+d*n[1], xx=*x+d*n[2];
    f32 c=trilin(ct,nz,ny,nx,zz,yy,xx);
    if(c>best){best=c;bz=zz;by=yy;bx=xx;}
  }
  if(best<i_min) return 0;
  if(trilin(sheet,nz,ny,nx,bz,by,bx) < s_min) return 0;
  *z=bz;*y=by;*x=bx; return 1;
}

typedef struct { f32 z,y,x,nz,ny,nx; } frontpt;

int sheet_trace(const f32 *vol, int nz, int ny, int nx, const trace_params *pp, u8 *out){
  trace_params p = pp?*pp:trace_default_params();
  size_t nynx=(size_t)ny*nx, N=(size_t)nz*nynx;
  memset(out,0,N);
  f32 *sheet=malloc(N*sizeof(f32)), *normal=malloc(3*N*sizeof(f32));
  if(!sheet||!normal){free(sheet);free(normal);return -1;}
  if(st_sheet_detect(vol,nz,ny,nx,&p.st,sheet,normal)!=0){free(sheet);free(normal);return -1;}

  u8 *visited=calloc(N,1);
  frontpt *q=malloc(N*sizeof(frontpt)); // bounded ring not needed; cap by visited
  int nsheets=0;

  // seeds in raster order (skip visited); each unvisited strong sheet voxel starts a front
  for (int sz=0; sz<nz; sz++) for(int sy=0; sy<ny; sy++) for(int sx=0; sx<nx; sx++){
    size_t si=IDX(sz,sy,sx);
    if (visited[si] || sheet[si] < p.seed_thresh || vol[si] < p.i_min) continue;

    // start a front (visited marked at push time => queue bounded by N, no dups)
    size_t head=0, tail=0;
    f32 n0[3]; normal_at(normal,nz,ny,nx,sz,sy,sx,n0);
    f32 fz=sz,fy=sy,fx=sx;
    if(!snap(vol,sheet,nz,ny,nx,p.snap_radius,p.i_min,p.seed_thresh,n0,&fz,&fy,&fx)) { visited[si]=1; continue; }
    { int wi_z=(int)lroundf(fz),wi_y=(int)lroundf(fy),wi_x=(int)lroundf(fx);
      if(wi_z<0||wi_z>=nz||wi_y<0||wi_y>=ny||wi_x<0||wi_x>=nx){ visited[si]=1; continue; }
      size_t wi=IDX(wi_z,wi_y,wi_x); visited[wi]=1; out[wi]=1;
      q[tail++]=(frontpt){fz,fy,fx,n0[0],n0[1],n0[2]}; }

    while(head<tail){
      frontpt cp=q[head++];
      f32 n[3]={cp.nz,cp.ny,cp.nx};
      f32 a[3]={1,0,0}; if(fabsf(n[0])>0.9f){a[0]=0;a[1]=1;}
      f32 t1[3]={a[1]*n[2]-a[2]*n[1], a[2]*n[0]-a[0]*n[2], a[0]*n[1]-a[1]*n[0]};
      f32 L=sqrtf(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]); if(L<1e-6f)continue; t1[0]/=L;t1[1]/=L;t1[2]/=L;
      f32 t2[3]={n[1]*t1[2]-n[2]*t1[1], n[2]*t1[0]-n[0]*t1[2], n[0]*t1[1]-n[1]*t1[0]};
      for(int k=0;k<8;k++){
        f32 ang=k*0.7853981634f, cu=cosf(ang), cv=sinf(ang);
        f32 dz=cu*t1[0]+cv*t2[0], dy=cu*t1[1]+cv*t2[1], dx=cu*t1[2]+cv*t2[2];
        f32 qz=cp.z+p.step*dz, qy=cp.y+p.step*dy, qx=cp.x+p.step*dx;
        f32 nn[3]; normal_at(normal,nz,ny,nx,qz,qy,qx,nn);
        if(!snap(vol,sheet,nz,ny,nx,p.snap_radius,p.i_min,p.cont_thresh,nn,&qz,&qy,&qx)) continue;
        if(fabsf(n[0]*nn[0]+n[1]*nn[1]+n[2]*nn[2]) < p.normal_cos) continue; // wrap-jump guard
        int wz=(int)lroundf(qz),wy=(int)lroundf(qy),wx=(int)lroundf(qx);
        if(wz<0||wz>=nz||wy<0||wy>=ny||wx<0||wx>=nx) continue;
        draw_line(out,nz,ny,nx,cp.z,cp.y,cp.x,qz,qy,qx);   // always connect
        size_t wi=IDX(wz,wy,wx);
        if(visited[wi]) continue;
        visited[wi]=1; out[wi]=1;
        if(tail<(size_t)N) q[tail++]=(frontpt){qz,qy,qx,nn[0],nn[1],nn[2]};
      }
    }
    nsheets++;
  }
  free(sheet);free(normal);free(visited);free(q);
  return nsheets;
}

/* ---- labeled variant: identical front growth, paints a per-front s32 id instead of a union ---- */
static void set_vox_lab(s32 *o,int nz,int ny,int nx,int z,int y,int x,s32 id){ size_t nynx=(size_t)ny*nx;
  if(z>=0&&z<nz&&y>=0&&y<ny&&x>=0&&x<nx) o[IDX(z,y,x)]=id; }
static void draw_line_lab(s32 *o,int nz,int ny,int nx, f32 z0,f32 y0,f32 x0,f32 z1,f32 y1,f32 x1,s32 id){
  f32 dz=z1-z0,dy=y1-y0,dx=x1-x0; f32 L=fabsf(dz); if(fabsf(dy)>L)L=fabsf(dy); if(fabsf(dx)>L)L=fabsf(dx);
  int m=(int)L+1; for(int i=0;i<=m;i++){f32 t=m?(f32)i/m:0; set_vox_lab(o,nz,ny,nx,(int)lroundf(z0+t*dz),(int)lroundf(y0+t*dy),(int)lroundf(x0+t*dx),id);}
}
int sheet_trace_lab(const f32 *vol, int nz, int ny, int nx, const trace_params *pp, s32 *lab){
  trace_params p = pp?*pp:trace_default_params();
  size_t nynx=(size_t)ny*nx, N=(size_t)nz*nynx;
  for(size_t i=0;i<N;i++) lab[i]=0;
  f32 *sheet=malloc(N*sizeof(f32)), *normal=malloc(3*N*sizeof(f32));
  if(!sheet||!normal){free(sheet);free(normal);return -1;}
  if(st_sheet_detect(vol,nz,ny,nx,&p.st,sheet,normal)!=0){free(sheet);free(normal);return -1;}
  u8 *visited=calloc(N,1); frontpt *q=malloc(N*sizeof(frontpt)); int nsheets=0;
  int dbg = getenv("TR_DBG")!=NULL;
  long att=0,fsnap=0,fnorm=0,fbound=0,fvis=0,acc=0;
  for (int sz=0; sz<nz; sz++) for(int sy=0; sy<ny; sy++) for(int sx=0; sx<nx; sx++){
    size_t si=IDX(sz,sy,sx);
    if (visited[si] || sheet[si] < p.seed_thresh || vol[si] < p.i_min) continue;
    s32 id=nsheets+1; size_t head=0, tail=0;
    f32 n0[3]; normal_at(normal,nz,ny,nx,sz,sy,sx,n0);
    f32 fz=sz,fy=sy,fx=sx;
    if(!snap(vol,sheet,nz,ny,nx,p.snap_radius,p.i_min,p.seed_thresh,n0,&fz,&fy,&fx)) { visited[si]=1; continue; }
    { int wz=(int)lroundf(fz),wy=(int)lroundf(fy),wx=(int)lroundf(fx);
      if(wz<0||wz>=nz||wy<0||wy>=ny||wx<0||wx>=nx){ visited[si]=1; continue; }
      size_t wi=IDX(wz,wy,wx); visited[wi]=1; lab[wi]=id; q[tail++]=(frontpt){fz,fy,fx,n0[0],n0[1],n0[2]}; }
    while(head<tail){
      frontpt cp=q[head++]; f32 n[3]={cp.nz,cp.ny,cp.nx};
      f32 a[3]={1,0,0}; if(fabsf(n[0])>0.9f){a[0]=0;a[1]=1;}
      f32 t1[3]={a[1]*n[2]-a[2]*n[1], a[2]*n[0]-a[0]*n[2], a[0]*n[1]-a[1]*n[0]};
      f32 L=sqrtf(t1[0]*t1[0]+t1[1]*t1[1]+t1[2]*t1[2]); if(L<1e-6f)continue; t1[0]/=L;t1[1]/=L;t1[2]/=L;
      f32 t2[3]={n[1]*t1[2]-n[2]*t1[1], n[2]*t1[0]-n[0]*t1[2], n[0]*t1[1]-n[1]*t1[0]};
      for(int k=0;k<8;k++){
        f32 ang=k*0.7853981634f, cu=cosf(ang), cv=sinf(ang);
        f32 dz=cu*t1[0]+cv*t2[0], dy=cu*t1[1]+cv*t2[1], dx=cu*t1[2]+cv*t2[2];
        f32 qz=cp.z+p.step*dz, qy=cp.y+p.step*dy, qx=cp.x+p.step*dx;
        att++;
        f32 nn[3]; normal_at(normal,nz,ny,nx,qz,qy,qx,nn);
        if(!snap(vol,sheet,nz,ny,nx,p.snap_radius,p.i_min,p.cont_thresh,nn,&qz,&qy,&qx)) { fsnap++; continue; }
        if(fabsf(n[0]*nn[0]+n[1]*nn[1]+n[2]*nn[2]) < p.normal_cos) { fnorm++; continue; }
        int wz=(int)lroundf(qz),wy=(int)lroundf(qy),wx=(int)lroundf(qx);
        if(wz<0||wz>=nz||wy<0||wy>=ny||wx<0||wx>=nx) { fbound++; continue; }
        draw_line_lab(lab,nz,ny,nx,cp.z,cp.y,cp.x,qz,qy,qx,id);
        size_t wi=IDX(wz,wy,wx); if(visited[wi]) { fvis++; continue; }
        visited[wi]=1; lab[wi]=id; acc++;
        if(tail<(size_t)N) q[tail++]=(frontpt){qz,qy,qx,nn[0],nn[1],nn[2]};
      }
    }
    nsheets++;
  }
  if(dbg) fprintf(stderr,"[TR_DBG] growth attempts=%ld  fail: snap=%ld(%.0f%%) normcos=%ld(%.0f%%) bounds=%ld visited=%ld(%.0f%%)  accepted=%ld(%.0f%%)\n",
    att,fsnap,100.0*fsnap/(att?att:1),fnorm,100.0*fnorm/(att?att:1),fbound,fvis,100.0*fvis/(att?att:1),acc,100.0*acc/(att?att:1));
  free(sheet);free(normal);free(visited);free(q);
  return nsheets;
}
