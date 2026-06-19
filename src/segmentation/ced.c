/* ced.c — see ced.h. Explicit anisotropic (sheet) diffusion. */
#include "segmentation/ced.h"

#include <stdlib.h>
#include <string.h>

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))

ced_params ced_default_params(void) {
  ced_params p;
  p.st = st_default_params();
  p.st.sigma_grad = 0.5f; p.st.sigma_tensor = 1.0f;
  p.c_norm = 0.05f; p.dt = 0.12f; p.iters = 12;
  return p;
}

/* clamped index (Neumann/reflect at borders) */
static inline int cl(int i, int n) { return i<0?0:(i>=n?n-1:i); }

int ced_diffuse(f32 *vol, int nz, int ny, int nx, const ced_params *p) {
  ced_params q = p ? *p : ced_default_params();
  size_t nynx=(size_t)ny*nx, N=(size_t)nz*nynx;

  // normal field (interleaved x,y,z) + sheetness, computed once from the input
  f32 *normal=malloc(3*N*sizeof(f32)), *sheet=malloc(N*sizeof(f32));
  if(!normal||!sheet){ free(normal);free(sheet); return 1; }
  if(st_sheet_detect(vol,nz,ny,nx,&q.st,sheet,normal)!=0){ free(normal);free(sheet); return 1; }

  f32 *u=vol, *flux=malloc(3*N*sizeof(f32));
  if(!flux){ free(normal);free(sheet); return 1; }

  for(int it=0; it<q.iters; it++){
    // flux = D grad(u),  D = c_perp*I - (c_perp-c_norm) n n^T,  c_perp tied to sheetness
    #pragma omp parallel for schedule(static)
    for(int z=0;z<nz;z++) for(int y=0;y<ny;y++) for(int x=0;x<nx;x++){
      size_t i=IDX(z,y,x);
      f32 ux=0.5f*(u[IDX(z,y,cl(x+1,nx))]-u[IDX(z,y,cl(x-1,nx))]);
      f32 uy=0.5f*(u[IDX(z,cl(y+1,ny),x)]-u[IDX(z,cl(y-1,ny),x)]);
      f32 uz=0.5f*(u[IDX(cl(z+1,nz),y,x)]-u[IDX(cl(z-1,nz),y,x)]);
      f32 nx_=normal[3*i+0], ny_=normal[3*i+1], nz_=normal[3*i+2];
      f32 cper = q.c_norm + (1.0f-q.c_norm)*sheet[i];   // strong diffusion where sheet is confident
      f32 beta = cper - q.c_norm;
      f32 ndg = nx_*ux + ny_*uy + nz_*uz;               // n . grad u
      // D grad u = cper*grad u - beta*(n.gradu)*n
      flux[3*i+0]=cper*ux - beta*ndg*nx_;
      flux[3*i+1]=cper*uy - beta*ndg*ny_;
      flux[3*i+2]=cper*uz - beta*ndg*nz_;
    }
    // u += dt * div(flux)
    #pragma omp parallel for schedule(static)
    for(int z=0;z<nz;z++) for(int y=0;y<ny;y++) for(int x=0;x<nx;x++){
      size_t i=IDX(z,y,x);
      f32 dfx=0.5f*(flux[3*IDX(z,y,cl(x+1,nx))+0]-flux[3*IDX(z,y,cl(x-1,nx))+0]);
      f32 dfy=0.5f*(flux[3*IDX(z,cl(y+1,ny),x)+1]-flux[3*IDX(z,cl(y-1,ny),x)+1]);
      f32 dfz=0.5f*(flux[3*IDX(cl(z+1,nz),y,x)+2]-flux[3*IDX(cl(z-1,nz),y,x)+2]);
      u[i]+= q.dt*(dfx+dfy+dfz);
    }
  }
  free(normal);free(sheet);free(flux);
  return 0;
}
