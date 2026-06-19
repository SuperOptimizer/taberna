/* nsd.c — see nsd.h. Faithful port of surface_distance for spacing (1,1,1). */
#include "eval/nsd.h"

#include <stdlib.h>
#include <math.h>

/* marching-cubes surfel area per 2x2x2 neighbour code, spacing (1,1,1).
 * Generated from surface_distance.lookup_tables.create_table_neighbour_code_to_surface_area. */
static const float NSD_AREA[256] = {
  0.000000f, 0.216506f, 0.216506f, 0.707107f, 0.216506f, 0.707107f, 0.433013f, 1.149519f,
  0.216506f, 0.433013f, 0.707107f, 1.149519f, 0.707107f, 1.149519f, 1.149519f, 1.000000f,
  0.216506f, 0.707107f, 0.433013f, 1.149519f, 0.433013f, 1.149519f, 0.649519f, 1.299038f,
  0.433013f, 0.923613f, 0.923613f, 1.573132f, 0.923613f, 1.573132f, 1.366025f, 1.149519f,
  0.216506f, 0.433013f, 0.707107f, 1.149519f, 0.433013f, 0.923613f, 0.923613f, 1.573132f,
  0.433013f, 0.649519f, 1.149519f, 1.299038f, 0.923613f, 1.366025f, 1.573132f, 1.149519f,
  0.707107f, 1.149519f, 1.149519f, 1.000000f, 0.923613f, 1.573132f, 1.366025f, 1.149519f,
  0.923613f, 1.366025f, 1.573132f, 1.149519f, 1.414214f, 0.923613f, 0.923613f, 0.707107f,
  0.216506f, 0.433013f, 0.433013f, 0.923613f, 0.707107f, 1.149519f, 0.923613f, 1.573132f,
  0.433013f, 0.649519f, 0.923613f, 1.366025f, 1.149519f, 1.299038f, 1.573132f, 1.149519f,
  0.707107f, 1.149519f, 0.923613f, 1.573132f, 1.149519f, 1.000000f, 1.366025f, 1.149519f,
  0.923613f, 1.366025f, 1.414214f, 0.923613f, 1.573132f, 1.149519f, 0.923613f, 0.707107f,
  0.433013f, 0.649519f, 0.923613f, 1.366025f, 0.923613f, 1.366025f, 1.414214f, 0.923613f,
  0.649519f, 0.866025f, 1.366025f, 0.649519f, 1.366025f, 0.649519f, 0.923613f, 0.433013f,
  1.149519f, 1.299038f, 1.573132f, 1.149519f, 1.573132f, 1.149519f, 0.923613f, 0.707107f,
  1.366025f, 0.649519f, 0.923613f, 0.433013f, 0.923613f, 0.433013f, 0.433013f, 0.216506f,
  0.216506f, 0.433013f, 0.433013f, 0.923613f, 0.433013f, 0.923613f, 0.649519f, 1.366025f,
  0.707107f, 0.923613f, 1.149519f, 1.573132f, 1.149519f, 1.573132f, 1.299038f, 1.149519f,
  0.433013f, 0.923613f, 0.649519f, 1.366025f, 0.649519f, 1.366025f, 0.866025f, 0.649519f,
  0.923613f, 1.414214f, 1.366025f, 0.923613f, 1.366025f, 0.923613f, 0.649519f, 0.433013f,
  0.707107f, 0.923613f, 1.149519f, 1.573132f, 0.923613f, 1.414214f, 1.366025f, 0.923613f,
  1.149519f, 1.366025f, 1.000000f, 1.149519f, 1.573132f, 0.923613f, 1.149519f, 0.707107f,
  1.149519f, 1.573132f, 1.299038f, 1.149519f, 1.366025f, 0.923613f, 0.649519f, 0.433013f,
  1.573132f, 0.923613f, 1.149519f, 0.707107f, 0.923613f, 0.433013f, 0.433013f, 0.216506f,
  0.707107f, 0.923613f, 0.923613f, 1.414214f, 1.149519f, 1.573132f, 1.366025f, 0.923613f,
  1.149519f, 1.366025f, 1.573132f, 0.923613f, 1.000000f, 1.149519f, 1.149519f, 0.707107f,
  1.149519f, 1.573132f, 1.366025f, 0.923613f, 1.299038f, 1.149519f, 0.649519f, 0.433013f,
  1.573132f, 0.923613f, 0.923613f, 0.433013f, 1.149519f, 0.707107f, 0.433013f, 0.216506f,
  1.149519f, 1.366025f, 1.573132f, 0.923613f, 1.573132f, 0.923613f, 0.923613f, 0.433013f,
  1.299038f, 0.649519f, 1.149519f, 0.433013f, 1.149519f, 0.433013f, 0.707107f, 0.216506f,
  1.000000f, 1.149519f, 1.149519f, 0.707107f, 1.149519f, 0.707107f, 0.433013f, 0.216506f,
  1.149519f, 0.433013f, 0.707107f, 0.216506f, 0.707107f, 0.216506f, 0.216506f, 0.000000f,
};

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))

/* 2x2x2 neighbour code at each voxel-corner (kernel weights 128..1 over the
 * block [z-1..z][y-1..y][x-1..x], out-of-bounds = 0). */
static void code_map(const u8 *m, int nz, int ny, int nx, u8 *code) {
  size_t nynx = (size_t)ny * nx;
  #define MM(z,y,x) (((z)>=0&&(y)>=0&&(x)>=0)? (m[IDX(z,y,x)]!=0):0)
  for (int z=0; z<nz; z++) for (int y=0; y<ny; y++) for (int x=0; x<nx; x++) {
    int c = 128*MM(z-1,y-1,x-1) + 64*MM(z-1,y-1,x) + 32*MM(z-1,y,x-1) + 16*MM(z-1,y,x)
          +   8*MM(z,  y-1,x-1) +  4*MM(z,  y-1,x) +  2*MM(z,  y,x-1) +  1*MM(z,  y,x);
    code[IDX(z,y,x)] = (u8)c;
  }
  #undef MM
}

/* 1D squared-distance transform (Felzenszwalb-Huttenlocher) over f -> d. */
static void dt1d(const float *f, int n, float *d, int *v, float *zb) {
  int k = 0; v[0] = 0; zb[0] = -1e20f; zb[1] = 1e20f;
  for (int q=1; q<n; q++) {
    float s;
    for (;;) {
      s = ((f[q]+(float)q*q) - (f[v[k]]+(float)v[k]*v[k])) / (float)(2*q - 2*v[k]);
      if (s > zb[k]) break;
      k--;
    }
    k++; v[k]=q; zb[k]=s; zb[k+1]=1e20f;
  }
  k = 0;
  for (int q=0; q<n; q++) {
    while (zb[k+1] < (float)q) k++;
    float dq = (float)(q - v[k]);
    d[q] = dq*dq + f[v[k]];
  }
}

/* exact squared EDT: distance^2 from each voxel to the nearest site (sites[i]!=0). */
static void edt2(const u8 *sites, int nz, int ny, int nx, float *out) {
  size_t nynx = (size_t)ny*nx, N = (size_t)nz*nynx;
  for (size_t i=0;i<N;i++) out[i] = sites[i] ? 0.0f : 1e20f;
  int maxn = nz>ny?(nz>nx?nz:nx):(ny>nx?ny:nx);
  float *f = malloc(maxn*sizeof(float)), *d = malloc(maxn*sizeof(float)), *zb = malloc((maxn+1)*sizeof(float));
  int *v = malloc(maxn*sizeof(int));
  // along x
  for (int z=0;z<nz;z++) for (int y=0;y<ny;y++) {
    for (int x=0;x<nx;x++) f[x]=out[IDX(z,y,x)];
    dt1d(f,nx,d,v,zb);
    for (int x=0;x<nx;x++) out[IDX(z,y,x)]=d[x];
  }
  // along y
  for (int z=0;z<nz;z++) for (int x=0;x<nx;x++) {
    for (int y=0;y<ny;y++) f[y]=out[IDX(z,y,x)];
    dt1d(f,ny,d,v,zb);
    for (int y=0;y<ny;y++) out[IDX(z,y,x)]=d[y];
  }
  // along z
  for (int y=0;y<ny;y++) for (int x=0;x<nx;x++) {
    for (int z=0;z<nz;z++) f[z]=out[IDX(z,y,x)];
    dt1d(f,nz,d,v,zb);
    for (int z=0;z<nz;z++) out[IDX(z,y,x)]=d[z];
  }
  free(f);free(d);free(zb);free(v);
}

double surface_dice_nsd(const u8 *gt, const u8 *pred, int nz, int ny, int nx, double tol) {
  size_t nynx = (size_t)ny*nx, N = (size_t)nz*nynx;
  // empty-mask conventions (match leaderboard)
  size_t ng=0, npx=0;
  for (size_t i=0;i<N;i++){ ng += gt[i]!=0; npx += pred[i]!=0; }
  if (ng==0 && npx==0) return 1.0;
  if ((ng==0) != (npx==0)) return 0.0;

  u8 *cg = malloc(N), *cp = malloc(N);
  code_map(gt, nz,ny,nx, cg);
  code_map(pred, nz,ny,nx, cp);
  u8 *bg = malloc(N), *bp = malloc(N);
  for (size_t i=0;i<N;i++){ bg[i] = (cg[i]!=0 && cg[i]!=255); bp[i] = (cp[i]!=0 && cp[i]!=255); }

  float *dist_to_pred = malloc(N*sizeof(float)), *dist_to_gt = malloc(N*sizeof(float));
  edt2(bp, nz,ny,nx, dist_to_pred);   // sq dist to nearest pred-surface
  edt2(bg, nz,ny,nx, dist_to_gt);     // sq dist to nearest gt-surface

  double t2 = tol*tol, num = 0.0, den = 0.0;
  for (size_t i=0;i<N;i++) {
    if (bg[i]) { double a = NSD_AREA[cg[i]]; den += a; if (dist_to_pred[i] <= t2) num += a; }
    if (bp[i]) { double a = NSD_AREA[cp[i]]; den += a; if (dist_to_gt[i]   <= t2) num += a; }
  }
  free(cg);free(cp);free(bg);free(bp);free(dist_to_pred);free(dist_to_gt);
  return den > 0 ? num/den : 1.0;
}
