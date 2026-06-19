/* morph.c — see morph.h. */
#include "postproc/morph.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

void majority_filter_thresh(const u8 *in, u8 *out, int nz, int ny, int nx,
                            int iters, int on_thresh) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  u8 *a = (u8 *)malloc(n), *b = (u8 *)malloc(n);
  memcpy(a, in, n);
  for (int it = 0; it < iters; it++) {
    #pragma omp parallel for schedule(static)
    for (int z = 0; z < nz; z++)
      for (int y = 0; y < ny; y++)
        for (int x = 0; x < nx; x++) {
          int cnt = 0;
          for (int dz = -1; dz <= 1; dz++) {
            int zz = z + dz; if (zz < 0 || zz >= nz) continue;
            for (int dy = -1; dy <= 1; dy++) {
              int yy = y + dy; if (yy < 0 || yy >= ny) continue;
              for (int dx = -1; dx <= 1; dx++) {
                int xx = x + dx; if (xx < 0 || xx >= nx) continue;
                cnt += a[IDX(zz, yy, xx)] != 0;
              }
            }
          }
          b[IDX(z, y, x)] = cnt >= on_thresh;
        }
    u8 *t = a; a = b; b = t;   // a holds latest
  }
  memcpy(out, a, n);
  free(a); free(b);
}

void majority_filter(const u8 *in, u8 *out, int nz, int ny, int nx, int iters) {
  majority_filter_thresh(in, out, nz, ny, nx, iters, 14);
}

size_t remove_small_components(u8 *mask, int nz, int ny, int nx,
                               topo_conn conn, size_t min_voxels) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  u32 *lab = (u32 *)malloc(n * sizeof(u32));
  u32 nc = cc_label(mask, nz, ny, nx, conn, lab);
  if (!nc) { free(lab); return 0; }
  size_t *sz = (size_t *)calloc(nc + 1, sizeof(size_t));
  for (size_t i = 0; i < n; i++) sz[lab[i]]++;
  u8 *keep = (u8 *)calloc(nc + 1, 1);
  size_t removed = 0;
  for (u32 c = 1; c <= nc; c++) { if (sz[c] >= min_voxels) keep[c] = 1; else removed++; }
  for (size_t i = 0; i < n; i++) if (lab[i] && !keep[lab[i]]) mask[i] = 0;
  free(lab); free(sz); free(keep);
  return removed;
}

void fill_holes(u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  // flood-fill background (6-conn) from the border; any bg not reached is a cavity.
  u8 *reached = (u8 *)calloc(n, 1);
  size_t *stack = (size_t *)malloc(n * sizeof(size_t));
  size_t sp = 0;
  #define PUSH(z, y, x) do { size_t _i = IDX(z, y, x); \
    if (!mask[_i] && !reached[_i]) { reached[_i] = 1; stack[sp++] = _i; } } while (0)
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++)
        if (z == 0 || z == nz-1 || y == 0 || y == ny-1 || x == 0 || x == nx-1)
          PUSH(z, y, x);
  while (sp) {
    size_t i = stack[--sp];
    int z = (int)(i / nynx), r = (int)(i % nynx), y = r / nx, x = r % nx;
    if (z > 0)    PUSH(z-1, y, x);
    if (z < nz-1) PUSH(z+1, y, x);
    if (y > 0)    PUSH(z, y-1, x);
    if (y < ny-1) PUSH(z, y+1, x);
    if (x > 0)    PUSH(z, y, x-1);
    if (x < nx-1) PUSH(z, y, x+1);
  }
  #undef PUSH
  for (size_t i = 0; i < n; i++) if (!mask[i] && !reached[i]) mask[i] = 1;
  free(reached); free(stack);
}

size_t plug_pinholes(u8 *mask, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  size_t plugged = 0;
  // interior only (a face-incomplete border voxel can't be fully 6-enclosed)
  for (int z = 1; z < nz-1; z++)
    for (int y = 1; y < ny-1; y++)
      for (int x = 1; x < nx-1; x++) {
        size_t i = IDX(z, y, x);
        if (mask[i]) continue;
        if (mask[IDX(z-1,y,x)] && mask[IDX(z+1,y,x)] && mask[IDX(z,y-1,x)] &&
            mask[IDX(z,y+1,x)] && mask[IDX(z,y,x-1)] && mask[IDX(z,y,x+1)]) {
          mask[i] = 1; plugged++;
        }
      }
  return plugged;
}

/* ---- ball morphology ------------------------------------------------------ */
typedef struct { int dz, dy, dx; } off3;

static off3 *ball_offsets(int r, int *count) {
  int cap = 0;
  for (int dz = -r; dz <= r; dz++)
    for (int dy = -r; dy <= r; dy++)
      for (int dx = -r; dx <= r; dx++)
        if (dz*dz + dy*dy + dx*dx <= r*r) cap++;
  off3 *o = (off3 *)malloc(cap * sizeof(off3));
  int k = 0;
  for (int dz = -r; dz <= r; dz++)
    for (int dy = -r; dy <= r; dy++)
      for (int dx = -r; dx <= r; dx++)
        if (dz*dz + dy*dy + dx*dx <= r*r) { o[k].dz=dz; o[k].dy=dy; o[k].dx=dx; k++; }
  *count = cap;
  return o;
}

static void morph_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r, int erode) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  int no; off3 *o = ball_offsets(r, &no);
  u8 *res = (u8 *)malloc(n);
  #pragma omp parallel for schedule(static)
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        int val;
        if (erode) {  // erode: off iff any SE voxel is bg (border counts as bg)
          val = 1;
          for (int k = 0; k < no; k++) {
            int zz=z+o[k].dz, yy=y+o[k].dy, xx=x+o[k].dx;
            if (zz<0||zz>=nz||yy<0||yy>=ny||xx<0||xx>=nx || !in[IDX(zz,yy,xx)]) { val=0; break; }
          }
        } else {       // dilate: on iff any SE voxel is fg
          val = 0;
          for (int k = 0; k < no; k++) {
            int zz=z+o[k].dz, yy=y+o[k].dy, xx=x+o[k].dx;
            if (zz<0||zz>=nz||yy<0||yy>=ny||xx<0||xx>=nx) continue;
            if (in[IDX(zz,yy,xx)]) { val=1; break; }
          }
        }
        res[IDX(z,y,x)] = (u8)val;
      }
  memcpy(out, res, n);
  free(res); free(o);
}

void dilate_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r) {
  morph_ball(in, out, nz, ny, nx, r, 0);
}
void erode_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r) {
  morph_ball(in, out, nz, ny, nx, r, 1);
}
void closing_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r) {
  size_t n = (size_t)nz * ny * nx;
  u8 *tmp = (u8 *)malloc(n);
  morph_ball(in, tmp, nz, ny, nx, r, 0);   // dilate
  morph_ball(tmp, out, nz, ny, nx, r, 1);  // erode
  free(tmp);
}
void opening_ball(const u8 *in, u8 *out, int nz, int ny, int nx, int r) {
  size_t n = (size_t)nz * ny * nx;
  u8 *tmp = (u8 *)malloc(n);
  morph_ball(in, tmp, nz, ny, nx, r, 1);   // erode
  morph_ball(tmp, out, nz, ny, nx, r, 0);  // dilate
  free(tmp);
}

/* one in-plane step: dilate (grow=1) turns on a bg voxel reached from a fg neighbor
 * through an in-plane offset; erode (grow=0) turns off a fg voxel that has an
 * in-plane bg neighbor. `normal` indexed per the *moving* voxel's anchor: for
 * dilation the fg source neighbor, for erosion the fg voxel itself. */
static void inplane_step(const u8 *in, u8 *out, const f32 *normal,
                         int nz, int ny, int nx, f32 plane_tol, int grow) {
  size_t nynx = (size_t)ny * nx, n = (size_t)nz * nynx;
  for (size_t i = 0; i < n; i++) out[i] = in[i];
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (grow ? (in[i] != 0) : (in[i] == 0)) continue;  // dilate scans bg, erode scans fg
        for (int dz = -1; dz <= 1; dz++)
          for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
              if (!dz && !dy && !dx) continue;
              int zz = z+dz, yy = y+dy, xx = x+dx;
              if (zz<0||zz>=nz||yy<0||yy>=ny||xx<0||xx>=nx) continue;
              size_t j = IDX(zz, yy, xx);
              // dilate: neighbor must be fg (source); erode: neighbor must be bg
              if (grow ? (in[j] == 0) : (in[j] != 0)) continue;
              // in-plane test using the fg anchor's normal (j for dilate, i for erode)
              size_t a = grow ? j : i;
              f32 nx_ = normal[3*a+0], ny_ = normal[3*a+1], nz_ = normal[3*a+2];
              f32 dlen = sqrtf((f32)(dz*dz + dy*dy + dx*dx));
              f32 dot = (dx*nx_ + dy*ny_ + dz*nz_) / dlen;
              if (fabsf(dot) <= plane_tol) { out[i] = grow ? 1 : 0; break; }
            }
      }
}

void inplane_close(u8 *mask, const f32 *normal, int nz, int ny, int nx,
                   f32 plane_tol, int iters) {
  size_t n = (size_t)nz * ny * nx;
  u8 *tmp = (u8 *)malloc(n);
  for (int it = 0; it < iters; it++) { inplane_step(mask, tmp, normal, nz, ny, nx, plane_tol, 1); memcpy(mask, tmp, n); }
  for (int it = 0; it < iters; it++) { inplane_step(mask, tmp, normal, nz, ny, nx, plane_tol, 0); memcpy(mask, tmp, n); }
  free(tmp);
}
