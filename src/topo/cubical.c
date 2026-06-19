/* cubical.c — see cubical.h. Generic Z/2 boundary-matrix reduction over the
 * sublevel T-construction cubical complex. Correctness-first (small volumes). */
#include "topo/cubical.h"

#include <stdlib.h>
#include <string.h>

/* global cell id layout over N = nz*ny*nx voxels:
 *   vertex v            : v                       in [0, N)
 *   edge(axis a, v)     : N + a*N + v             a: 0=x,1=y,2=z
 *   square(normal a, v) : 4N + a*N + v            a: 0=x(yz),1=y(xz),2=z(xy)
 *   cube(v)             : 7N + v
 * (ids may name non-existent cells; only generated cells get a filtration slot.) */

typedef struct { float val; int dim; int gid; } cell;

static int cell_cmp(const void *pa, const void *pb) {
  const cell *a = pa, *b = pb;
  if (a->val < b->val) return -1;
  if (a->val > b->val) return 1;
  if (a->dim != b->dim) return a->dim - b->dim;     // faces before cofaces
  return (a->gid > b->gid) - (a->gid < b->gid);
}

/* growable int vector (a sparse Z/2 column = sorted row positions) */
typedef struct { int *a; int n, cap; } ivec;
static void iv_push(ivec *v, int x) {
  if (v->n == v->cap) { v->cap = v->cap ? v->cap*2 : 8; v->a = realloc(v->a, v->cap*sizeof(int)); }
  v->a[v->n++] = x;
}
/* symmetric difference of two SORTED vectors into dst (overwrites dst) */
static void iv_xor(ivec *dst, const ivec *src) {
  int i = 0, j = 0; ivec out = {0,0,0};
  while (i < dst->n && j < src->n) {
    if (dst->a[i] < src->a[j]) iv_push(&out, dst->a[i++]);
    else if (dst->a[i] > src->a[j]) iv_push(&out, src->a[j++]);
    else { i++; j++; }                                  // equal -> cancels
  }
  while (i < dst->n) iv_push(&out, dst->a[i++]);
  while (j < src->n) iv_push(&out, src->a[j++]);
  free(dst->a); *dst = out;
}

pers_pair *cubical_persistence(const f32 *field, int nz, int ny, int nx, int *npairs) {
  *npairs = 0;
  size_t N = (size_t)nz * ny * nx;
  const int xs = 1, ys = nx, zs = nx * ny;
  #define V(z,y,x) ((size_t)(z)*zs + (size_t)(y)*ys + (x))
  #define MAX2(a,b) ((a)>(b)?(a):(b))

  // --- 1) enumerate existing cells with filtration values ---
  // full cubical complex has up to 8N cells: N vertices + ~3N edges + ~3N squares + ~N cubes
  size_t cap = 8 * N + 16;
  cell *cells = malloc(cap * sizeof(cell));
  size_t nc = 0;
  #define ADD(gid_,dim_,val_) do{ cells[nc].gid=(int)(gid_); cells[nc].dim=(dim_); cells[nc].val=(val_); nc++; }while(0)
  for (int z=0; z<nz; z++) for (int y=0; y<ny; y++) for (int x=0; x<nx; x++) {
    size_t v = V(z,y,x); float fv = field[v];
    ADD(v, 0, fv);                                                  // vertex
    if (x+1<nx) ADD(N + 0*N + v, 1, MAX2(fv, field[v+xs]));         // edge x
    if (y+1<ny) ADD(N + 1*N + v, 1, MAX2(fv, field[v+ys]));         // edge y
    if (z+1<nz) ADD(N + 2*N + v, 1, MAX2(fv, field[v+zs]));         // edge z
    if (y+1<ny&&x+1<nx){ float m=MAX2(MAX2(fv,field[v+xs]),MAX2(field[v+ys],field[v+xs+ys])); ADD(4*N+2*N+v,2,m);} // sq normal z
    if (z+1<nz&&x+1<nx){ float m=MAX2(MAX2(fv,field[v+xs]),MAX2(field[v+zs],field[v+xs+zs])); ADD(4*N+1*N+v,2,m);} // sq normal y
    if (z+1<nz&&y+1<ny){ float m=MAX2(MAX2(fv,field[v+ys]),MAX2(field[v+zs],field[v+ys+zs])); ADD(4*N+0*N+v,2,m);} // sq normal x
    if (z+1<nz&&y+1<ny&&x+1<nx){
      float m=fv; m=MAX2(m,field[v+xs]); m=MAX2(m,field[v+ys]); m=MAX2(m,field[v+zs]);
      m=MAX2(m,field[v+xs+ys]); m=MAX2(m,field[v+xs+zs]); m=MAX2(m,field[v+ys+zs]); m=MAX2(m,field[v+xs+ys+zs]);
      ADD(7*N+v,3,m); }
  }

  // --- 2) filtration order; pos[gid] = order index ---
  qsort(cells, nc, sizeof(cell), cell_cmp);
  int *pos = malloc(8*N*sizeof(int));
  for (size_t i=0;i<8*N;i++) pos[i] = -1;
  for (size_t i=0;i<nc;i++) pos[cells[i].gid] = (int)i;

  // --- 3) reduce ---
  ivec *col = calloc(nc, sizeof(ivec));        // reduced column per cell (if pivot)
  int  *low2col = malloc(nc*sizeof(int));      // row position -> column that pivots there
  for (size_t i=0;i<nc;i++) low2col[i] = -1;
  char *killed = calloc(nc, 1);
  pers_pair *pairs = malloc(nc * sizeof(pers_pair)); int np = 0;

  for (size_t j=0;j<nc;j++) {
    int gid = cells[j].gid, dim = cells[j].dim;
    size_t base = gid % N;            // anchor voxel of the cell (gid layout is +k*N)
    ivec c = {0,0,0};
    // build boundary (faces) as filtration positions
    #define BD(faceid) do{ int p = pos[(faceid)]; if(p>=0) iv_push(&c,p); }while(0)
    if (dim==1) {
      int axis = (gid - (int)N) / (int)N;  // 0..2
      size_t v = base;
      BD(v);                                            // vertex v
      BD(v + (axis==0?xs:axis==1?ys:zs));               // vertex v+e_axis
    } else if (dim==2) {
      int normal = (gid - 4*(int)N) / (int)N;           // 0=x,1=y,2=z
      size_t v = base;
      if (normal==2) { BD(N+0*N+v); BD(N+0*N+v+ys); BD(N+1*N+v); BD(N+1*N+v+xs); }       // xy: x,y edges
      else if (normal==1){ BD(N+0*N+v); BD(N+0*N+v+zs); BD(N+2*N+v); BD(N+2*N+v+xs); }   // xz: x,z edges
      else { BD(N+1*N+v); BD(N+1*N+v+zs); BD(N+2*N+v); BD(N+2*N+v+ys); }                 // yz: y,z edges
    } else if (dim==3) {
      size_t v = base;
      BD(4*N+2*N+v); BD(4*N+2*N+v+zs);    // z-normal faces
      BD(4*N+1*N+v); BD(4*N+1*N+v+ys);    // y-normal faces
      BD(4*N+0*N+v); BD(4*N+0*N+v+xs);    // x-normal faces
    }
    #undef BD
    // sort the (small) boundary list
    for (int p=1;p<c.n;p++){ int t=c.a[p],q=p-1; while(q>=0&&c.a[q]>t){c.a[q+1]=c.a[q];q--;} c.a[q+1]=t; }

    // reduce
    while (c.n) {
      int l = c.a[c.n-1];
      if (low2col[l] < 0) break;
      iv_xor(&c, &col[low2col[l]]);
    }
    if (c.n) {
      int l = c.a[c.n-1];
      low2col[l] = (int)j;
      col[j] = c;                        // keep for future XORs
      killed[l] = 1;
      if (cells[j].val > cells[l].val) {  // off-diagonal => real feature
        pairs[np].dim = cells[l].dim; pairs[np].birth = cells[l].val; pairs[np].death = cells[j].val; np++;
      }
    } else {
      free(c.a);                          // birth (empty column)
    }
  }
  // essential classes: a cell is a birth iff its reduced column was empty; it is
  // essential iff it is a birth and never killed (never a death's pivot row).
  for (size_t i=0;i<nc;i++)
    if (col[i].n==0 && !killed[i]) {
      pairs[np].dim = cells[i].dim; pairs[np].birth = cells[i].val; pairs[np].death = TOPO_INF; np++;
    }

  for (size_t i=0;i<nc;i++) free(col[i].a);
  free(col); free(low2col); free(killed); free(pos); free(cells);
  #undef V
  #undef MAX2
  #undef ADD
  *npairs = np;
  return pairs;
}
