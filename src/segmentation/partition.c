/* partition.c — see partition.h. Mutex Watershed. */
#include "segmentation/partition.h"

#include <math.h>
#include <stdlib.h>

// ---- union-find -------------------------------------------------------------

typedef struct {
  int *parent;
  int *size;
} ufind;

static int uf_find(ufind *u, int x) {
  while (u->parent[x] != x) { u->parent[x] = u->parent[u->parent[x]]; x = u->parent[x]; }
  return x;
}

// ---- per-root mutex lists ---------------------------------------------------

typedef struct {
  int *ids;
  int  n, cap;
} mxlist;

static void mx_push(mxlist *m, int v) {
  if (m->n == m->cap) {
    m->cap = m->cap ? m->cap * 2 : 4;
    m->ids = (int *)realloc(m->ids, (size_t)m->cap * sizeof(int));
  }
  m->ids[m->n++] = v;
}

// is there a mutex between roots ra and rb? (entries may be stale -> resolve)
static int mx_blocked(ufind *u, mxlist *mx, int ra, int rb) {
  mxlist *m = &mx[ra];
  for (int i = 0; i < m->n; i++)
    if (uf_find(u, m->ids[i]) == rb) return 1;
  // check the (usually shorter) other side too, in case of asymmetry
  m = &mx[rb];
  for (int i = 0; i < m->n; i++)
    if (uf_find(u, m->ids[i]) == ra) return 1;
  return 0;
}

// ---- edge sort by descending |w| -------------------------------------------

static int edge_cmp(const void *pa, const void *pb) {
  float wa = fabsf(((const sg_edge *)pa)->w);
  float wb = fabsf(((const sg_edge *)pb)->w);
  if (wa < wb) return 1;   // descending
  if (wa > wb) return -1;
  return 0;
}

int mws_partition(const sgraph *g, u32 *labels) {
  int nn = g->nnodes;
  if (nn <= 0) return 0;

  ufind u;
  u.parent = (int *)malloc((size_t)nn * sizeof(int));
  u.size = (int *)malloc((size_t)nn * sizeof(int));
  mxlist *mx = (mxlist *)calloc((size_t)nn, sizeof(mxlist));
  if (!u.parent || !u.size || !mx) { free(u.parent); free(u.size); free(mx); return -1; }
  for (int i = 0; i < nn; i++) { u.parent[i] = i; u.size[i] = 1; }

  // sort a private copy of the edges by |w| descending
  sg_edge *e = (sg_edge *)malloc((size_t)g->nedges * sizeof(sg_edge));
  if (!e && g->nedges) { free(u.parent); free(u.size); free(mx); return -1; }
  for (int i = 0; i < g->nedges; i++) e[i] = g->edges[i];
  qsort(e, (size_t)g->nedges, sizeof(sg_edge), edge_cmp);

  for (int i = 0; i < g->nedges; i++) {
    int ra = uf_find(&u, (int)e[i].a);
    int rb = uf_find(&u, (int)e[i].b);
    if (ra == rb) continue;
    if (e[i].w > 0.0f) {                 // attractive: merge unless mutexed
      if (mx_blocked(&u, mx, ra, rb)) continue;
      // union by size; keep the larger root, inherit the smaller's mutex list
      int big = ra, small = rb;
      if (u.size[big] < u.size[small]) { big = rb; small = ra; }
      u.parent[small] = big;
      u.size[big] += u.size[small];
      mxlist *ms = &mx[small];
      for (int k = 0; k < ms->n; k++) mx_push(&mx[big], ms->ids[k]);
      free(ms->ids); ms->ids = NULL; ms->n = ms->cap = 0;
    } else if (e[i].w < 0.0f) {          // repulsive: plant a mutex
      mx_push(&mx[ra], rb);
      mx_push(&mx[rb], ra);
    }
  }

  // relabel roots to dense [0, ncluster)
  int *remap = (int *)malloc((size_t)nn * sizeof(int));
  for (int i = 0; i < nn; i++) remap[i] = -1;
  int ncluster = 0;
  for (int i = 0; i < nn; i++) {
    int r = uf_find(&u, i);
    if (remap[r] < 0) remap[r] = ncluster++;
    labels[i] = (u32)remap[r];
  }

  for (int i = 0; i < nn; i++) free(mx[i].ids);
  free(mx); free(remap); free(e); free(u.parent); free(u.size);
  return ncluster;
}
