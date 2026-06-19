/* persist0.c — see persist0.h. Union-find merge tree over ascending voxel order. */
#include "topo/persist0.h"

#include <stdlib.h>

#define IDX(z,y,x) ((size_t)(z)*nynx + (size_t)(y)*nx + (x))

typedef struct { f32 v; int idx; } vk;
static int vk_cmp(const void *a, const void *b) {
  const vk *x = a, *y = b;
  if (x->v < y->v) return -1;
  if (x->v > y->v) return 1;
  return (x->idx > y->idx) - (x->idx < y->idx);   // stable tie-break by index
}

static int uf_find(int *p, int a) { while (p[a]!=a){ p[a]=p[p[a]]; a=p[a]; } return a; }

/* shared core: process voxels ascending, union 6-neighbors already active.
 * If `pairs` != NULL, record off-diagonal dim-0 deaths. If `field_io` != NULL,
 * simplify in place (raise basins with persistence < tau). */
static pers_pair *dim0_core(f32 *field_io, const f32 *field_ro, int nz, int ny, int nx,
                            f32 tau, int *npairs) {
  size_t nynx = (size_t)ny*nx, N = (size_t)nz*nynx;
  const f32 *F = field_io ? field_io : field_ro;
  vk *ord = malloc(N*sizeof(vk));
  for (size_t i=0;i<N;i++){ ord[i].v=F[i]; ord[i].idx=(int)i; }
  qsort(ord, N, sizeof(vk), vk_cmp);

  int *parent = malloc(N*sizeof(int));
  f32 *birth  = malloc(N*sizeof(f32));   // per-root: min value of the basin
  int *bvtx   = malloc(N*sizeof(int));   // per-root: birth voxel (for tie-break)
  char *active= calloc(N,1);
  // for simplification: per-root linked list of member voxels (head/tail/next) to raise
  int *head = field_io ? malloc(N*sizeof(int)) : NULL;
  int *tail = field_io ? malloc(N*sizeof(int)) : NULL;
  int *next = field_io ? malloc(N*sizeof(int)) : NULL;

  pers_pair *pairs = npairs ? malloc((N+1)*sizeof(pers_pair)) : NULL;
  int np = 0;

  const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};

  for (size_t o=0;o<N;o++) {
    int vi = ord[o].idx;
    f32 s = F[vi];
    int z = vi/(int)nynx, r = vi%(int)nynx, y = r/nx, x = r%nx;
    parent[vi]=vi; birth[vi]=s; bvtx[vi]=vi; active[vi]=1;
    if (head){ head[vi]=vi; tail[vi]=vi; next[vi]=-1; }

    for (int k=0;k<6;k++) {
      int zz=z+off[k][0], yy=y+off[k][1], xx=x+off[k][2];
      if (zz<0||zz>=nz||yy<0||yy>=ny||xx<0||xx>=nx) continue;
      int ui = (int)IDX(zz,yy,xx);
      if (!active[ui]) continue;
      int ra = uf_find(parent, vi), rb = uf_find(parent, ui);
      if (ra==rb) continue;
      // elder rule: older = smaller (birth, bvtx)
      int older, younger;
      if (birth[ra] < birth[rb] || (birth[ra]==birth[rb] && bvtx[ra]<bvtx[rb])) { older=ra; younger=rb; }
      else { older=rb; younger=ra; }
      if (pairs && s > birth[younger])   // off-diagonal death of the younger basin
        pairs[np++] = (pers_pair){0, birth[younger], s};
      if (head && (s - birth[younger]) < tau) {           // simplify: raise dying basin
        for (int t=head[younger]; t!=-1; t=next[t]) field_io[t] = s;
      }
      // union younger -> older, keep older's birth; splice member lists (O(1) via tail)
      parent[younger]=older;
      if (head) { next[tail[older]] = head[younger]; tail[older] = tail[younger]; }
      // (vi may merge several neighbors in this loop; root tracked via find)
    }
  }
  // essential classes (surviving roots)
  if (pairs) {
    char *seen = calloc(N,1);
    for (size_t i=0;i<N;i++){ int rt=uf_find(parent,(int)i); if(!seen[rt]){ seen[rt]=1;
      pairs[np++] = (pers_pair){0, birth[rt], TOPO_INF}; } }
    free(seen);
  }
  free(ord); free(parent); free(birth); free(bvtx); free(active);
  free(head); free(tail); free(next);
  if (npairs) *npairs = np;
  return pairs;
}

pers_pair *dim0_persistence(const f32 *field, int nz, int ny, int nx, int *npairs) {
  return dim0_core(NULL, field, nz, ny, nx, 0.0f, npairs);
}

void simplify_dim0(f32 *field, int nz, int ny, int nx, f32 tau) {
  dim0_core(field, NULL, nz, ny, nx, tau, NULL);
}
