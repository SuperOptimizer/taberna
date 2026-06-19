/* affinity.c — see affinity.h. */
#include "segmentation/affinity.h"

#include <math.h>
#include <stdlib.h>

#define IDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

affinity_params affinity_default_params(void) {
  return (affinity_params){.k_attract = 1.0f, .k_repel = 1.0f};
}

void sgraph_free(sgraph *g) {
  if (g && g->edges) { free(g->edges); g->edges = NULL; g->nedges = 0; }
}

int affinity_build_voxel(const u8 *mask, const f32 *sheet, const f32 *normal,
                         int nz, int ny, int nx, const affinity_params *p_in,
                         sgraph *g, s32 **node_of_out, u32 **voxel_of_out) {
  affinity_params p = p_in ? *p_in : affinity_default_params();
  size_t n = (size_t)nz * ny * nx;
  size_t nynx = (size_t)ny * nx;

  s32 *node_of = (s32 *)malloc(n * sizeof(s32));
  if (!node_of) return -1;

  // assign dense node ids to foreground voxels
  u32 nnodes = 0;
  for (size_t i = 0; i < n; i++) node_of[i] = mask[i] ? (s32)(nnodes++) : -1;

  u32 *voxel_of = (u32 *)malloc((size_t)(nnodes ? nnodes : 1) * sizeof(u32));
  if (!voxel_of) { free(node_of); return -1; }
  for (size_t i = 0; i < n; i++) if (node_of[i] >= 0) voxel_of[node_of[i]] = (u32)i;

  // edges: +x/+y/+z neighbor contacts among foreground voxels
  int cap = 1024, ne = 0;
  sg_edge *edges = (sg_edge *)malloc((size_t)cap * sizeof(sg_edge));
  if (!edges) { free(node_of); free(voxel_of); return -1; }

  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        size_t i = IDX(z, y, x);
        if (node_of[i] < 0) continue;
        // (neighbor index in +dir, normal component axis)
        for (int ax = 0; ax < 3; ax++) {
          int zz = z + (ax == 2), yy = y + (ax == 1), xx = x + (ax == 0);
          if (xx >= nx || yy >= ny || zz >= nz) continue;
          size_t j = IDX(zz, yy, xx);
          if (node_of[j] < 0) continue;
          double across = 0.5 * (fabs((double)normal[3 * i + ax]) +
                                 fabs((double)normal[3 * j + ax]));
          double sb = 0.5 * ((double)sheet[i] + (double)sheet[j]);
          double w = p.k_attract * (1.0 - across) - p.k_repel * (across * sb);
          if (ne == cap) {
            cap *= 2;
            sg_edge *t = (sg_edge *)realloc(edges, (size_t)cap * sizeof(sg_edge));
            if (!t) { free(edges); free(node_of); free(voxel_of); return -1; }
            edges = t;
          }
          edges[ne].a = (u32)node_of[i];
          edges[ne].b = (u32)node_of[j];
          edges[ne].w = (f32)w;
          ne++;
        }
      }
    }
  }

  g->edges = edges;
  g->nedges = ne;
  g->nnodes = (int)nnodes;
  if (node_of_out) *node_of_out = node_of; else free(node_of);
  if (voxel_of_out) *voxel_of_out = voxel_of; else free(voxel_of);
  return 0;
}
