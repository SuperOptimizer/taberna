/* snic.c — see snic.h.
 *
 * Ported from https://github.com/spelufo/stabia/tree/main (MIT licensed,
 * modified by VerditeLabs for usage with Hraun), itself based on the reference
 * SNIC implementation and paper:
 *   - https://www.epfl.ch/labs/ivrl/research/snic-superpixels/
 *   - https://github.com/achanta/SNIC/
 *
 * MIT License
 *
 * Copyright (c) 2023 Santiago Pelufo
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "snic.h"

#include <math.h>
#include <stdlib.h>

// HEAP ////////////////////////////////////////////////////////////////////////

// A binary max-heap priority queue. The only SNIC-specific parts are HeapNode
// and heap_node_val: the value is the SNIC distance negated, because this is a
// max-heap and SNIC wants the smallest distance popped first.

typedef struct HeapNode {
  f32 d;
  u32 k;
  u16 x, y, z;
  u16 pad;
} HeapNode;

#define heap_node_val(n) (-n.d)

typedef struct Heap {
  int len, size;
  HeapNode *nodes;
} Heap;

#define heap_left(i)   (2 * (i))
#define heap_right(i)  (2 * (i) + 1)
#define heap_parent(i) ((i) / 2)
#define heap_fix_edge(heap, i, j)                                  \
  if (heap_node_val(heap->nodes[j]) > heap_node_val(heap->nodes[i])) { \
    HeapNode tmp = heap->nodes[j];                                 \
    heap->nodes[j] = heap->nodes[i];                               \
    heap->nodes[i] = tmp;                                          \
  }

static Heap heap_alloc(int size) {
  return (Heap){.len = 0,
                .size = size,
                .nodes = (HeapNode *)calloc((size_t)size * 2 + 1, sizeof(HeapNode))};
}

static void heap_free(Heap *heap) {
  free(heap->nodes);
}

static void heap_push(Heap *heap, HeapNode node) {
  // assert(heap->len <= heap->size);
  heap->len++;
  heap->nodes[heap->len] = node;
  for (int i = heap->len, j = 0; i > 1; i = j) {
    j = heap_parent(i);
    heap_fix_edge(heap, j, i) else break;
  }
}

static HeapNode heap_pop(Heap *heap) {
  // assert(heap->len > 0);
  HeapNode node = heap->nodes[1];
  heap->len--;
  heap->nodes[1] = heap->nodes[heap->len + 1];
  for (int i = 1, j = 0; i <= heap->len; i = j) {
    int l = heap_left(i);
    int r = heap_right(i);
    if (l > heap->len) {
      break;
    }
    j = l;
    if (r <= heap->len && heap_node_val(heap->nodes[l]) < heap_node_val(heap->nodes[r])) {
      j = r;
    }
    heap_fix_edge(heap, i, j) else break;
  }
  return node;
}

#undef heap_left
#undef heap_right
#undef heap_parent
#undef heap_fix_edge

// SNIC ////////////////////////////////////////////////////////////////////////

int snic_superpixel_max_neighs(void) {
  return SUPERPIXEL_MAX_NEIGHS;
}

int superpixel_add_neighbors(Superpixel *superpixels, u32 k1, u32 k2) {
  for (int i = 0; i < SUPERPIXEL_MAX_NEIGHS; i++) {
    if (superpixels[k1].neighs[i] == 0) {
      superpixels[k1].neighs[i] = k2;
      return 0;
    } else if (superpixels[k1].neighs[i] == k2) {
      return 0;
    }
  }
  return 1;
}

int snic_superpixel_count(int nz, int ny, int nx, int d_seed) {
  int cz = (nz - d_seed / 2 + d_seed - 1) / d_seed;
  int cy = (ny - d_seed / 2 + d_seed - 1) / d_seed;
  int cx = (nx - d_seed / 2 + d_seed - 1) / d_seed;
  return cx * cy * cz;
}

// STRUCTURE-TENSOR EXTENSION SEAM: the spiral-undeformation work wants a
// dominant sheet orientation on each supervoxel. The clean place to add it is
// here — accumulate a per-supervoxel 6-vector (Jxx,Jyy,Jzz,Jxy,Jxz,Jyz) of the
// structure tensor alongside c/x/y/z in the pop loop below, then eigendecompose
// per supervoxel in the finalization loop. That keeps the proven distance
// metric untouched; orientation becomes a node attribute, not a clustering
// term. Left out until there's a consumer.
int snic(f32 *img, int nz, int ny, int nx, int d_seed, f32 compactness,
         f32 lowmid, f32 midhig, u32 *labels, Superpixel *superpixels) {
  int neigh_overflow = 0; // Number of adjacencies that couldn't be recorded.
  int nynx = ny * nx;
  int img_size = nynx * nz;
  // taberna layout: z-major, x-fastest. (z,y,x) -> linear, +1=x +nx=y +nynx=z.
  #define idx(z, y, x) ((z) * nynx + (y) * nx + (x))
  #define sqr(x) ((x) * (x))

  // Seed the priority queue on a grid with step d_seed.
  Heap pq = heap_alloc(img_size * 16);
  u32 numk = 0;
  for (u16 iz = d_seed / 2; iz < nz; iz += d_seed) {
    for (u16 iy = d_seed / 2; iy < ny; iy += d_seed) {
      for (u16 ix = d_seed / 2; ix < nx; ix += d_seed) {
        numk++;
        // Nudge the seed toward the lowest-gradient voxel in its 3x3x3
        // neighborhood. Not essential but improves results.
        u16 x = ix, y = iy, z = iz;
        f32 grad = INFINITY;
        for (s16 dz = -1; dz <= 1; dz++) {
          for (s16 dy = -1; dy <= 1; dy++) {
            for (s16 dx = -1; dx <= 1; dx++) {
              int jx = ix + dx, jy = iy + dy, jz = iz + dz;
              if (0 < jx && jx < nx - 1 && 0 < jy && jy < ny - 1 && 0 < jz && jz < nz - 1) {
                f32 gx = img[idx(jz, jy, jx + 1)] - img[idx(jz, jy, jx - 1)];
                f32 gy = img[idx(jz, jy + 1, jx)] - img[idx(jz, jy - 1, jx)];
                f32 gz = img[idx(jz + 1, jy, jx)] - img[idx(jz - 1, jy, jx)];
                f32 jgrad = sqr(gx) + sqr(gy) + sqr(gz);
                if (jgrad < grad) {
                  x = (u16)jx; y = (u16)jy; z = (u16)jz;
                  grad = jgrad;
                }
              }
            }
          }
        }
        heap_push(&pq, (HeapNode){.d = 0.0f, .k = numk, .x = x, .y = y, .z = z});
      }
    }
  }
  // assert(numk == snic_superpixel_count(nz, ny, nx, d_seed));
  if (numk == 0) {
    heap_free(&pq);
    return 0;
  }

  f32 invwt = (compactness * compactness * numk) / (f32)(img_size);

  while (pq.len > 0) {
    HeapNode n = heap_pop(&pq);
    int i = idx(n.z, n.y, n.x);
    if (labels[i] > 0) continue;

    u32 k = n.k;
    labels[i] = k;
    superpixels[k].c += img[i];
    superpixels[k].x += n.x;
    superpixels[k].y += n.y;
    superpixels[k].z += n.z;
    superpixels[k].n += 1;
    if      (img[i] <= lowmid) superpixels[k].nlow += 1;
    else if (img[i] <= midhig) superpixels[k].nmid += 1;
    else                       superpixels[k].nhig += 1;

    #define do_neigh(ndz, ndy, ndx, ioffset)                                  \
      {                                                                        \
        int xx = n.x + (ndx); int yy = n.y + (ndy); int zz = n.z + (ndz);      \
        if (0 <= xx && xx < nx && 0 <= yy && yy < ny && 0 <= zz && zz < nz) {  \
          int ii = i + (ioffset);                                             \
          if (labels[ii] <= 0) {                                              \
            f32 ksize = (f32)superpixels[k].n;                                \
            f32 dc = sqr(100.0f * (superpixels[k].c - (img[ii] * ksize)));     \
            f32 dx = superpixels[k].x - xx * ksize;                           \
            f32 dy = superpixels[k].y - yy * ksize;                           \
            f32 dz = superpixels[k].z - zz * ksize;                           \
            f32 dpos = sqr(dx) + sqr(dy) + sqr(dz);                           \
            f32 d = (dc + dpos * invwt) / (ksize * ksize);                    \
            heap_push(&pq, (HeapNode){.d = d, .k = k, .x = (u16)xx, .y = (u16)yy, .z = (u16)zz}); \
          } else if (k != labels[ii]) {                                       \
            neigh_overflow += superpixel_add_neighbors(superpixels, k, labels[ii]); \
            neigh_overflow += superpixel_add_neighbors(superpixels, labels[ii], k); \
          }                                                                    \
        }                                                                      \
      }

    do_neigh( 0,  0,  1,  1);     // +x
    do_neigh( 0,  0, -1, -1);     // -x
    do_neigh( 0,  1,  0,  nx);    // +y
    do_neigh( 0, -1,  0, -nx);    // -y
    do_neigh( 1,  0,  0,  nynx);  // +z
    do_neigh(-1,  0,  0, -nynx);  // -z
    #undef do_neigh
  }

  for (u32 k = 1; k <= numk; k++) {
    f32 ksize = (f32)superpixels[k].n;
    superpixels[k].c /= ksize;
    superpixels[k].x /= ksize;
    superpixels[k].y /= ksize;
    superpixels[k].z /= ksize;
  }

  #undef sqr
  #undef idx
  heap_free(&pq);
  return neigh_overflow;
}
