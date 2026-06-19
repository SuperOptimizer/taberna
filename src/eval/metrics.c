/* metrics.c — see metrics.h. */
#include "eval/metrics.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MIDX(z, y, x) ((size_t)(z) * nynx + (size_t)(y) * nx + (size_t)(x))

static int is_boundary(const u8 *m, int z, int y, int x, int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  if (!m[MIDX(z, y, x)]) return 0;
  if (x == 0 || x == nx - 1 || y == 0 || y == ny - 1 || z == 0 || z == nz - 1) return 1;
  return !m[MIDX(z, y, x - 1)] || !m[MIDX(z, y, x + 1)] || !m[MIDX(z, y - 1, x)] ||
         !m[MIDX(z, y + 1, x)] || !m[MIDX(z - 1, y, x)] || !m[MIDX(z + 1, y, x)];
}

static int near_boundary(const u8 *other, int z, int y, int x, int tol,
                         int nz, int ny, int nx) {
  size_t nynx = (size_t)ny * nx;
  for (int dz = -tol; dz <= tol; dz++)
    for (int dy = -tol; dy <= tol; dy++)
      for (int dx = -tol; dx <= tol; dx++) {
        int zz = z + dz, yy = y + dy, xx = x + dx;
        if (xx < 0 || xx >= nx || yy < 0 || yy >= ny || zz < 0 || zz >= nz) continue;
        if (is_boundary(other, zz, yy, xx, nz, ny, nx)) return 1;
      }
  return 0;
}

double surface_dice_at_tol(const u8 *seg_mask, const u8 *gt_mask,
                           int nz, int ny, int nx, int tol) {
  size_t ms = 0, mg = 0, hit_s = 0, hit_g = 0;
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
        if (is_boundary(seg_mask, z, y, x, nz, ny, nx)) {
          ms++;
          if (near_boundary(gt_mask, z, y, x, tol, nz, ny, nx)) hit_s++;
        }
        if (is_boundary(gt_mask, z, y, x, nz, ny, nx)) {
          mg++;
          if (near_boundary(seg_mask, z, y, x, tol, nz, ny, nx)) hit_g++;
        }
      }
  double denom = (double)(ms + mg);
  return denom > 0 ? (double)(hit_s + hit_g) / denom : 1.0;
}

surface_eval eval_surface(const u8 *pred, const u8 *label, int nz, int ny, int nx,
                          int tol, u8 surface_value, u8 ignore_value) {
  surface_eval r;
  memset(&r, 0, sizeof r);
  size_t n = (size_t)nz * ny * nx;

  // valid-restricted prediction and GT-surface masks
  u8 *p = (u8 *)malloc(n), *g = (u8 *)malloc(n);
  size_t inter = 0;
  for (size_t i = 0; i < n; i++) {
    int valid = label[i] != ignore_value;
    u8 pi = (valid && pred[i]) ? 1 : 0;
    u8 gi = (label[i] == surface_value) ? 1 : 0;   // surface implies valid
    p[i] = pi; g[i] = gi;
    r.valid_vox += valid;
    r.pred_vox += pi;
    r.gt_vox += gi;
    inter += (pi & gi);
  }
  double denom = (double)(r.pred_vox + r.gt_vox);
  r.dice = denom > 0 ? (2.0 * (double)inter) / denom : 1.0;
  r.surface_dice = surface_dice_at_tol(p, g, nz, ny, nx, tol);

  topo_betti pb = betti_numbers(p, nz, ny, nx);
  r.pred_b0 = pb.b0; r.pred_b1 = pb.b1; r.pred_b2 = pb.b2;
  topo_betti gb = betti_numbers(g, nz, ny, nx);
  r.gt_b0 = gb.b0;

  free(p); free(g);
  return r;
}

// ---- minimal open-addressing hash: u64 key -> u64 count ---------------------

#define H_EMPTY (~(u64)0)

typedef struct {
  u64   *keys;
  u64   *vals;
  size_t cap;   // power of two
  size_t len;
} h64;

static void h64_init(h64 *h, size_t hint) {
  size_t cap = 16;
  while (cap < hint * 2) cap <<= 1;
  h->keys = (u64 *)malloc(cap * sizeof(u64));
  h->vals = (u64 *)calloc(cap, sizeof(u64));
  h->cap = cap;
  h->len = 0;
  for (size_t i = 0; i < cap; i++) h->keys[i] = H_EMPTY;
}
static void h64_free(h64 *h) { free(h->keys); free(h->vals); }

static u64 h64_mix(u64 x) {  // splitmix64 finalizer
  x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27; x *= 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return x;
}

static void h64_grow(h64 *h);

static void h64_inc(h64 *h, u64 key, u64 by) {
  if ((h->len + 1) * 4 >= h->cap * 3) h64_grow(h);
  size_t mask = h->cap - 1;
  size_t i = h64_mix(key) & mask;
  while (h->keys[i] != H_EMPTY && h->keys[i] != key) i = (i + 1) & mask;
  if (h->keys[i] == H_EMPTY) { h->keys[i] = key; h->len++; }
  h->vals[i] += by;
}

static void h64_grow(h64 *h) {
  h64 n;
  size_t cap = h->cap << 1;
  n.keys = (u64 *)malloc(cap * sizeof(u64));
  n.vals = (u64 *)calloc(cap, sizeof(u64));
  n.cap = cap; n.len = 0;
  for (size_t i = 0; i < cap; i++) n.keys[i] = H_EMPTY;
  for (size_t i = 0; i < h->cap; i++)
    if (h->keys[i] != H_EMPTY) h64_inc(&n, h->keys[i], h->vals[i]);
  h64_free(h);
  *h = n;
}

// ---- evaluation -------------------------------------------------------------

eval_result eval_seg(const u32 *seg, const u32 *gt, size_t n) {
  eval_result r = {0, 0, 0, 0, 0, 0};
  if (n == 0) return r;

  h64 joint, ms, mg;
  h64_init(&joint, 1024);
  h64_init(&ms, 256);
  h64_init(&mg, 256);

  for (size_t i = 0; i < n; i++) {
    u64 a = seg[i], b = gt[i];
    h64_inc(&joint, (a << 32) | b, 1);
    h64_inc(&ms, a, 1);
    h64_inc(&mg, b, 1);
  }

  double N = (double)n;
  double I = 0.0, Hs = 0.0, Hg = 0.0;
  double S = 0.0, A2 = 0.0, B2 = 0.0;  // adapted-Rand sums

  for (size_t i = 0; i < joint.cap; i++) {
    if (joint.keys[i] == H_EMPTY) continue;
    double nij = (double)joint.vals[i];
    u32 a = (u32)(joint.keys[i] >> 32);
    u32 b = (u32)(joint.keys[i] & 0xffffffffu);
    // look up marginals
    double ai = 0, bj = 0;
    {
      size_t mask = ms.cap - 1, k = h64_mix(a) & mask;
      while (ms.keys[k] != H_EMPTY && ms.keys[k] != a) k = (k + 1) & mask;
      ai = (double)ms.vals[k];
    }
    {
      size_t mask = mg.cap - 1, k = h64_mix(b) & mask;
      while (mg.keys[k] != H_EMPTY && mg.keys[k] != b) k = (k + 1) & mask;
      bj = (double)mg.vals[k];
    }
    I += (nij / N) * log2((nij * N) / (ai * bj));
    S += nij * nij;
  }
  for (size_t i = 0; i < ms.cap; i++) {
    if (ms.keys[i] == H_EMPTY) continue;
    double a = (double)ms.vals[i];
    Hs -= (a / N) * log2(a / N);
    A2 += a * a;
  }
  for (size_t i = 0; i < mg.cap; i++) {
    if (mg.keys[i] == H_EMPTY) continue;
    double b = (double)mg.vals[i];
    Hg -= (b / N) * log2(b / N);
    B2 += b * b;
  }

  r.voi_split = Hs - I;  // H(seg|gt)
  r.voi_merge = Hg - I;  // H(gt|seg)
  if (r.voi_split < 0) r.voi_split = 0;  // clamp tiny negatives from fp error
  if (r.voi_merge < 0) r.voi_merge = 0;
  r.voi = r.voi_split + r.voi_merge;
  r.precision = S / B2;
  r.recall = S / A2;
  r.are = 1.0 - 2.0 * r.precision * r.recall / (r.precision + r.recall);

  h64_free(&joint);
  h64_free(&ms);
  h64_free(&mg);
  return r;
}
