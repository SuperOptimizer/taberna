/* stitch.c — see stitch.h. */
#include "segmentation/stitch.h"

#include <stdlib.h>

int stitch_overlap(const u32 *a_over, const u32 *b_over, size_t n,
                   int na_labels, int nb_labels, u32 *remap_b) {
  if (nb_labels <= 0) return na_labels;

  // For each B label, find the A label it most co-occurs with in the overlap.
  // best_a[b] = argmax_a count(a,b); best_cnt tracks the max.
  int *best_a = (int *)malloc((size_t)nb_labels * sizeof(int));
  size_t *best_cnt = (size_t *)calloc((size_t)nb_labels, sizeof(size_t));
  size_t *cur_cnt = (size_t *)calloc((size_t)nb_labels, sizeof(size_t));
  int *cur_a = (int *)malloc((size_t)nb_labels * sizeof(int));
  for (int b = 0; b < nb_labels; b++) { best_a[b] = -1; cur_a[b] = -1; }

  // Simple approach: for each overlap voxel, tally; track per-B the running best.
  // (A full contingency table would be na*nb; we keep it light with a running
  // best per B, recomputed by scanning — fine for v0 / modest label counts.)
  for (size_t i = 0; i < n; i++) {
    int a = (int)a_over[i], b = (int)b_over[i];
    if (b < 0 || b >= nb_labels || a < 0 || a >= na_labels) continue;
    if (cur_a[b] == a) {
      cur_cnt[b]++;
    } else {
      // restart running count for this B label's current candidate
      cur_a[b] = a;
      cur_cnt[b] = 1;
    }
    if (cur_cnt[b] > best_cnt[b]) { best_cnt[b] = cur_cnt[b]; best_a[b] = a; }
  }

  int next = na_labels;
  for (int b = 0; b < nb_labels; b++) {
    if (best_a[b] >= 0) remap_b[b] = (u32)best_a[b];   // merge into the A label
    else remap_b[b] = (u32)(next++);                   // B-only label -> fresh id
  }

  free(best_a); free(best_cnt); free(cur_cnt); free(cur_a);
  return next;
}
