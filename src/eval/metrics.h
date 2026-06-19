/* metrics.h — segmentation evaluation for the unwrapping experiments.
 *
 * Operates on two equal-length integer label volumes (seg = our result, gt =
 * ground truth). Labels are arbitrary u32 (label 0 is just another label, e.g.
 * background); the metrics are permutation-invariant.
 *
 *   Variation of Information (in bits): voi_split = H(seg|gt) (false SPLITS, one
 *   region cut into many) and voi_merge = H(gt|seg) (false MERGES, distinct
 *   regions fused). For the scroll a merge across adjacent wraps is catastrophic
 *   and a split within a wrap is recoverable, so report them separately and weight
 *   voi_merge more heavily downstream (see docs/unwrapping-plan.md §6).
 *
 *   Adapted-Rand (CREMI style): precision (normalized by GT cluster sizes),
 *   recall (by seg cluster sizes), and are = 1 - 2*p*r/(p+r).
 */
#ifndef TABERNA_METRICS_H
#define TABERNA_METRICS_H

#include <stddef.h>
#include "common/types.h"
#include "eval/topo.h"

typedef struct {
  double voi_split;   // H(seg|gt), bits  (false splits / over-segmentation)
  double voi_merge;   // H(gt|seg), bits  (false merges / under-segmentation)
  double voi;         // voi_split + voi_merge
  double precision;   // sum n_ij^2 / sum gt_j^2
  double recall;      // sum n_ij^2 / sum seg_i^2
  double are;         // adapted Rand error, 1 - 2pr/(p+r); 0 = perfect
} eval_result;

eval_result eval_seg(const u32 *seg, const u32 *gt, size_t n);

// Surface Dice at tolerance: agreement of the two foreground-surface (boundary)
// voxel sets, counting a boundary voxel as matched if a boundary voxel of the
// other mask lies within Chebyshev distance `tol`. Returns [0,1], 1 = perfect.
// (v0 brute-force within a (2tol+1)^3 window — fine for experiment cubes.)
double surface_dice_at_tol(const u8 *seg_mask, const u8 *gt_mask,
                           int nz, int ny, int nx, int tol);

/* Ignore-aware surface-detection evaluation against a 3-class competition label
 * (`label`: 0 = background, `surface_value` = papyrus surface, `ignore_value` =
 * not-evaluated). Voxels equal to ignore_value are excluded from every measure
 * (predictions there are zeroed first). Reports the volumetric Dice and the
 * boundary Surface Dice@tol over the valid region, plus the topology of the
 * (valid-restricted) prediction — the things the Kaggle metric blends. The exact
 * official composite additionally uses persistent-homology Betti matching and
 * instance VOI; those are future work (see docs/surface-detection-kaggle.md). */
typedef struct {
  double dice;          // 2|P∩G| / (|P|+|G|) over valid voxels
  double surface_dice;  // boundary agreement @tol over valid region
  long   pred_b0, pred_b1, pred_b2;  // topology of the prediction
  long   gt_b0;         // # connected components of the GT surface (context)
  size_t pred_vox, gt_vox, valid_vox;
} surface_eval;

surface_eval eval_surface(const u8 *pred, const u8 *label, int nz, int ny, int nx,
                          int tol, u8 surface_value, u8 ignore_value);

#endif // TABERNA_METRICS_H
