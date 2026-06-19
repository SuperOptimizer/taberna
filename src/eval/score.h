/* score.h — the official Kaggle Surface-Detection composite score.
 *
 *   Score = 0.30 * TopoScore + 0.35 * SurfaceDice@tau + 0.35 * VOI_score   (tau=2)
 *
 * a linear blend of three bounded [0,1] scores: surface similarity (Surface Dice at
 * tolerance), instance consistency (Variation of Information over surface
 * connected-component instances), and topological correctness (Betti agreement).
 * Computed only over the VALID region (label != ignore_value); predictions in the
 * ignore region are dropped.
 *
 * EXACTNESS: the blend weights (0.30/0.35/0.35) and tau (2.0) are confirmed from
 * the competition. Two sub-score transforms still need the official Evaluation page
 * to be bit-exact and are documented + isolated below:
 *   - VOI_score: raw VOI (bits) -> [0,1]. We use exp(-VOI). (PLACEHOLDER)
 *   - TopoScore: we use a per-dimension Betti-count agreement, a non-spatial
 *     stand-in for the official spatial Betti-matching. (APPROXIMATION)
 * Both are swappable in one place once the official definitions are confirmed.
 */
#ifndef TABERNA_SCORE_H
#define TABERNA_SCORE_H

#include "common/types.h"

typedef struct {
  double surface_dice;   // SurfaceDice@tol, [0,1]
  double voi;            // raw VOI over instances (bits)
  double voi_score;      // bounded [0,1]
  double topo_score;     // [0,1]
  double score;          // 0.30*topo + 0.35*surface_dice + 0.35*voi_score
  // components for inspection
  long pred_b0, pred_b1, pred_b2, gt_b0, gt_b1, gt_b2;
} comp_score;

comp_score competition_score(const u8 *pred, const u8 *label, int nz, int ny, int nx,
                             int tol, u8 surface_value, u8 ignore_value);

#endif // TABERNA_SCORE_H
