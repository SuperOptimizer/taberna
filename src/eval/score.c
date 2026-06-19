/* score.c — see score.h. */
#include "eval/score.h"
#include "eval/metrics.h"
#include "eval/topo.h"

#include <math.h>
#include <stdlib.h>

comp_score competition_score(const u8 *pred, const u8 *label, int nz, int ny, int nx,
                             int tol, u8 surface_value, u8 ignore_value) {
  comp_score r = {0};
  const double VOI_ALPHA = 0.3;   // official: voi_transform='one_over_one_plus', alpha=0.3
  size_t n = (size_t)nz * ny * nx;

  // valid-restricted prediction and GT-surface masks
  u8 *p = (u8 *)malloc(n), *g = (u8 *)malloc(n);
  for (size_t i = 0; i < n; i++) {
    int valid = label[i] != ignore_value;
    p[i] = (valid && pred[i]) ? 1 : 0;
    g[i] = (label[i] == surface_value) ? 1 : 0;
  }

  // --- SurfaceDice@tol (confirmed component) --------------------------------
  r.surface_dice = surface_dice_at_tol(p, g, nz, ny, nx, tol);

  // --- VOI: cc3d(26) instances, scored over UNION-FG voxels only (matches
  //     topometrics voi.py: use_union_mask=True). VOI total is symmetric, so
  //     argument order doesn't matter for the score. -------------------------
  u32 *pl = (u32 *)malloc(n * sizeof(u32)), *gl = (u32 *)malloc(n * sizeof(u32));
  cc_label(p, nz, ny, nx, TOPO_CONN26, pl);
  cc_label(g, nz, ny, nx, TOPO_CONN26, gl);
  u32 *ps = (u32 *)malloc(n * sizeof(u32)), *gs = (u32 *)malloc(n * sizeof(u32));
  size_t m = 0;
  for (size_t i = 0; i < n; i++) {
    if (pl[i] == 0 && gl[i] == 0) continue;   // union of foreground only
    ps[m] = pl[i]; gs[m] = gl[i]; m++;
  }
  eval_result er = eval_seg(gs, ps, m);   // (GT, PR): split=H(GT|PR), merge=H(PR|GT)
  r.voi = er.voi;
  r.voi_score = 1.0 / (1.0 + VOI_ALPHA * er.voi);   // official 'one_over_one_plus'

  // --- TopoScore PROXY: official is weighted Topo-F1 = 2*m_k/(p_k+g_k) over
  //     dims {0,1,2}, where m_k is the SPATIAL Betti-matching matched count from
  //     the Betti-Matching-3D lib (per 8 octant tiles, masks inverted). We have
  //     only Betti *counts*, so we proxy m_k = min(p_k,g_k) (optimistic upper
  //     bound) and apply the official F1 + weights. Use the bridge (scripts/
  //     official_score.py) for the real value. ------------------------------
  topo_betti pb = betti_numbers(p, nz, ny, nx);
  topo_betti gb = betti_numbers(g, nz, ny, nx);
  r.pred_b0=pb.b0; r.pred_b1=pb.b1; r.pred_b2=pb.b2;
  r.gt_b0=gb.b0;   r.gt_b1=gb.b1;   r.gt_b2=gb.b2;
  long pbn[3] = {pb.b0, pb.b1, pb.b2}, gbn[3] = {gb.b0, gb.b1, gb.b2};
  static const double tw[3] = {0.34, 0.33, 0.33};
  double tnum = 0.0, tden = 0.0;
  for (int d = 0; d < 3; d++) {
    long pk = pbn[d], gk = gbn[d], mk = pk < gk ? pk : gk;  // proxy match
    if (pk + gk > 0) {
      double f1 = gk != 0 ? (2.0*mk)/(double)(pk+gk) : 0.5/((double)(pk+gk)+0.5);
      tnum += tw[d] * f1; tden += tw[d];
    }
  }
  r.topo_score = tden > 0 ? tnum / tden : 1.0;

  r.score = 0.30 * r.topo_score + 0.35 * r.surface_dice + 0.35 * r.voi_score;

  free(p); free(g); free(pl); free(gl); free(ps); free(gs);
  return r;
}
