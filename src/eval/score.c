/* score.c — see score.h. */
#include "eval/score.h"
#include "eval/metrics.h"
#include "eval/topo.h"
#include "eval/nsd.h"

#include <math.h>
#include <stdlib.h>

/* reduced 6-conn Betti (rb0=b0-1 clamped, rb1=b1, rb2=b2) of a cropped sub-box */
static void reduced_betti6_crop(const u8 *mask, int nz, int ny, int nx,
                                int z0, int z1, int y0, int y1, int x0, int x1,
                                long rb[3]) {
  int dz = z1-z0, dy = y1-y0, dx = x1-x0;
  size_t nynx = (size_t)ny * nx;
  u8 *sub = (u8 *)malloc((size_t)dz*dy*dx);
  for (int z=0; z<dz; z++)
    for (int y=0; y<dy; y++)
      for (int x=0; x<dx; x++)
        sub[((size_t)z*dy + y)*dx + x] = mask[(size_t)(z0+z)*nynx + (size_t)(y0+y)*nx + (x0+x)];
  topo_betti b = betti_numbers_6(sub, dz, dy, dx);
  rb[0] = b.b0 > 0 ? b.b0 - 1 : 0;  // reduced H0 (drop the essential component)
  rb[1] = b.b1;
  rb[2] = b.b2;
  free(sub);
}

/* PROXY TopoScore. The official matched count m_k is the IMAGE-INDUCED Betti
 * matching (features of pred and gt that map to the same class in the comparison
 * image) — NOT inclusion-exclusion on Betti numbers: rb(pred)+rb(gt)-rb(union)
 * overcounts whenever an induced map has a kernel (a tunnel in pred filled by the
 * union), which is common on dense masks. We therefore CLAMP the estimate to the
 * valid range [0, min(p_k,g_k)]; this is a bounded proxy (exact only when features
 * are well separated). For the exact value use scripts/official_score.py, or
 * taberna's real Betti matching once built on src/topo/cubical.c. */
/* EXACT dim-0 matched count on a cropped tile via union-find:
 * m_0 = |{union-comps hit by pred} ∩ {union-comps hit by gt}| - 1  (clamped >=0).
 * Validated against the Betti-Matching oracle (60/60). */
static long matched0_crop(const u8 *pred, const u8 *gt, const u8 *u,
                          int nz, int ny, int nx,
                          int z0, int z1, int y0, int y1, int x0, int x1) {
  int dz=z1-z0, dy=y1-y0, dx=x1-x0;
  size_t nynx=(size_t)ny*nx, m=(size_t)dz*dy*dx;
  u8 *su=malloc(m), *sp=malloc(m), *sg=malloc(m);
  for (int z=0;z<dz;z++) for (int y=0;y<dy;y++) for (int x=0;x<dx;x++) {
    size_t s=((size_t)z*dy+y)*dx+x, g=(size_t)(z0+z)*nynx+(size_t)(y0+y)*nx+(x0+x);
    su[s]=u[g]; sp[s]=pred[g]; sg[s]=gt[g];
  }
  u32 *ul=malloc(m*sizeof(u32));
  u32 cu=cc_label(su,dz,dy,dx,TOPO_CONN6,ul);
  u8 *inA=calloc(cu+1,1), *inB=calloc(cu+1,1);
  for (size_t i=0;i<m;i++){ if(sp[i]&&ul[i]) inA[ul[i]]=1; if(sg[i]&&ul[i]) inB[ul[i]]=1; }
  long inter=0; for(u32 c=1;c<=cu;c++) if(inA[c]&&inB[c]) inter++;
  free(su);free(sp);free(sg);free(ul);free(inA);free(inB);
  return inter>0 ? inter-1 : 0;
}

double toposcore_native(const u8 *pred, const u8 *gt, int nz, int ny, int nx) {
  size_t n = (size_t)nz * ny * nx;
  u8 *u = (u8 *)malloc(n);
  for (size_t i = 0; i < n; i++) u[i] = (pred[i] || gt[i]) ? 1 : 0;

  // official _axis_cuts(n,2): first chunk gets the ceil half
  int zc = (nz+1)/2, yc = (ny+1)/2, xc = (nx+1)/2;
  int zb[3] = {0, zc, nz}, yb[3] = {0, yc, ny}, xb[3] = {0, xc, nx};

  long M[3] = {0,0,0}, P[3] = {0,0,0}, G[3] = {0,0,0};
  for (int zi=0; zi<2; zi++) for (int yi=0; yi<2; yi++) for (int xi=0; xi<2; xi++) {
    int z0=zb[zi], z1=zb[zi+1], y0=yb[yi], y1=yb[yi+1], x0=xb[xi], x1=xb[xi+1];
    if (z1<=z0 || y1<=y0 || x1<=x0) continue;
    long rp[3], rg[3], ru[3];
    reduced_betti6_crop(pred, nz,ny,nx, z0,z1,y0,y1,x0,x1, rp);
    reduced_betti6_crop(gt,   nz,ny,nx, z0,z1,y0,y1,x0,x1, rg);
    reduced_betti6_crop(u,    nz,ny,nx, z0,z1,y0,y1,x0,x1, ru);
    P[0]+=rp[0]; G[0]+=rg[0];
    M[0]+= matched0_crop(pred,gt,u, nz,ny,nx, z0,z1,y0,y1,x0,x1);   // EXACT dim-0
    for (int k=1;k<3;k++){                                         // dim 1,2 PROXY
      P[k]+=rp[k]; G[k]+=rg[k];
      long m = rp[k]+rg[k]-ru[k];
      long hi = rp[k]<rg[k]?rp[k]:rg[k];
      if (m < 0) m = 0; if (m > hi) m = hi;
      M[k]+= m;
    }
  }
  free(u);

  static const double w[3] = {0.34, 0.33, 0.33};
  double num = 0.0, den = 0.0;
  for (int k=0;k<3;k++) {
    long denom = P[k] + G[k];
    if (denom <= 0) continue;  // inactive dim
    double f1 = G[k] != 0 ? (2.0*(double)M[k])/(double)denom
                          : 0.5/((double)denom + 0.5);
    num += w[k]*f1; den += w[k];
  }
  return den > 0 ? num/den : 1.0;
}

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

  // --- SurfaceDice@tol — exact native NSD (Google surface_distance port) -----
  r.surface_dice = surface_dice_nsd(g, p, nz, ny, nx, (double)tol);

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

  // --- TopoScore: EXACT native (tiled binary Betti reduction; matches the
  //     official Betti-Matching TopoScore for binary masks) -------------------
  r.topo_score = toposcore_native(p, g, nz, ny, nx);
  topo_betti pb = betti_numbers_6(p, nz, ny, nx);
  topo_betti gb = betti_numbers_6(g, nz, ny, nx);
  r.pred_b0=pb.b0; r.pred_b1=pb.b1; r.pred_b2=pb.b2;
  r.gt_b0=gb.b0;   r.gt_b1=gb.b1;   r.gt_b2=gb.b2;

  r.score = 0.30 * r.topo_score + 0.35 * r.surface_dice + 0.35 * r.voi_score;

  free(p); free(g); free(pl); free(gl); free(ps); free(gs);
  return r;
}
