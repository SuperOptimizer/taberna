#!/usr/bin/env python3
"""Crop-overlap winding consistency -- a ground-truth-free validation + tiling metric.

Two sheet_sep3d _vol.f32 dumps that cover overlapping world regions must agree on winding
(up to a single global integer offset) wherever they share material voxels. Strong agreement
means the winding is a property of the DATA, not the crop -- and is the precondition for
stitching crops into a full-volume solve.

_vol.f32 layout: int32 header {dz,dy,dx,z0,y0,x0} then dz*dy*dx float32 (NAN off-material).

Usage: overlap_consistency.py A_vol.f32 B_vol.f32 [inner_core_frac=0.15]
The inner umbilicus core (winding ill-defined as r->0, over-sampled) is reported separately.
"""
import sys, numpy as np

def load(fn):
    with open(fn, 'rb') as f:
        dz, dy, dx, z0, y0, x0 = np.fromfile(f, np.int32, 6)
        v = np.fromfile(f, np.float32, dz*dy*dx).reshape(dz, dy, dx)
    return v, (int(z0), int(y0), int(x0)), (int(dz), int(dy), int(dx))

def main():
    A, oA, sA = load(sys.argv[1]); B, oB, sB = load(sys.argv[2])
    core_frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.15
    z0, y0, x0 = (max(oA[i], oB[i]) for i in range(3))
    z1 = min(oA[0]+sA[0], oB[0]+sB[0]); y1 = min(oA[1]+sA[1], oB[1]+sB[1]); x1 = min(oA[2]+sA[2], oB[2]+sB[2])
    if z1 <= z0 or y1 <= y0 or x1 <= x0:
        print("no world overlap"); return
    sa = A[z0-oA[0]:z1-oA[0], y0-oA[1]:y1-oA[1], x0-oA[2]:x1-oA[2]]
    sb = B[z0-oB[0]:z1-oB[0], y0-oB[1]:y1-oB[1], x0-oB[2]:x1-oB[2]]
    m = np.isfinite(sa) & np.isfinite(sb)
    d = sa[m] - sb[m]
    off = np.median(d); r = d - off
    print(f"overlap world z[{z0},{z1}) y[{y0},{y1}) x[{x0},{x1}); {m.sum()} shared material voxels")
    print(f"  global offset={off:.2f} wraps")
    # RAW agreement (no median removed) -- the median subtraction below absorbs ANY constant additive
    # bias including a uniform integer-wrap mislabel, so report raw FIRST. A large |offset| with high raw
    # disagreement but high median-subtracted agreement = the two crops differ by a constant wrap count.
    print(f"  RAW (no-median): agree <0.25: {np.mean(np.abs(d)<0.25):.3f}   <0.5: {np.mean(np.abs(d)<0.5):.3f}   |offset|={abs(off):.2f} wraps")
    print(f"  median-subtracted: agree <0.25 wrap: {np.mean(np.abs(r)<0.25):.3f}   <0.5 wrap: {np.mean(np.abs(r)<0.5):.3f}   tail-std={r.std():.2f}")
    # WARNING: under a shared coarse prior (GLAM) both crops are pulled to the same field, so this metric
    # measures prior coupling, not independent data consistency. Only trust it on INDEPENDENTLY-solved crops.
    # split out the umbilicus core of the overlap region (its own centroid of shared material)
    zz, yy, xx = np.where(m)
    cy, cx = yy.mean(), xx.mean(); rad = np.hypot(yy-cy, xx-cx); rmax = rad.max()
    outer = rad > core_frac*rmax
    ro = r[outer]
    if len(ro):
        print(f"  EXCLUDING inner {core_frac:.0%} core: agree <0.25 wrap: {np.mean(np.abs(ro)<0.25):.3f}   <0.5: {np.mean(np.abs(ro)<0.5):.3f}   tail-std={ro.std():.2f}")

if __name__ == "__main__":
    main()
