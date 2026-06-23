#!/usr/bin/env python3
"""leak_metric.py — quantify DEFECT 2 (leaked/fused touches) in a winding field, the measurement the
touching-sheets plan lacked. A LEAK is a continuous-material span ~one pitch long along the outward
normal whose winding stays flat (the field failed to count the turn THROUGH a gapless fused touch) —
distinct from a healthy wrap boundary, where an air gap interrupts the material path.

Method (robust): material-masked in-plane grad U -> smoothed outward normal n; MARCH outward in 2px
steps up to `pitch`, requiring CONTINUOUS material (a fused span); a leak is a fused span whose
end-to-start winding increment dW < ~0.5 (should be ~1). Reports overall + per-radius prevalence.

  python3 scripts/leak_metric.py FIELD_vol.f32 [pitch=40] [umb_cy umb_cx]

Finding (2026-06-23, L1 region z7872 y1436 x3160, TV field): fused spans = 3.2% of material, LEAKS =
~1% (0.66% at dW<0.3, 1.1% at dW<0.5), ~flat across radius. => after Phase 3a (wind_tv) the residual
touch-leak is a ~1% LONG TAIL, not a dominant defect; the heavy LOGISMOS min-separation cut is low-yield
on accessible mid-scroll data. See docs/touching-sheets-plan.md.
"""
import sys, numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

def loadW(f, zc=None):
    h = np.fromfile(f, dtype=np.int32, count=6); dz, dy, dx = int(h[0]), int(h[1]), int(h[2])
    W = np.fromfile(f, dtype=np.float32, offset=24).reshape(dz, dy, dx)
    return W[dz // 2 if zc is None else zc]

def leak_map(U, pitch=40.0):
    mat = np.isfinite(U); Uf = np.nan_to_num(U)
    gy = np.zeros_like(Uf); gx = np.zeros_like(Uf)
    gy[1:-1] = np.where(mat[2:] & mat[:-2], (Uf[2:] - Uf[:-2]) * 0.5, 0)
    gx[:, 1:-1] = np.where(mat[:, 2:] & mat[:, :-2], (Uf[:, 2:] - Uf[:, :-2]) * 0.5, 0)
    m = mat.astype(float)
    sgx = gaussian_filter(gx * m, 3) / (gaussian_filter(m, 3) + 1e-6)
    sgy = gaussian_filter(gy * m, 3) / (gaussian_filter(m, 3) + 1e-6)
    mag = np.hypot(sgx, sgy) + 1e-9; nx, ny = sgx / mag, sgy / mag
    H, W = U.shape; ys, xs = np.mgrid[0:H, 0:W].astype(np.float64)
    yst, xst = ys.copy(), xs.copy(); allmat = mat.copy()
    for _ in range(int(pitch // 2)):
        yst += 2 * ny; xst += 2 * nx
        allmat &= (map_coordinates(m, [yst, xst], order=1, mode='nearest') > 0.9)
    dW = map_coordinates(Uf, [yst, xst], order=1, mode='nearest') - U
    cand = mat & allmat & (np.abs(sgx) + np.abs(sgy) > 0.005)   # in a continuous-material full-pitch span
    return mat, cand, dW

if __name__ == "__main__":
    f = sys.argv[1]; pitch = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
    U = loadW(f); mat, cand, dW = leak_map(U, pitch)
    nm = mat.sum()
    print(f"material={nm}  fused spans={100*cand.sum()/nm:.1f}%")
    for thr in (0.3, 0.5):
        print(f"  LEAK (fused & dW<{thr}): {100*(cand & (dW < thr)).sum()/nm:.2f}% of material")
    if len(sys.argv) > 4:
        cy, cx = float(sys.argv[3]), float(sys.argv[4])
        H, W = U.shape; ys, xs = np.mgrid[0:H, 0:W]; R = np.hypot(ys - cy, xs - cx)
        leak = cand & (dW < 0.5)
        print("  per-radius leak%:")
        for lo, hi in [(0, 300), (300, 500), (500, 700), (700, 1000), (1000, 1600)]:
            b = mat & (R >= lo) & (R < hi)
            if b.sum() > 1000:
                print(f"    r[{lo},{hi}): leak={100*(leak & b).sum()/b.sum():.2f}%")
