#!/usr/bin/env python3
"""Multi-resolution tiled winding -- the WORKING tiling pipeline (cf. tile_winding.py, the negative
diagnostic that proved independent-tile stitching fails).

The fix (sheet_sep3d commit 746e457): a SHARED GLOBAL REFERENCE. Solve a coarse winding ONCE over
the whole region (one crop = one connected component = internally consistent), then run every fine
tile with that coarse field as its prior:
  - priorvol + priorlod + LPRI~8  -> tiles share an ABSOLUTE scale (floating integer offset dies)
  - GLAM~0.3                       -> the dense fill is anchored per-voxel to the coarse field, so
                                      adjacent tiles' DENSE windings agree (~0.97 within 0.25 wrap)
Merge is then TRIVIAL: no per-tile integer offset, no BFS -- just average the overlaps.

Usage: tile_winding_mr.py ARC FINELOD Z0 Y0 X0 DZ TS NY NX STEP COARSELOD OUTDIR [NZ=1 STEPZ=DZ]
  All of Z0/Y0/X0/DZ/TS/STEP are in the FINELOD grid (same convention as sheet_sep3d).
  Tiles: NZ x NY x NX grid; each tile is DZ x TS x TS, stepping STEPZ in z and STEP in y/x.
  COARSELOD: the (coarser) LOD at which the one global field is solved, e.g. FINELOD+1.
  With NZ>1 the pipeline stacks z-bands too -> a tall merged volume for a readable unroll.
"""
import sys, os, subprocess, numpy as np

TOOL = os.path.join(os.path.dirname(__file__), "..", "build", "tools", "sheet_sep3d")
LPRI, GLAM = "8", "0.3"          # node-prior strength / dense-fill global-anchor strength

def load(fn):
    with open(fn, "rb") as f:
        dz, dy, dx, z0, y0, x0 = np.fromfile(f, np.int32, 6)
        v = np.fromfile(f, np.float32, int(dz)*int(dy)*int(dx)).reshape(int(dz), int(dy), int(dx))
    return v, (int(z0), int(y0), int(x0))

def run(arc, lod, z0, y0, x0, dz, dy, dx, out, pitch, prior=None, plod=None):
    # ARC OUT lod z0 y0 x0 dz dy dx minseg pitch LAM WEQ ztie recto LPRI radw barclose zmed SP
    #   umbref wdrescue robsig usesheet slicecv cxoff [priorvol priorlod GLAM]
    a = [TOOL, arc, out, str(lod), str(z0), str(y0), str(x0), str(dz), str(dy), str(dx),
         "40", str(pitch), "0.6", "3", "0", "0",
         LPRI if prior else "0.15", "1.0", "1", "0", "0", "1", "0", "0.6", "1", "0", "0"]
    if prior:
        a += [prior, str(plod), GLAM]
    env = dict(os.environ, OMP_NUM_THREADS="8", ASAN_OPTIONS="detect_leaks=0")
    p = subprocess.run(a, capture_output=True, text=True, env=env)
    if not os.path.exists(f"{out}_vol.f32"):
        sys.stderr.write(p.stderr[-2000:]); raise SystemExit(f"tile {out} failed")
    return p.stderr

def main():
    a = sys.argv[1:]
    arc, flod, z0, y0, x0, dz, ts, ny, nx, step, clod, outdir = a[:12]
    flod, z0, y0, x0, dz, ts, ny, nx, step, clod = map(int, (flod, z0, y0, x0, dz, ts, ny, nx, step, clod))
    nz   = int(a[12]) if len(a) > 12 else 1
    stepz = int(a[13]) if len(a) > 13 else dz
    os.makedirs(outdir, exist_ok=True)
    # union bbox in FINELOD coords
    uz1, uy1, ux1 = z0 + (nz-1)*stepz + dz, y0 + (ny-1)*step + ts, x0 + (nx-1)*step + ts
    pitch_f = 20                                  # native-L2 pitch; scaled per-LOD below
    # --- 1) ONE coarse global solve over the FULL 3D union, at COARSELOD ---
    s = 1 << (clod - flod)                         # fine->coarse downsample factor
    cz0, cy0, cx0 = z0 // s, y0 // s, x0 // s
    cdz, cdy, cdx = max(2, (uz1 - z0) // s), (uy1 - y0) // s, (ux1 - x0) // s
    cg = f"{outdir}/coarse"
    print(f"coarse solve @ LOD{clod}: z{cz0} y{cy0} x{cx0} + {cdz}x{cdy}x{cdx}")
    run(arc, clod, cz0, cy0, cx0, cdz, cdy, cdx, cg, max(4, pitch_f // s))
    cgvol = f"{cg}_vol.f32"
    # --- 2) fine tiles (NZ x NY x NX), each anchored to the coarse field ---
    tiles = {}
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                tz, ty, tx = z0 + iz*stepz, y0 + iy*step, x0 + ix*step
                out = f"{outdir}/t_{iz}_{iy}_{ix}"
                run(arc, flod, tz, ty, tx, dz, ts, ts, out, pitch_f, prior=cgvol, plod=clod)
                tiles[(iz, iy, ix)] = load(f"{out}_vol.f32")
        print(f"z-band {iz} done")
    # --- 3) seam agreement (DIRECT, no integer offset), across all 3 axes ---
    seams = []
    for (iz, iy, ix), (v0, o0) in tiles.items():
        for diz, diy, dix in ((0, 0, 1), (0, 1, 0), (1, 0, 0)):
            nb = (iz+diz, iy+diy, ix+dix)
            if nb not in tiles:
                continue
            v1, o1 = tiles[nb]
            za, zb = max(o0[0], o1[0]), min(o0[0]+v0.shape[0], o1[0]+v1.shape[0])
            ya, yb = max(o0[1], o1[1]), min(o0[1]+v0.shape[1], o1[1]+v1.shape[1])
            xa, xb = max(o0[2], o1[2]), min(o0[2]+v0.shape[2], o1[2]+v1.shape[2])
            if zb <= za or yb <= ya or xb <= xa:
                continue
            aa = v0[za-o0[0]:zb-o0[0], ya-o0[1]:yb-o0[1], xa-o0[2]:xb-o0[2]]
            bb = v1[za-o1[0]:zb-o1[0], ya-o1[1]:yb-o1[1], xa-o1[2]:xb-o1[2]]
            m = np.isfinite(aa) & np.isfinite(bb)
            if m.sum() < 100:
                continue
            d = aa[m] - bb[m]; md = float(np.median(d))
            seams.append((f"{iz},{iy},{ix}->{nb[0]},{nb[1]},{nb[2]}", md,
                          float(np.mean(np.abs(d - md) < 0.25)), float((d - md).std()), int(m.sum())))
    print("\nseams (median-offset, agree<0.25, tail-std, voxels):")
    for s_ in seams:
        print(f"  {s_[0]}: off={s_[1]:+.2f} agree={s_[2]:.3f} tailstd={s_[3]:.2f} n={s_[4]}")
    if seams:
        print(f"MEAN seam agreement <0.25 wrap = {np.mean([s_[2] for s_ in seams]):.3f}  "
              f"mean|offset|={np.mean([abs(s_[1]) for s_ in seams]):.3f}")
    # --- 4) merge: direct average over overlaps (offsets are ~0) ---
    Z0 = min(o[0] for _, o in tiles.values()); Y0 = min(o[1] for _, o in tiles.values()); X0 = min(o[2] for _, o in tiles.values())
    Z1 = max(o[0]+v.shape[0] for v, o in tiles.values())
    Y1 = max(o[1]+v.shape[1] for v, o in tiles.values()); X1 = max(o[2]+v.shape[2] for v, o in tiles.values())
    D, H, W = Z1-Z0, Y1-Y0, X1-X0
    acc = np.zeros((D, H, W), np.float64); cnt = np.zeros((D, H, W), np.int32)
    for v, o in tiles.values():
        fin = np.isfinite(v); sz, sy, sx = o[0]-Z0, o[1]-Y0, o[2]-X0
        acc[sz:sz+v.shape[0], sy:sy+v.shape[1], sx:sx+v.shape[2]] += np.where(fin, v, 0)
        cnt[sz:sz+v.shape[0], sy:sy+v.shape[1], sx:sx+v.shape[2]] += fin
    merged = np.where(cnt > 0, acc/np.maximum(cnt, 1), np.nan).astype(np.float32)
    print(f"\nmerged {merged.shape}, covered frac {np.mean(cnt>0):.2f}, range "
          f"{np.nanmin(merged):.1f}..{np.nanmax(merged):.1f} wraps")
    with open(f"{outdir}/merged.f32", "wb") as f:
        np.array([D, H, W, Z0, Y0, X0], np.int32).tofile(f); merged.tofile(f)
    print(f"wrote {outdir}/merged.f32")

if __name__ == "__main__":
    main()
