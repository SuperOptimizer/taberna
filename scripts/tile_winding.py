#!/usr/bin/env python3
"""Tile a region into overlapping crops, solve winding per tile, and stitch into one field.

The per-tile winding is consistent across tiles ONLY if they share (a) one global umbilicus and
(b) one pitch. So we run sheet_sep3d with umbref=0 (every crop uses the same coarse-LOD5 world
center) and a fixed pitch. Stitching then needs only a single INTEGER offset per tile, found from
the median winding difference in each overlap (relative to a reference tile, propagated by BFS).
Merge = offset-aligned average over all tiles covering each world voxel.

FINDING (this driver's first result): independent-tile stitching DOES NOT WORK. On properly-filled
fields adjacent crops agree only ~6% within 0.25 wrap (the earlier "85%" was measured on the FLAT
fields from the robust-diffusion bug -- agreement on zero). A single integer offset per tile can't
align them: each crop assigns DIFFERENT integer offsets to the same physical wraps (floating graph
components + per-crop rmin/dense-fill), and a stronger radius prior barely helps (LPRI 4 -> 11%).
The winding is NOT crop-invariant. The correct fix is a SHARED GLOBAL REFERENCE: a coarse global
winding field used as the per-node prior (multi-resolution -- already the design of the older
tools/unroll pipeline), or a single global solve. Kept as the diagnostic that proves the problem.

Usage: tile_winding.py ARC LOD Z0 Y0 X0 DZ TS NY NX STEP OUTDIR
  tiles a NY x NX grid of TS-square crops (single z-band of DZ) stepping STEP in y and x.
"""
import sys, os, subprocess, numpy as np

TOOL = os.path.join(os.path.dirname(__file__), "..", "build", "tools", "sheet_sep3d")

def load(fn):
    with open(fn, "rb") as f:
        dz, dy, dx, z0, y0, x0 = np.fromfile(f, np.int32, 6)
        v = np.fromfile(f, np.float32, int(dz)*int(dy)*int(dx)).reshape(int(dz), int(dy), int(dx))
    return v, (int(z0), int(y0), int(x0))

def run_tile(arc, lod, z0, y0, x0, dz, ts, pitch, out):
    # ...minseg pitch LAM WEQ ztie recto LPRI radw barclose zmed SP umbref wdrescue robsig usesheet slicecv cxoff
    args = [TOOL, arc, out, str(lod), str(z0), str(y0), str(x0), str(dz), str(ts), str(ts),
            "40", str(pitch), "0.6", "3", "0", "0", "0.15", "1.0", "1", "0", "0",
            "0",       # umbref OFF -> global consistent center
            "0", "0.6", "1", "0", "0"]
    env = dict(os.environ, OMP_NUM_THREADS="8", ASAN_OPTIONS="detect_leaks=0")
    p = subprocess.run(args, capture_output=True, text=True, env=env)
    pit = None
    for line in p.stderr.splitlines():
        if line.startswith("pitch="):
            pit = float(line.split()[0].split("=")[1])
    return pit

def main():
    arc, lod, z0, y0, x0, dz, ts, ny, nx, step, outdir = sys.argv[1:12]
    lod, z0, y0, x0, dz, ts, ny, nx, step = map(int, (lod, z0, y0, x0, dz, ts, ny, nx, step))
    os.makedirs(outdir, exist_ok=True)
    # 1) global pitch from the centre tile (auto), reused for all
    cy0, cx0 = y0 + (ny//2)*step, x0 + (nx//2)*step
    pitch = run_tile(arc, lod, z0, cy0, cx0, dz, ts, 0, f"{outdir}/probe")
    print(f"global pitch = {pitch}")
    # 2) run all tiles with the fixed global pitch + global center
    tiles = {}
    for iy in range(ny):
        for ix in range(nx):
            ty, tx = y0 + iy*step, x0 + ix*step
            out = f"{outdir}/t_{iy}_{ix}"
            run_tile(arc, lod, z0, ty, tx, dz, ts, pitch, out)
            v, org = load(f"{out}_vol.f32")
            tiles[(iy, ix)] = (v, org)
            print(f"tile {iy},{ix} origin {org} dims {v.shape}")
    # 3) stitch: integer offset per tile via BFS from (0,0), median diff in overlap
    off = {(0, 0): 0}
    order = [(0, 0)]
    seams = []
    while order:
        iy, ix = order.pop(0)
        v0, o0 = tiles[(iy, ix)]
        for diy, dix in ((0, 1), (1, 0), (0, -1), (-1, 0)):
            nb = (iy+diy, ix+dix)
            if nb not in tiles or nb in off:
                continue
            v1, o1 = tiles[nb]
            # world overlap
            za, zb = max(o0[0], o1[0]), min(o0[0]+v0.shape[0], o1[0]+v1.shape[0])
            ya, yb = max(o0[1], o1[1]), min(o0[1]+v0.shape[1], o1[1]+v1.shape[1])
            xa, xb = max(o0[2], o1[2]), min(o0[2]+v0.shape[2], o1[2]+v1.shape[2])
            if zb <= za or yb <= ya or xb <= xa:
                continue
            a = v0[za-o0[0]:zb-o0[0], ya-o0[1]:yb-o0[1], xa-o0[2]:xb-o0[2]]
            b = v1[za-o1[0]:zb-o1[0], ya-o1[1]:yb-o1[1], xa-o1[2]:xb-o1[2]]
            m = np.isfinite(a) & np.isfinite(b)
            if m.sum() < 100:
                continue
            d = (a[m] + off[(iy, ix)]) - b[m]
            o = int(round(np.median(d)))
            off[nb] = o
            r = d - o
            agree = float(np.mean(np.abs(r) < 0.25))
            seams.append((f"{iy},{ix}->{nb[0]},{nb[1]}", o, agree, int(m.sum())))
            order.append(nb)
    print("\nseams (offset, agree<0.25, voxels):")
    for s in seams:
        print(f"  {s[0]}: off={s[1]:+d} agree={s[2]:.3f} n={s[3]}")
    if seams:
        print(f"mean seam agreement <0.25 wrap = {np.mean([s[2] for s in seams]):.3f}")
    # 4) merge into one world array (offset-aligned average)
    Y0 = min(o[1] for _, o in tiles.values()); X0 = min(o[2] for _, o in tiles.values())
    Y1 = max(o[1]+v.shape[1] for v, o in tiles.values()); X1 = max(o[2]+v.shape[2] for v, o in tiles.values())
    H, W = Y1-Y0, X1-X0
    acc = np.zeros((dz, H, W), np.float64); cnt = np.zeros((dz, H, W), np.int32)
    for k, (v, o) in tiles.items():
        if k not in off:
            continue
        fin = np.isfinite(v)
        sy, sx = o[1]-Y0, o[2]-X0
        acc[:, sy:sy+v.shape[1], sx:sx+v.shape[2]] += np.where(fin, v + off[k], 0)
        cnt[:, sy:sy+v.shape[1], sx:sx+v.shape[2]] += fin
    merged = np.where(cnt > 0, acc/np.maximum(cnt, 1), np.nan).astype(np.float32)
    cov = float(np.mean(cnt > 0))
    print(f"\nmerged {merged.shape}, covered-of-bbox frac {cov:.2f}, tiles stitched {len(off)}/{len(tiles)}")
    with open(f"{outdir}/merged.f32", "wb") as f:
        np.array([dz, H, W, z0, Y0, X0], np.int32).tofile(f); merged.tofile(f)
    print(f"wrote {outdir}/merged.f32")

if __name__ == "__main__":
    main()
