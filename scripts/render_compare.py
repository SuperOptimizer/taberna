#!/usr/bin/env python3
"""Render xy/xz/yz center-slice 3-way comparisons (raw | conservative | aggressive)
from the local region.zarr + conservative.zarr + aggressive.zarr produced by the
tuning run. One PNG per plane, plus a combined contact sheet.

Usage: render_compare.py <scroll>   (default PHerc0332)
"""
import json, sys
import numpy as np
from pathlib import Path
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent.parent
scroll = sys.argv[1] if len(sys.argv) > 1 else "PHerc0332"
TUNE = ROOT / "data/tune" / scroll
# optional argv[2]: comma list of label=subdir (e.g. "raw=region.zarr,aggr=aggressive.zarr")
if len(sys.argv) > 2:
    VERSIONS = [(s.split("=")[0], s.split("=")[1]) for s in sys.argv[2].split(",")]
else:
    VERSIONS = [("raw", "region.zarr"), ("conservative", "conservative.zarr"),
                ("aggressive", "aggressive.zarr")]

def load_zarr(d):
    za = json.loads((d / "0" / ".zarray").read_text())
    sz, sy, sx = za["shape"]; cz, cy, cx = za["chunks"]
    vol = np.zeros((sz, sy, sx), np.uint8)
    nz, ny, nx = sz // cz, sy // cy, sx // cx
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                f = d / "0" / str(iz) / str(iy) / str(ix)
                if f.exists():
                    vol[iz*cz:(iz+1)*cz, iy*cy:(iy+1)*cy, ix*cx:(ix+1)*cx] = \
                        np.frombuffer(f.read_bytes(), np.uint8).reshape(cz, cy, cx)
    return vol

def center_slices(vol):
    z, y, x = (s // 2 for s in vol.shape)
    return {"xy": vol[z, :, :], "xz": vol[:, y, :], "yz": vol[:, :, x]}

# load one volume at a time (peak ~1 GB), keep only the 3 center slices
slices = {}
for name, sub in VERSIONS:
    print(f"loading {sub} ...", flush=True)
    vol = load_zarr(TUNE / sub)
    slices[name] = center_slices(vol)
    del vol

def to_img(arr, lo, hi):
    a = np.clip((arr.astype(np.float32) - lo) / max(hi - lo, 1e-6), 0, 1)
    return (a * 255).astype(np.uint8)

BAR = 28
def panel(plane):
    # common display window per plane from the FIRST version (honest before/after)
    raw = slices[VERSIONS[0][0]][plane]
    lo, hi = np.percentile(raw, 1), np.percentile(raw, 99)
    imgs = []
    for name, _ in VERSIONS:
        g = to_img(slices[name][plane], lo, hi)
        im = Image.fromarray(g, "L").convert("RGB")
        canvas = Image.new("RGB", (im.width, im.height + BAR), (0, 0, 0))
        canvas.paste(im, (0, BAR))
        d = ImageDraw.Draw(canvas)
        d.text((6, 7), f"{plane}  {name}  [win {lo:.0f}-{hi:.0f}]", fill=(255, 255, 0))
        imgs.append(canvas)
    sep = 4
    W = sum(i.width for i in imgs) + sep * (len(imgs) - 1)
    H = imgs[0].height
    row = Image.new("RGB", (W, H), (40, 40, 40))
    x = 0
    for i in imgs:
        row.paste(i, (x, 0)); x += i.width + sep
    return row

rows = []
for plane in ("xy", "xz", "yz"):
    row = panel(plane)
    out = TUNE / f"compare_{plane}.png"
    row.save(out)
    print(f"wrote {out}  ({row.width}x{row.height})")
    rows.append(row)

# combined contact sheet (3 planes stacked)
W = max(r.width for r in rows); H = sum(r.height for r in rows) + 4 * (len(rows) - 1)
sheet = Image.new("RGB", (W, H), (40, 40, 40))
y = 0
for r in rows:
    sheet.paste(r, (0, y)); y += r.height + 4
sheet.save(TUNE / "compare_all.png")
print(f"wrote {TUNE/'compare_all.png'}  ({sheet.width}x{sheet.height})")
