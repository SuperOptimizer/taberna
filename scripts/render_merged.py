#!/usr/bin/env python3
"""Render a merged winding volume (merged.f32): mid-z spiral + a z-y reslice (z-coherence check).
Usage: render_merged.py merged.f32 OUTDIR"""
import sys, numpy as np
from PIL import Image

def jet(t, fin):
    img = np.zeros(t.shape + (3,), np.uint8)
    img[..., 0] = (255*np.clip(1.5-np.abs(4*t-3), 0, 1)).astype(np.uint8)
    img[..., 1] = (255*np.clip(1.5-np.abs(4*t-2), 0, 1)).astype(np.uint8)
    img[..., 2] = (255*np.clip(1.5-np.abs(4*t-1), 0, 1)).astype(np.uint8)
    img[~fin] = 0
    return img

def main():
    fn, outd = sys.argv[1], sys.argv[2]
    with open(fn, "rb") as f:
        D, H, W, Z0, Y0, X0 = np.fromfile(f, np.int32, 6)
        v = np.fromfile(f, np.float32, int(D)*int(H)*int(W)).reshape(int(D), int(H), int(W))
    lo, hi = np.nanpercentile(v, 1), np.nanpercentile(v, 99)
    print(f"merged {v.shape} z0={Z0} winding {np.nanmin(v):.1f}..{np.nanmax(v):.1f} (clip {lo:.1f}..{hi:.1f})")
    # mid-z spiral
    mid = v[int(D)//2]; fin = np.isfinite(mid)
    t = np.clip((np.where(fin, mid, lo)-lo)/(hi-lo), 0, 1)
    Image.fromarray(jet(t, fin)).save(f"{outd}/merged_midz.png")
    # z-y reslice at mid-x: vertical=z (stacked bands), horizontal=y. A clean winding is smooth in z
    # -> horizontal bands continuous top-to-bottom across z-tile seams (no horizontal cut lines).
    sx = int(W)//2; resl = v[:, :, sx]; finr = np.isfinite(resl)
    tr = np.clip((np.where(finr, resl, lo)-lo)/(hi-lo), 0, 1)
    img = jet(tr, finr)
    scale = max(1, 600//int(D))
    Image.fromarray(img).resize((int(H), int(D)*scale), Image.NEAREST).save(f"{outd}/merged_zy.png")
    print(f"wrote {outd}/merged_midz.png ({H}x{W}) and {outd}/merged_zy.png (z-y reslice, {D} z rows x{scale})")

if __name__ == "__main__":
    main()
