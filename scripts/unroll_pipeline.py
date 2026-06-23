#!/usr/bin/env python3
"""End-to-end classical unroll pipeline: region -> coarse winding -> fine tiles -> merge -> unroll.

One reproducible entry point chaining the validated pieces:
  1. scripts/tile_winding_mr.py : coarse global winding + GLAM-anchored fine tiles -> OUTDIR/merged.f32
  2. tools/unroll_wind (auto umbilicus center) : merged winding -> OUTDIR/unroll.pgm (spiral flatten)
     optionally FINE mode: sample the winding against a finer archive for native-res output.

Usage:
  unroll_pipeline.py ARC FINELOD Z0 Y0 X0 DZ TS NY NX STEP COARSELOD OUTDIR
                     [NZ STEPZ] [--samp N] [--fine FINEARC SCALE ZA ZN]
All region args are the tile_winding_mr convention (FINELOD-grid).
  --fine: render the (L2/L1) winding against FINEARC (e.g. the L1 or L0 archive) at SCALE x, over the
          winding-z sub-band [ZA, ZA+ZN), for ink-scale resolution the winding LOD doesn't carry.
"""
import sys, os, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
TILER = os.path.join(HERE, "tile_winding_mr.py")
UNROLL = os.path.join(HERE, "..", "build", "tools", "unroll_wind")

def main():
    a = sys.argv[1:]
    # split off optional flags
    samp = "160"; fine = None
    if "--samp" in a:
        i = a.index("--samp"); samp = a[i+1]; del a[i:i+2]
    if "--fine" in a:
        i = a.index("--fine"); fine = a[i+1:i+5]; del a[i:i+5]   # FINEARC SCALE ZA ZN
    # positional tiling args (12 required; NZ STEPZ optional)
    if len(a) < 12:
        print(__doc__); sys.exit(2)
    arc, flod, outdir = a[0], a[1], a[11]
    os.makedirs(outdir, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "8"),
               ASAN_OPTIONS="detect_leaks=0")

    print("=== [1/2] tiling -> merged winding ===", flush=True)
    rc = subprocess.run([sys.executable, TILER] + a, env=env).returncode
    merged = os.path.join(outdir, "merged.f32")
    if rc != 0 or not os.path.exists(merged):
        sys.exit(f"tiling failed (rc={rc})")

    print("\n=== [2/2] unroll (auto-center spiral) ===", flush=True)
    out = os.path.join(outdir, "unroll.pgm")
    cmd = [UNROLL, arc, flod, merged, out, samp, "auto"]
    if fine:
        cmd += fine                                   # FINEARC SCALE ZA ZN -> FINE mode
        out = os.path.join(outdir, "unroll_fine.pgm"); cmd[3] = out
    rc = subprocess.run(cmd, env=env).returncode
    if rc != 0:
        sys.exit(f"unroll failed (rc={rc})")
    print(f"\nDONE -> {out}")

if __name__ == "__main__":
    main()
