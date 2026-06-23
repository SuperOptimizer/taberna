#!/usr/bin/env python3
"""End-to-end classical unroll pipeline: region -> coarse winding -> fine tiles -> merge -> unroll.

One reproducible entry point chaining the validated pieces:
  1. scripts/tile_winding_mr.py : coarse global winding + GLAM-anchored fine tiles -> OUTDIR/merged.f32
  1.5 tools/wind_tv (global weighted-TV regularizer) : merged.f32 -> OUTDIR/merged_tv.f32. ON by default
     (--no-tv to skip). Validated: the TV-cleaned field unrolls measurably sharper (+73% sheet-edge
     sharpness, +33% contrast) because the smoother field bins each sheet into a tighter render-coord
     range. See docs/touching-sheets-plan.md Phase 3a.
  2. tools/unroll_wind (auto umbilicus center) : (cleaned) winding -> OUTDIR/unroll.pgm (spiral flatten)
     optionally FINE mode: sample the winding against a finer archive for native-res output.

Usage:
  unroll_pipeline.py ARC FINELOD Z0 Y0 X0 DZ TS NY NX STEP COARSELOD OUTDIR
                     [NZ STEPZ] [--samp N] [--fine FINEARC SCALE ZA ZN] [--no-tv] [--tv-lambda L]
All region args are the tile_winding_mr convention (FINELOD-grid).
  --fine: render the (L2/L1) winding against FINEARC (e.g. the L1 or L0 archive) at SCALE x, over the
          winding-z sub-band [ZA, ZA+ZN), for ink-scale resolution the winding LOD doesn't carry.
  --no-tv: unroll the RAW merged field (skip the wind_tv regularization).
"""
import sys, os, subprocess, struct

HERE = os.path.dirname(os.path.abspath(__file__))
TILER = os.path.join(HERE, "tile_winding_mr.py")
UNROLL = os.path.join(HERE, "..", "build", "tools", "unroll_wind")
WIND_TV = os.path.join(HERE, "..", "build", "tools", "wind_tv")

def main():
    a = sys.argv[1:]
    # split off optional flags
    samp = "160"; fine = None; do_tv = True; tv_lambda = "0.3"
    if "--samp" in a:
        i = a.index("--samp"); samp = a[i+1]; del a[i:i+2]
    if "--fine" in a:
        i = a.index("--fine"); fine = a[i+1:i+5]; del a[i:i+5]   # FINEARC SCALE ZA ZN
    if "--no-tv" in a:
        do_tv = False; a.remove("--no-tv")
    if "--tv-lambda" in a:
        i = a.index("--tv-lambda"); tv_lambda = a[i+1]; del a[i:i+2]
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

    field = merged
    if do_tv:
        print("\n=== [1.5] wind_tv regularize (cleaner field -> sharper unroll) ===", flush=True)
        # merged.f32 header = int32 {dz,dy,dx,z0,y0,x0} on the FINELOD grid; wind_tv reads the matching
        # ARC region for sheetness and the prior winding from merged.f32 (same grid -> priorlod=FINELOD).
        with open(merged, "rb") as fh:
            dz_, dy_, dx_, mz0, my0, mx0 = struct.unpack("<6i", fh.read(24))
        base = os.path.join(outdir, "mtv")
        cmd = [WIND_TV, arc, base, flod, str(mz0), str(my0), str(mx0),
               str(dz_), str(dy_), str(dx_), merged, flod, tv_lambda]
        rc = subprocess.run(cmd, env=env).returncode
        tvfield = base + "_tv.f32"
        if rc != 0 or not os.path.exists(tvfield):
            print(f"*** wind_tv failed (rc={rc}); falling back to the RAW merged field", flush=True)
        else:
            field = tvfield

    print("\n=== [2/2] unroll (auto-center spiral) ===", flush=True)
    out = os.path.join(outdir, "unroll.pgm")
    cmd = [UNROLL, arc, flod, field, out, samp, "auto"]
    if fine:
        cmd += fine                                   # FINEARC SCALE ZA ZN -> FINE mode
        out = os.path.join(outdir, "unroll_fine.pgm"); cmd[3] = out
    rc = subprocess.run(cmd, env=env).returncode
    if rc != 0:
        sys.exit(f"unroll failed (rc={rc})")
    print(f"\nDONE -> {out}")

if __name__ == "__main__":
    main()
