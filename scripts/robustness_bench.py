#!/usr/bin/env python3
"""Multi-location robustness benchmark for the sheet_sep3d winding solver.

Runs the STANDALONE solver (auto-pitch, no coarse prior) at a set of diverse L2 locations and
parses its diagnostic suite into one table, so we can see WHERE the pipeline is robust and where it
breaks -- instead of validating at a single tuned spot. No ground truth: we report the
ground-truth-free checks (scale gate, edge consistency, crossing-count agreement, fill climb,
backward-switch, z-coherence) and FLAG the locations that fail any of them.

Usage: robustness_bench.py [ARC] [OUTDIR] [JOBS]
  Locations are defined in LOCATIONS below (edit to add coverage). Each is a standalone crop.

NOTE: this measures STANDALONE (no coarse prior) robustness = the worst case. In production every
fine tile gets a coarse prior + GLAM, which RESCUES the standalone-weak regimes: off-center tile
backward-switch 134->62 (halved), inner-core scale graph 20->14 (snaps to the direct sheet count).
So a standalone FLAG on a partial/off-center/core crop is NOT a production failure -- the centered
coarse solve (scale 0-10% at all centered locations here) sets the global reference and the prior
propagates it. Trust the centered-crop rows for the global scale; treat partial-crop flags as
"needs the prior" not "broken".
"""
import sys, os, re, subprocess, csv
from concurrent.futures import ThreadPoolExecutor

ARC    = sys.argv[1] if len(sys.argv) > 1 else "data/exports/pherc0332_L2.mca"
OUTDIR = sys.argv[2] if len(sys.argv) > 2 else "/home/forrest/mcdebug/bench"
JOBS   = int(sys.argv[3]) if len(sys.argv) > 3 else 3
TOOL   = os.path.join(os.path.dirname(__file__), "..", "build", "tools", "sheet_sep3d")
os.makedirs(OUTDIR, exist_ok=True)

# (name, lod, z0, y0, x0, dz, dy, dx) -- in the archive's lod-grid. Diverse by z, position, size.
# L2 dims: nz=7936 ny=2560 nx=3584. The y718 x1580 1024-crop contains the umbilicus at z~3900.
LOCATIONS = [
    # --- centered crops (contain umbilicus) at different z: tests z-stability of the scale ---
    ("z2500_c", 0, 2500, 718, 1580, 32, 1024, 1024),
    ("z3200_c", 0, 3200, 718, 1580, 32, 1024, 1024),
    ("z3936_c", 0, 3936, 718, 1580, 32, 1024, 1024),
    ("z4600_c", 0, 4600, 718, 1580, 32, 1024, 1024),
    ("z5300_c", 0, 5300, 718, 1580, 32, 1024, 1024),
    # --- bigger centered crop: full cross-section span (more wraps, outer delaminated) ---
    ("z3936_big", 0, 3936, 350, 1100, 32, 1792, 1792),
    # --- off-center / partial tiles (umbilicus near edge or outside): tests partial regime ---
    ("z3936_offA", 0, 3936, 718, 1580, 32, 640, 640),
    ("z3936_offB", 0, 3936, 1300, 2200, 32, 640, 640),
    # --- inner compressed core only (small radius band near center) ---
    ("z3936_core", 0, 3936, 980, 1840, 32, 512, 512),
    # --- thicker z-band: tests z-coherence over more slices ---
    ("z3936_tallz", 0, 3872, 718, 1580, 96, 1024, 1024),
    # --- different angular sector of the scroll ---
    ("z3936_quad2", 0, 3936, 350, 1580, 32, 1024, 1024),
    ("z3936_quad3", 0, 3936, 718, 2092, 32, 1024, 1024),
]

# regex -> (field name) for the lines sheet_sep3d prints
PATTERNS = {
    "pitch":      re.compile(r"^pitch=([\d.]+)"),
    "pitch_cal":  re.compile(r"pitch calibrated.*-> pitch=([\d.]+)"),
    "edge_in_pct":re.compile(r"edge consistency: \d+/\d+ \(([\d.]+)%\)"),
    "wraps":      re.compile(r"GLOBAL winding:.*\(([\d.]+) wraps\), inlier-resid=([\d.]+)"),
    "crossing":   re.compile(r"radial-crossing estimate.*raw-valley=([\d.]+), SHEETNESS=([\d.]+), graph-solved=([\d.]+)"),
    "crossing2":  re.compile(r"radial-crossing estimate.*raw-valley=([\d.]+), graph-solved=([\d.]+)"),
    "scale_ok":   re.compile(r"scale check OK: graph ([\d.]+) vs crossing ([\d.]+) \(([\d.]+)%\)"),
    "scale_warn": re.compile(r"SCALE WARNING: graph-solved ([\d.]+) wraps vs crossing-count ([\d.]+) \(([\d.]+)%"),
    "scale_skip": re.compile(r"scale gate skipped \((.*?) --"),
    "backsw":     re.compile(r"backward-switch=([\d.]+)/1000"),
    "climb":      re.compile(r"fill coverage:.*climb=([\d.-]+) wraps.*frac=([\d.]+)"),
    "zcoh":       re.compile(r"z-coherence: std across z = ([\d.]+) wraps \(raw\).*adj-z jump ([\d.]+)"),
}

def run_one(loc):
    name, lod, z0, y0, x0, dz, dy, dx = loc
    out = f"{OUTDIR}/{name}"
    args = [TOOL, ARC, out, str(lod), str(z0), str(y0), str(x0), str(dz), str(dy), str(dx), "40", "0"]
    env = dict(os.environ, OMP_NUM_THREADS="6", ASAN_OPTIONS="detect_leaks=0")
    try:
        p = subprocess.run(args, capture_output=True, text=True, env=env, timeout=900)
    except subprocess.TimeoutExpired:
        return {"name": name, "status": "TIMEOUT"}
    r = {"name": name, "z": z0, "size": f"{dz}x{dy}x{dx}", "status": "ok"}
    # sheet_sep3d prints some diagnostics to stdout (printf) and some to stderr (fprintf) -- parse both
    for line in (p.stdout + "\n" + p.stderr).splitlines():
        for key, rx in PATTERNS.items():
            m = rx.search(line)
            if not m:
                continue
            if key == "pitch":      r["pitch"] = float(m.group(1))
            elif key == "pitch_cal":r["pitch_cal"] = float(m.group(1))
            elif key == "edge_in_pct": r["edge_in%"] = float(m.group(1))
            elif key == "wraps":    r["wraps"] = float(m.group(1)); r["resid"] = float(m.group(2))
            elif key == "crossing": r["raw_valley"], r["sheetness"], r["graph"] = map(float, m.groups())
            elif key == "crossing2" and "sheetness" not in r:
                r["raw_valley"], r["graph"] = float(m.group(1)), float(m.group(2))
            elif key == "scale_ok":   r["scale"] = f"OK {m.group(3)}%"; r["scale_off%"] = float(m.group(3))
            elif key == "scale_warn": r["scale"] = f"WARN {m.group(3)}%"; r["scale_off%"] = float(m.group(3))
            elif key == "scale_skip": r["scale"] = f"skip({m.group(1)[:14]})"
            elif key == "backsw":   r["backsw"] = float(m.group(1))
            elif key == "climb":    r["climb"] = float(m.group(1)); r["unfilled"] = float(m.group(2))
            elif key == "zcoh":     r["zstd"] = float(m.group(1)); r["zjump"] = float(m.group(2))
    if "wraps" not in r:
        r["status"] = "SOLVE-FAIL"
        r["stderr_tail"] = " | ".join(p.stderr.strip().splitlines()[-3:])
    return r

def flags(r):
    """ground-truth-free failure flags."""
    f = []
    if r.get("status") != "ok": f.append(r.get("status"))
    if r.get("scale", "").startswith("WARN"): f.append("SCALE")
    if r.get("edge_in%", 0) > 10: f.append("EDGE")          # off-center / umbilicus bias
    if r.get("backsw", 0) > 120: f.append("BACKSW")          # leaky fill
    if r.get("climb", 99) < 3: f.append("FLAT")              # collapsed fill
    if r.get("unfilled", 0) > 0.30: f.append("UNFILLED")
    if r.get("zstd", 0) > 0.5: f.append("ZINCOH")            # z-incoherent
    if r.get("resid", 0) > 0.20: f.append("RESID")
    return ",".join(f) if f else "-"

def main():
    print(f"benchmark: {len(LOCATIONS)} locations, {JOBS} parallel, arc={ARC}")
    with ThreadPoolExecutor(max_workers=JOBS) as ex:
        results = list(ex.map(run_one, LOCATIONS))
    cols = ["name","z","size","pitch_cal","wraps","raw_valley","sheetness","graph","scale",
            "edge_in%","backsw","climb","unfilled","zstd","zjump","resid","status"]
    # print table
    w = {c: max(len(c), *(len(f"{r.get(c,'')}") for r in results)) for c in cols}
    hdr = "  ".join(c.ljust(w[c]) for c in cols) + "  FLAGS"
    print("\n" + hdr); print("-"*len(hdr))
    nfail = 0
    for r in results:
        fl = flags(r)
        if fl != "-": nfail += 1
        print("  ".join(f"{r.get(c,'')}".ljust(w[c]) for c in cols) + "  " + fl)
    print(f"\n{nfail}/{len(results)} locations flagged. CSV -> {OUTDIR}/bench.csv")
    with open(f"{OUTDIR}/bench.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols+["flags"]); wr.writeheader()
        for r in results:
            row = {c: r.get(c, "") for c in cols}; row["flags"] = flags(r); wr.writerow(row)

if __name__ == "__main__":
    main()
