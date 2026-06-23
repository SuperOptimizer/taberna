#!/usr/bin/env python3
"""Parameter sensitivity sweep (P3): vary one solver param at a time around its default at a fixed
centered location, measure backward-switch / fill-climb / wrap-count, and report which params
actually move the result (fragile = lock tightly) vs which are inert (safe to leave).

Usage: param_sweep.py [ARC] [z0 y0 x0 dz dy dx]
"""
import sys, os, re, subprocess
ARC = sys.argv[1] if len(sys.argv) > 1 else "data/exports/pherc0332_L2.mca"
Z0,Y0,X0,DZ,DY,DX = (sys.argv[2:8] + ["3936","718","1580","32","1024","1024"][len(sys.argv)-2:]) if len(sys.argv)>2 else ("3936","718","1580","32","1024","1024")
TOOL = os.path.join(os.path.dirname(__file__), "..", "build", "tools", "sheet_sep3d")

# positional tail after z0 y0 x0 dz dy dx: minseg pitch LAM WEQ ztie recto LPRI radw barclose zmed SP
#   umbref wdr robsig usesheet slicecv cxoff [priorvol priorlod GLAM NF]
# defaults (standalone, auto-pitch):
D = dict(minseg="40", pitch="0", LAM="0.6", WEQ="3", ztie="0", recto="0", LPRI="0.15", radw="1.0",
         barclose="1", zmed="0", SP="0", umbref="1", wdr="0", robsig="0.6", usesheet="1",
         slicecv="0", cxoff="0", priorvol="", priorlod="", GLAM="0", NF="0")
ORDER = ["minseg","pitch","LAM","WEQ","ztie","recto","LPRI","radw","barclose","zmed","SP","umbref",
         "wdr","robsig","usesheet","slicecv","cxoff","priorvol","priorlod","GLAM","NF"]
# param -> values to sweep
SWEEP = {"LAM":["0.3","0.6","1.0"], "WEQ":["1","3","6"], "LPRI":["0.1","0.15","0.4"],
         "radw":["0.6","1.0"], "robsig":["0.4","0.6","0.9"], "NF":["0","0.4","0.7"],
         "barclose":["0","1"]}
RX = {"backsw": re.compile(r"backward-switch=([\d.]+)/1000"),
      "climb":  re.compile(r"climb=([\d.-]+) wraps"),
      "wraps":  re.compile(r"GLOBAL winding:.*\(([\d.]+) wraps\)"),
      "align":  re.compile(r"grad\(w\)-normal alignment=([\d.]+)")}

def run(overrides):
    d = dict(D); d.update(overrides)
    args = [TOOL, ARC, "/tmp/psw", "0", Z0,Y0,X0,DZ,DY,DX] + [d[k] for k in ORDER]
    env = dict(os.environ, OMP_NUM_THREADS="8", ASAN_OPTIONS="detect_leaks=0")
    p = subprocess.run(args, capture_output=True, text=True, env=env, timeout=600)
    out = p.stdout + "\n" + p.stderr; r = {}
    for k,rx in RX.items():
        m = rx.search(out); r[k] = float(m.group(1)) if m else None
    return r

def main():
    base = run({}); print(f"location z{Z0} y{Y0} x{X0} {DZ}x{DY}x{DX}")
    print(f"BASELINE: backsw={base['backsw']} climb={base['climb']} wraps={base['wraps']} align={base['align']}\n")
    print(f"{'param':8} {'value':6} {'backsw':>7} {'climb':>7} {'wraps':>6} {'align':>6}  sensitivity")
    print("-"*60)
    for param, vals in SWEEP.items():
        rows = []
        for v in vals:
            r = run({param: v}); rows.append((v, r))
        bs = [r['backsw'] for _,r in rows if r['backsw'] is not None]
        spread = (max(bs)-min(bs)) if len(bs)>1 else 0
        tag = "FRAGILE" if spread>25 else ("mild" if spread>10 else "inert")
        for v,r in rows:
            print(f"{param:8} {v:6} {str(r['backsw']):>7} {str(r['climb']):>7} {str(r['wraps']):>6} {str(r['align']):>6}")
        print(f"  -> backsw spread {spread:.0f} = {tag}\n")

if __name__ == "__main__":
    main()
