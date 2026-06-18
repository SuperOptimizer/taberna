#!/usr/bin/env python3
"""Download a chunk-aligned 1024^3 region from the middle of an open-data scroll
volume (LOD0) and write it as a local raw-u8 zarr (0/.zarray + 0/z/y/x chunks,
compressor null) that fysics-process can consume. Copies metadata.json alongside.

Usage: tune_region.py <scroll> [size]   (default scroll PHerc0332, size 1024)
"""
import json, sys, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
scroll = sys.argv[1] if len(sys.argv) > 1 else "PHerc0332"
SIZE = int(sys.argv[2]) if len(sys.argv) > 2 else 1024

man = json.load(open(ROOT / "data/opendata/manifest.json"))
vols = man["scrolls"][scroll]["volumes"]
# pick the ~2.4um masked volume (the ESRF reconstruction)
vol = next((v for v in vols if v.get("voxel_um") and 2.3 <= v["voxel_um"] <= 2.5), vols[0])
zarr_key = vol["zarr"]                     # e.g. PHerc0332/volumes/...zarr
z0 = vol["lod0"]
SZ, SY, SX = z0["shape"]
CZ, CY, CX = z0["chunks"]
assert CZ == CY == CX == 128, f"unexpected chunk shape {z0['chunks']}"
assert SIZE % 128 == 0
nchunk = SIZE // 128

# chunk-aligned region centered in the volume
def start(dim):
    c = dim // 2
    s = (c // 128 - nchunk // 2) * 128
    return max(0, min(s, dim - SIZE))
oz, oy, ox = start(SZ), start(SY), start(SX)
print(f"{scroll}  vol={vol['name']}\n  volume {SZ}x{SY}x{SX}  region {SIZE}^3 at z{oz} y{oy} x{ox}")

outdir = ROOT / "data/tune" / scroll / "region.zarr"
(outdir / "0").mkdir(parents=True, exist_ok=True)
(outdir / "0" / ".zarray").write_text(json.dumps({
    "chunks": [128, 128, 128], "compressor": None, "dtype": "|u1", "fill_value": 0,
    "order": "C", "shape": [SIZE, SIZE, SIZE], "zarr_format": 2, "dimension_separator": "/"}))
# copy metadata.json (downloaded by the crawler) so fysics-process reads physics params
src_meta = ROOT / "data/opendata" / scroll / vol["name"] / "metadata.json"
if src_meta.exists():
    (outdir / "metadata.json").write_bytes(src_meta.read_bytes())

base = f"{BUCKET}/{zarr_key}/0"
import time
def fetch(cz, cy, cx):
    az, ay, ax = oz//128 + cz, oy//128 + cy, ox//128 + cx
    url = f"{base}/{az}/{ay}/{ax}"
    data = None
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read()
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return (cz, cy, cx, 0)    # absent -> fill 0 (omit file)
            if attempt == 5:
                raise
        except Exception:
            if attempt == 5:
                raise
        time.sleep(0.5 * (attempt + 1))   # backoff on reset/throttle
    d = outdir / "0" / str(cz) / str(cy)
    d.mkdir(parents=True, exist_ok=True)
    (d / str(cx)).write_bytes(data)
    return (cz, cy, cx, len(data))

coords = [(z, y, x) for z in range(nchunk) for y in range(nchunk) for x in range(nchunk)]
present = 0
with ThreadPoolExecutor(max_workers=12) as ex:
    for (cz, cy, cx, n) in ex.map(lambda c: fetch(*c), coords):
        if n:
            present += 1
print(f"  chunks: {present}/{len(coords)} present ({present*2}MB)  ->  {outdir}")
print(str(outdir))
