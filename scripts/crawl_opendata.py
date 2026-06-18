#!/usr/bin/env python3
"""Crawl the Vesuvius Challenge open-data S3 bucket (anonymous) and build a JSON
manifest of every scroll + volume: LOD0 dims from 0/.zarray, plus the
metadata.json / transform.json for volumes that have them (saved locally).

Hierarchy: <scroll>/volumes/<name>.zarr/{0..N/, metadata.json, transform.json}
Lists prefixes only (delimiter=/) so it never touches the millions of chunk files.
"""
import json, re, sys, time, urllib.request, urllib.error
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUCKET = "vesuvius-challenge-open-data"
REGION = "us-east-1"
BASE = f"https://{BUCKET}.s3.{REGION}.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"
OUT = Path(__file__).resolve().parent.parent / "data" / "opendata"

# VC scroll-name mapping for the well-known fragments (best-effort; others null).
SCROLL_NAMES = {
    "PHercParis4": "Scroll 1", "PHercParis3": "Scroll 2", "PHerc0332": "Scroll 3",
    "PHerc1667": "Scroll 4", "PHerc0172": "Scroll 5",
}

def _get(url, binary=False):
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            b = r.read()
            return b if binary else b.decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return None
        raise
    except Exception as e:
        sys.stderr.write(f"GET fail {url}: {e}\n")
        return None

def list_prefix(prefix, want="prefixes"):
    """List one level under `prefix` (delimiter=/). Returns CommonPrefixes (subdirs)
    or Contents keys depending on `want`. Handles continuation tokens."""
    out, token = [], None
    while True:
        url = f"{BASE}/?list-type=2&prefix={urllib.parse.quote(prefix)}&delimiter=/"
        if token:
            url += f"&continuation-token={urllib.parse.quote(token)}"
        xml = _get(url)
        if xml is None:
            break
        root = ET.fromstring(xml)
        if want == "prefixes":
            out += [p.find(f"{NS}Prefix").text for p in root.findall(f"{NS}CommonPrefixes")]
        else:
            out += [c.find(f"{NS}Key").text for c in root.findall(f"{NS}Contents")]
        tok = root.find(f"{NS}NextContinuationToken")
        if tok is not None and tok.text:
            token = tok.text
        else:
            break
    return out

def parse_name(name):
    info = {}
    m = re.search(r"(\d+\.?\d*)um", name);     info["voxel_um"] = float(m.group(1)) if m else None
    m = re.search(r"(\d+)keV", name);          info["energy_kev"] = int(m.group(1)) if m else None
    m = re.search(r"-(\d+\.?\d*)m-", name);     info["distance_m"] = float(m.group(1)) if m else None
    info["masked"] = "masked" in name
    return info

def zarray_dims(blob):
    if blob is None:
        return None
    try:
        z = json.loads(blob)
        return {"shape": z.get("shape"), "chunks": z.get("chunks"),
                "dtype": z.get("dtype"), "compressor": z.get("compressor"),
                "dimension_separator": z.get("dimension_separator"),
                "zarr_format": z.get("zarr_format")}
    except Exception:
        return None

def process_volume(scroll, vol_prefix):
    # vol_prefix e.g. "PHerc0332/volumes/2025...-masked.zarr/"
    name = vol_prefix.rstrip("/").split("/")[-1]
    base = f"{BASE}/{urllib.parse.quote(vol_prefix)}"
    # one list to discover levels + which json siblings exist
    entries = list_prefix(vol_prefix, "prefixes")           # 0/ 1/ ...
    keys = list_prefix(vol_prefix, "keys")                  # .zattrs, metadata.json, ...
    levels = sorted(int(p.rstrip("/").split("/")[-1]) for p in entries
                    if p.rstrip("/").split("/")[-1].isdigit())
    keyset = {k.split("/")[-1] for k in keys}
    # LOD0 dims: <vol>/0/.zarray, else flat <vol>/.zarray
    z0 = zarray_dims(_get(base + "0/.zarray")) if 0 in levels else None
    if z0 is None:
        z0 = zarray_dims(_get(base + ".zarray"))
    rec = {"name": name, "zarr": vol_prefix.rstrip("/"), "is_zarr": name.endswith(".zarr"),
           "levels": levels, "lod0": z0, **parse_name(name),
           "has_metadata": "metadata.json" in keyset, "has_transform": "transform.json" in keyset}
    # save metadata.json / transform.json locally for those that have them
    vdir = OUT / scroll / name
    for fn in ("metadata.json", "transform.json"):
        if fn in keyset:
            blob = _get(base + fn, binary=True)
            if blob is not None:
                vdir.mkdir(parents=True, exist_ok=True)
                (vdir / fn).write_bytes(blob)
    return rec

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    scrolls = [p.rstrip("/") for p in list_prefix("", "prefixes")]
    scrolls = [s for s in scrolls if not s.startswith("_")]
    sys.stderr.write(f"{len(scrolls)} scrolls\n")
    # collect (scroll, vol_prefix) pairs
    tasks = []
    for s in scrolls:
        vols = list_prefix(f"{s}/volumes/", "prefixes")
        for v in vols:
            tasks.append((s, v))
    sys.stderr.write(f"{len(tasks)} volumes; fetching dims + metadata...\n")
    manifest = {"bucket": BUCKET, "region": REGION,
                "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "scrolls": {}}
    results = {}
    with ThreadPoolExecutor(max_workers=24) as ex:
        futs = {ex.submit(process_volume, s, v): (s, v) for (s, v) in tasks}
        done = 0
        for f in futs:
            pass
        for f in list(futs):
            s, v = futs[f]
            try:
                rec = f.result()
            except Exception as e:
                rec = {"name": v, "error": str(e)}
            results.setdefault(s, []).append(rec)
            done += 1
            if done % 20 == 0:
                sys.stderr.write(f"  {done}/{len(tasks)}\n")
    nvol = 0
    for s in sorted(scrolls):
        vols = sorted(results.get(s, []), key=lambda r: r.get("name", ""))
        nvol += len(vols)
        manifest["scrolls"][s] = {"scroll_name": SCROLL_NAMES.get(s), "n_volumes": len(vols),
                                  "volumes": vols}
    manifest["scroll_count"] = len(scrolls)
    manifest["volume_count"] = nvol
    mpath = OUT / "manifest.json"
    mpath.write_text(json.dumps(manifest, indent=2))
    sys.stderr.write(f"wrote {mpath} ({nvol} volumes across {len(scrolls)} scrolls)\n")
    print(str(mpath))

if __name__ == "__main__":
    import urllib.parse
    main()
