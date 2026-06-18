#!/usr/bin/env python3
"""Recursively list a zarr level prefix on the (public) open-data bucket and report
present-chunk count + total bytes => actual occupancy of that level.

Usage: s3_occupancy.py <prefix ending in /0/>   (anonymous; paginated ListObjectsV2)
"""
import sys, urllib.request, urllib.parse
import xml.etree.ElementTree as ET

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
NS = "{http://s3.amazonaws.com/doc/2006-03-01/}"
prefix = sys.argv[1] if len(sys.argv) > 1 else \
    "PHerc0332/volumes/20251211183505-2.399um-0.2m-78keV-masked.zarr/0/"

n_chunks = 0          # data chunks (exclude .zarray etc.)
total_bytes = 0
token = None
pages = 0
while True:
    url = f"{BUCKET}/?list-type=2&prefix={urllib.parse.quote(prefix)}&max-keys=1000"
    if token:
        url += f"&continuation-token={urllib.parse.quote(token)}"
    import time
    data = None
    for attempt in range(8):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = r.read(); break
        except Exception:
            if attempt == 7: raise
            time.sleep(1.0 * (attempt + 1))
    root = ET.fromstring(data)
    for c in root.findall(f"{NS}Contents"):
        key = c.find(f"{NS}Key").text
        size = int(c.find(f"{NS}Size").text)
        if key.endswith(".zarray") or key.endswith(".zattrs"):
            continue
        n_chunks += 1
        total_bytes += size
    pages += 1
    if pages % 25 == 0:
        sys.stderr.write(f"  ...{pages} pages, {n_chunks} chunks, {total_bytes/1e9:.1f} GB\n")
    tok = root.find(f"{NS}NextContinuationToken")
    token = tok.text if (tok is not None and tok.text) else None
    if not token:
        break

# total possible chunks for the dense bounding box (from the known LOD0 shape)
SZ, SY, SX, C = 33592, 15761, 15761, 128
import math
poss = math.ceil(SZ/C)*math.ceil(SY/C)*math.ceil(SX/C)
dense_bytes = SZ*SY*SX
print(f"prefix: {prefix}")
print(f"present chunks : {n_chunks:,}")
print(f"possible chunks: {poss:,}   -> occupancy {100*n_chunks/poss:.1f}%")
print(f"actual bytes   : {total_bytes/1e9:.1f} GB  ({total_bytes/1e12:.2f} TB)")
print(f"count x 128^3  : {n_chunks*C*C*C/1e9:.1f} GB  (uniform-chunk estimate)")
print(f"dense bounding : {dense_bytes/1e12:.2f} TB   -> stored is {100*total_bytes/dense_bytes:.1f}% of dense")
