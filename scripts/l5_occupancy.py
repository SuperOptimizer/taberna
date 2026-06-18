#!/usr/bin/env python3
"""Estimate LOD0 occupancy of a scroll cheaply: download level 5/ (each voxel = a
32^3 LOD0 block), count non-zero voxels, multiply by 32^3. Robust to S3 connection
resets (retries + low concurrency)."""
import json, sys, math, time, urllib.request, urllib.error
import numpy as np
from concurrent.futures import ThreadPoolExecutor

BUCKET = "https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com"
ZARR = "PHerc0332/volumes/20251211183505-2.399um-0.2m-78keV-masked.zarr"
L = 5
DENSE = 33592 * 15761 * 15761      # LOD0 dense voxels

def get(url, tries=12):
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if a == tries - 1:
                raise
        except Exception:
            if a == tries - 1:
                raise
        time.sleep(min(1.5 * (a + 1), 12))
    return None

za = json.loads(get(f"{BUCKET}/{ZARR}/{L}/.zarray").decode())
SZ, SY, SX = za["shape"]; CZ, CY, CX = za["chunks"]
ncz, ncy, ncx = math.ceil(SZ/CZ), math.ceil(SY/CY), math.ceil(SX/CX)
sys.stderr.write(f"L{L} shape {SZ}x{SY}x{SX} chunks {CZ}x{CY}x{CX} grid {ncz}x{ncy}x{ncx} ({ncz*ncy*ncx} chunks)\n")

def fetch(c):
    cz, cy, cx = c
    b = get(f"{BUCKET}/{ZARR}/{L}/{cz}/{cy}/{cx}")
    if b is None:
        return (0, 0)
    a = np.frombuffer(b, np.uint8)
    return (int(np.count_nonzero(a)), 1)

coords = [(cz, cy, cx) for cz in range(ncz) for cy in range(ncy) for cx in range(ncx)]
nonzero = present = done = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    for nz, pr in ex.map(fetch, coords):
        nonzero += nz; present += pr; done += 1
        if done % 25 == 0:
            sys.stderr.write(f"  {done}/{len(coords)} chunks, nonzero so far {nonzero:,}\n")

total_l5 = SZ * SY * SX
mat = nonzero * (32**3)            # LOD0 material voxels = bytes (u8)
print(f"L5 non-zero voxels : {nonzero:,} / {total_l5:,}  ({100*nonzero/total_l5:.1f}% of L5)")
print(f"present L5 chunks   : {present}/{len(coords)}")
print(f"LOD0 material (x32^3): {mat:,} voxels = {mat/1e12:.2f} TB (u8 -> same in bytes)")
print(f"LOD0 dense          : {DENSE/1e12:.2f} TB")
print(f"==> occupancy       : {100*mat/DENSE:.1f}%")
print(f"est. .mca @ ~25x LOD0: {mat/25/1e9:.0f} GB")
