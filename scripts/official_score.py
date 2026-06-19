#!/usr/bin/env python3
"""Run the OFFICIAL Kaggle Surface-Detection metric on a taberna prediction.

This is the authoritative scorer: it calls the competition's `topometrics`
library (Betti-Matching-3D TopoScore + Google surface_distance NSD + cc3d/skimage
VOI), with the exact leaderboard parameters. Use it to validate taberna's fast
in-process proxy (`src/eval/score.c`).

Setup (one time) — see scripts/metric_setup.md. In short, from the team archive:
    cd topological-metrics-kaggle
    uv pip install --no-index --find-links=../wheels -r requirements.txt
    make build-betti
    uv pip install -e . --no-deps --no-index --no-build-isolation

Usage:
    python scripts/official_score.py GT.tif PRED.tif
    # PRED.tif is produced by:  tools/surface_predict IMAGE.tif GT.tif PRED.tif
"""
import sys
import numpy as np
from PIL import Image, ImageSequence


def load_volume(path):
    im = Image.open(path)
    return np.stack([np.array(p) for p in ImageSequence.Iterator(im)], axis=0)


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    gt = load_volume(sys.argv[1])
    pr = load_volume(sys.argv[2])

    import topometrics.leaderboard as lb

    rep = lb.compute_leaderboard_score(
        predictions=pr,
        labels=gt,
        dims=(0, 1, 2),
        spacing=(1.0, 1.0, 1.0),
        surface_tolerance=2.0,
        voi_connectivity=26,
        voi_transform="one_over_one_plus",
        voi_alpha=0.3,
        combine_weights=(0.30, 0.35, 0.35),  # (Topo, SurfaceDice, VOI)
        fg_threshold=None,                    # legacy "!= 0"
        ignore_label=2,
    )
    score = float(np.clip(rep.score, 0.0, 1.0))
    print(f"OFFICIAL score : {score:.4f}")
    print(f"  TopoScore    : {float(rep.topo.toposcore):.4f}  "
          f"(F1 by dim {rep.topo.topoF1_by_dim})")
    print(f"  SurfaceDice  : {float(rep.surface_dice):.4f}")
    print(f"  VOI_score    : {float(rep.voi.voi_score):.4f}  "
          f"(VOI={rep.voi.voi_total:.4f} split={rep.voi.voi_split:.4f} "
          f"merge={rep.voi.voi_merge:.4f})")


if __name__ == "__main__":
    main()
