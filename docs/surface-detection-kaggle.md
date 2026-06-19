# Kaggle "Vesuvius Challenge – Surface Detection" — solution analysis (for taberna)

Competition ran ~Dec 2025 → Feb 2026, $100K, hosted by Scroll Prize
(@giorgioangelotti, @seanjohnsonsp). Writeups dated Feb 27–28 2026. This doc
analyzes the top solutions through taberna's lens: **what is transferable to a
no-ML classical pipeline, and where taberna can actually differentiate.**

## The task

- **Input:** chunks of 3D micro-CT scroll volume.
- **Label:** binary mask of "smoothed sheet positions" (papyrus surfaces) winding
  through the volume. Label value 2 = ignore/uncertain voxels (notably a 3-voxel
  border on each volume face — a recurring nuisance everyone special-cased).
- **Output:** per-voxel binary surface mask.
- **Metric:** topology-aware linear blend = **surface proximity (Surface Dice)**
  + **instance consistency (VOI split/merge)** + **shape/topology correctness
  (TopoScore via persistent homology / Betti matching)**. Rewards: no gaps, no
  holes/tunnels (Betti-1), no internal cavities (Betti-2), no splits within a
  wrap, no mergers across adjacent wraps. **Component *location* matters** (it's
  persistence-homology based, not just counts). Top scores ~0.60–0.63.

This metric suite is essentially the one taberna already chose
([[taberna-segmentation]] eval: VOI, ARAND, SurfaceDice@τ, TopoScore/Betti).

## The dominant pattern (the big lesson)

**The model is commoditized; the leaderboard was won in classical
post-processing.** Nearly every top-10 team trained the *same* thing — an
nnU-Net ResEnc-UNet — with minor variations (patch size, epochs, ensemble
weights). Differences in the base model were worth hundredths; the
post-processing was worth the win. All of that post-processing is classical
image-processing / computational-topology — **exactly taberna's territory.**

Common model recipe (interchangeable):
- nnU-Net ResEnc, `fold: all` (train on everything, no val split), batch size 2.
- Long training: 1000 epochs (default) up to **4000** (1st place); "long epochs +
  small batch beats short epochs + large batch." Score keeps creeping up.
- Patch sizes 128–256, often fine-tuned upward (128 → 192 → 256).
- Ensemble 2–4 models by weighted **logit** or softmax-prob fusion + mirror/rotate
  TTA. (1st place note: logit fusion generalized better than prob fusion.)
- Threshold ~0.2 low; larger thresholds (0.35–0.4) generalized better to private LB.

## Post-processing toolkit (the actual competition, all classical)

Ranked roughly by score contribution (per 1st-place ablation: most of the gain is
**dust removal + small-hole plugging**, not the fancy large-hole repair):

0. **Iterated 3×3×3 binary MEDIAN filter** (10th place, from the 18th-place
   solution) — the best value-per-line in the whole competition. A binary 3×3×3
   median *is a 27-neighborhood majority vote* (voxel on iff ≥14/27 neighbors on);
   iterated ~6–10× it acts as a discrete surface-tension/curvature flow that shaves
   thin spurs and bridges (Betti-1 tunnels, false merges) and closes pinholes while
   leaving flat sheet interior intact. **Raises topo score without costing surface
   Dice**: 10th place got PB 0.614→0.624 (and CV +0.0074) from it, beating closing,
   hole-patching and cavity fill. Apply to your *highest-surface-Dice* model.
1. **Dust removal** — drop connected components < 1000–20000 voxels. *Biggest
   single jump* (1st place: .572 → .586).
2. **Small-hole plugging** — fill 1-voxel holes. 1st place uses a **256-entry
   lookup table over the 2×2×2 neighborhood** making each neighborhood 6-connected
   watertight (inspired by skimage `euler_number`'s neighborhood scheme). Cheaper
   and better-for-Dice than dilation. (.586 → .598.)
3. **`binary_closing`** (spherical footprint r=3) + **`binary_fill_holes`** on the
   whole mask for arbitrary cavities.
4. **Large-hole repair via PCA/height-map interpolation** (1st/2nd/4th all do a
   variant — see below). Satisfying but small score impact (few big holes, small
   area fraction).
5. **Topology detection/repair** — Euler number for cheap tunnel detection; Betti
   matching for the real thing (4th place).
6. **Metric hacking** via Betti-2 (4th place) — small, fragile.

### The PCA / height-map sheet-repair primitive (appears 3×, independently)

Per connected component (= one sheet patch):
```
coords = nonzero voxels;  mean-center
U,S,Vt = SVD(coords)            # PCA: tangent1, tangent2, normal = Vt[0],Vt[1],Vt[2]
(u,v) = coords·[t1,t2];  w = coords·normal     # project to best-fit plane + height
w_grid = griddata((u,v) -> w, regular u,v grid, method='linear')  # fill holes
remap w_grid back to 3D, dilate to 3 voxels
```
This is **local sheet flattening** — a micro version of taberna's global unwrap
(a sheet patch → 2D (u,v) parametrization + height). Gotchas they hit and fixed:
- Interpolation overshoots the component boundary → **flood-remove from grid
  edges** (delete interpolated voxels reachable from the border through original
  background) to respect the true sheet outline.
- Interpolation merges nearby sheets → **over-thicken to 5 voxels, remove
  overlap, erode back to 3** → touching sheets separate again.

## The open problem everyone named: touching / adhering sheets

- **1st place:** "we did not find an effective solution to the problem of touching
  sheets. We just relied on nnU-Net to minimize their occurrence."
- **4th place:** RBF sheet fitting gives 0-hole sheets, "However ... can not handle
  the adhering of multi sheets (RBF fitting would fail). I can not figure it out."
- **2nd place ("A postprocessing win") — the partial crack, and why they beat
  everyone on post-proc:** generate masks at **multiple thresholds (0.2 / 0.5 /
  0.6 / 0.7)**. Higher thresholds erode the weak "bridges" where two wraps touch,
  splitting a merged component into separate sheets. Pipeline: interpolate a
  component; score interpolation vs original by Dice + TP-coverage; **low score ⇒
  it's actually multiple merged sheets ⇒ redo on the next higher-threshold mask**,
  which now resolves into separate components, each interpolated individually.
  Accept oversegmentation over merged-with-holes (metric prefers it).
- **4th place** also used **betti_matching (C++)** persistent-homology barcodes on
  20³ chunks to localize Betti-1 tunnels, then PCA-fill locally. Dieter (5th, the
  one from-scratch non-nnU-Net solution) likewise used a sped-up betti-matching for
  tunnel detection + filling, plus a custom loss as his main differentiator.

## The one classical front-end: 7th place (most aligned with taberna)

Used **Coherence-Enhancing Diffusion (CED)** — PDE anisotropic diffusion where the
diffusion tensor is built from the **structure tensor** (Holoborodko derivative
kernels + Gaussian smoothing); diffuses strongly *along* the dominant sheet
orientation while preserving cross-sheet edges. This is taberna's structure-tensor
sheetness front-end ([[taberna-segmentation]] `sheet_tensor.c`) in a different
guise, and it placed 7th — evidence a classical front-end is viable. *(Need full
writeup to extract the segmentation step after CED.)*

## 10th place — diffeomorphic + median filter (notable)

Multi-stage learned pipeline (init seg → stacked ResEnc refine nets → a
**diffeomorphic** stage that predicts a stationary velocity field, scaling-and-
squaring to an `exp(v)` warp, plus an SDF "topology shift" channel to avoid folding
— it calibrates sheet *shape/thickness*, not topology). Interesting for taberna's
*unwrap* side (the SVF + Jacobian-log-barrier anti-fold is the same machinery as
[[taberna-pipeline]]'s deformation field / det-J fold check — but learned). Their
biggest *post*-processing win, though, was the classical median filter above. Also
notable: a host confirmed semi-supervised/extra-unlabeled-data approaches did *not*
beat plain supervised — i.e. data tricks weren't the edge either; post-proc was.

## Data format (downloaded)

`train_images/<id>.tif` + `train_labels/<id>.tif`, paired by `train.csv`
(`id,scroll_id`). Each is a **320³ u8 multi-page TIFF, LZW-compressed**
(Compression=5, stripped, Predictor=none, one strip/page, little-endian). Image =
CT 0–255. Label = {**0** background, **1** surface, **2** ignore/not-evaluated};
**label 2 is ~71% of the volume** (only ~29% is scored), surface (1) is ~4% of the
volume / ~15% of the valid region. The metric must mask out label-2 voxels.

## Implemented in taberna (this work)

- `third-party/tiff` (the vendored submodule) — added **TIFF-LZW decode** +
  **`tiff_read_volume`** (multi-IFD z-stack, raw+LZW). Verified bit-exact vs the
  Python reference on a real sample. *(Submodule edit — must be committed to
  SuperOptimizer/tiff and the pointer bumped for fresh clones.)*
- `src/io/tiff_vol.{h,c}` — `tiff_load_u8` → taberna z-major u8 volume.
- `src/eval/topo.{h,c}` — `cc_label` (6/26-conn), exact `euler_characteristic`
  (cubical V−E+F−C), `betti_numbers` (b0/b1/b2/χ), `region_b1`. Verified on
  ball/cavity/tunnel.
- `src/postproc/morph.{h,c}` — `majority_filter` (iterated median),
  `remove_small_components` (dust), `fill_holes`, `plug_pinholes`, ball
  dilate/erode/closing/opening.
- `src/eval/metrics.c` — `eval_surface` (ignore-aware Dice + Surface Dice@tol +
  prediction Betti, masking label-2).
- `tools/test_topo.c` — invariant unit tests (all pass, incl. majority unmerge
  b0 1→2). `tools/run_surface.c` — end-to-end classical baseline (CT → sheetness →
  threshold → dust/median/fill → ignore-aware score).

## First classical baseline (one 320³ cube, untuned)

`run_surface` with structure-tensor sheetness: thresh 0.5 → 80% of voxels fire
(Dice 0.27); thresh 0.95 → 6.8% but misplaced (Dice 0.11, b0=290). **Finding:**
raw sheetness is a *thick* response over whole sheet bodies and fires across the
(entirely layered) scroll, while the label is the thin *smoothed centerline*. The
front-end needs **ridge non-max-suppression to centerlines along the sheet normal
+ scale tuning (sigma_grad/sigma_tensor) + intensity gating**, not just a
threshold. The plumbing (IO/metric/postproc) is verified; the detector is the open
work. *Pending:* PCA height-map repair, multi-threshold unmerge, 256-entry 2×2×2
watertight LUT, the official composite metric (persistent-homology Betti matching).

## Implications for taberna (concrete)

1. **Adopt this competition as taberna's surface-detection benchmark.** Same data,
   same metric. It lets us score a *classical* pipeline head-to-head against the
   nnU-Net field (~0.60–0.63).
2. **The winning post-processing is front-end-agnostic and 100% classical →
   build it in C regardless of detector.** Port: dust removal, 2×2×2 LUT
   hole-plug, binary_closing/fill_holes, Euler-number tunnel detection, PCA
   height-map repair (with edge flood-remove + thicken/erode unmerge),
   multi-threshold unmerge. This is taberna's "win condition."
3. **Touching sheets is the real frontier, and it's exactly taberna's thesis.**
   The field attacks merges *post hoc* (threshold sweeps + interpolation). taberna
   attacks them *at the segmentation level* via **signed affinity + Mutex
   Watershed** — repulsive edges between adjacent wraps from structure-tensor
   normal disagreement ([[taberna-segmentation]] `affinity.c`, `partition.c`).
   Preventing the merge beats repairing it. This is the experiment where taberna
   could actually be better, not just equivalent.
4. **PCA height-map repair ≈ local unwrap.** Same primitive as the Archimedean
   spiral fit, at patch scale → shared code with [[taberna-pipeline]] unwrap.
5. **Front-end swap experiment:** structure-tensor sheetness as per-voxel surface
   probability → winning post-proc → score. A fully no-ML entry.

## Sources

- 1st place (Tony Li, OzanM., Yiheng Wang, PaulG): nnU-Net ensemble + 5-step
  post-proc (dust, 2×2×2 LUT plug, height-map patch, closing, fill_holes).
- 2nd place (Duong Nguyen, Marius Heuser) "A postprocessing win": 128³+160³
  ensemble + multi-threshold PCA-interpolation unmerge.
- 4th place (Starry): huge ResEnc (7-stage) + PCA hole-fill + betti_matching local
  repair + Betti-2 metric hacking.
- 7th place: Coherence-Enhancing Diffusion (structure-tensor) front-end.
- Scroll Prize recap: https://scrollprize.substack.com/p/back-to-the-challenge-100k-kaggle
- Discussion: https://www.kaggle.com/competitions/vesuvius-challenge-surface-detection/discussion
