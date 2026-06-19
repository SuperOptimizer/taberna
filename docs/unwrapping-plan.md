# Taberna — Scroll Unwrapping Plan

*Plan of record for taberna's scroll-undeformation / virtual-unwrapping subsystem.
No-ML, classical methods, C, block-wise to TB scale, self-contained (does **not**
depend on matter-compressor or any villa application code). Companion to the SOTA
survey in [`segmentation-sota.md`](./segmentation-sota.md) — read that for method
citations; this doc is the decisions, the experiment order, the metrics, and the
human-annotation strategy.*

Date: 2026-06-18.

---

## 0. Goal

Take a deformed scanned CT volume of a spiral-wound papyrus scroll and undeform it
into an ~ideal **Archimedean spiral**, recording the spiral parameters **and an
invertible deformation field**, so the sheet can be unrolled for downstream ink
detection / reading. Account for truly lost material, varying sheet thickness, and
crush/twist/bend/break.

## 1. Strategy

**Independent classical line** (decided 2026-06-18). We are building our *own*
unwrapping pipeline. We are **not** reimplementing villa / Volume Cartographer /
VC3D or trying to match how anything under villa works — villa (and Paul's
`spiral-v2`) is an **idea source only**: we take what helps (the annotation
*concepts*, the `grad_mag` winding estimator, the loop-closure check, the
`satisfied_*` evaluation idea) and leave the rest (their data formats, their tools,
their app architecture). Nothing below requires a villa format or tool.

We borrow ideas and an evaluation harness but design our own architecture:

- **Front-end:** structure tensor → signed affinity graph → per-wrap segmentation
  via signed-graph partitioning (GASP / Mutex Watershed).
- **Backbone:** a block-wise **Eulerian winding-number field** (not a global mesh).
- **Fit:** least-squares `r = a + b·θ` per cross-section + an invertible
  deformation field, guarded by det-J > 0.

**Relation to Paul Henderson's `spiral-v2`** (the live SOTA in `ScrollPrize/villa`,
which is itself *already essentially no-ML* — PyTorch is just an autodiff/SGD engine
over a flow field, the only ML being the upstream surface-prediction U-Net we
replace with the structure tensor): Paul is a **reference + benchmark**, not a port
target. We:
- **borrow** his annotation taxonomy (umbilicus + winding-annotated point
  collections), his `grad_mag` strip-integral winding estimator, his winding-graph
  loop-closure check, and his exact `satisfied_*` metrics;
- **benchmark** by diffing our per-winding output meshes against his `spiral-v2`
  result and the Grand-Prize ground-truth mesh on the same region;
- **diverge** on architecture: signed-graph segmentation + Eulerian winding field,
  which is block-wise/streamable in a way his single global ~19 h fit is not.

**Hard constraints:** no ML; no matter-compressor; no villa application dependency;
every stage must run per-brick (256³ + halo) and stitch, because the full volume is
tens of TB (never resident).

## 2. Pipeline architecture

Per 256³ brick (+ halo), embarrassingly parallel, then stitched:

```
CT brick (+halo)
  └─[1] SHEET DETECTION (two-channel, classical)
        • coherence-enhancing diffusion pre-filter (capped iterations)
        • structure tensor → sheet NORMAL field e0 + coherence
          (integration scale ρ < inter-wrap gap, ~1–2 voxels)
        • OOF / Descoteaux → sheetness SCALAR (+ width); phase-symmetry
          channel for contrast-invariance (prototype)
        →  per-voxel: sheetness scalar + sheet normal
  └─[2] SIGNED AFFINITY GRAPH
        • (optional) SNIC supervoxels with sheetness-driven distance = coarsening;
          nodes carry centroid, normal, sheetness   [snic.c already exists]
        • edges: + attractive in-plane (edge ⟂ normal, normals agree)
                 − repulsive across-sheet (edge ∥ normal, esp. across a gap)
        • human must-link / cannot-link injected as hard edges (§3)
  └─[3] PARTITION → per-wrap labels
        • GASP average-linkage + cannot-link (preferred; robust to noisy
          structure-tensor affinities) or Mutex Watershed
  └─[4] WINDING FIELD (Eulerian, per-brick)
        • integrate structure-tensor azimuthal orientation → winding-phase scalar
          (Poisson / heat-method); absolute seeds = Dirichlet conditions
        • level sets = sheets; field value = winding coordinate
        ↓ persist SMALL supervoxel RAG + per-brick winding (out-of-core)
GLOBAL (on the reduced graph / field, out-of-core)
  └─[5] STITCH
        • delayed-merge (Lu/Zlateski/Seung) → chunking-invariant labels
        • accumulate generalized winding number across bricks (halo-consistent)
        • loop-closure holonomy check flags inconsistencies → annotation targets
  └─[6] ARCHIMEDEAN FIT + DEFORMATION FIELD
        • LSQ fit r = a + b·θ (winding coord vs accumulated angle, per cross-section)
        • build remap (observed → ideal spiral); inverse = recorded deformation field
        • det-J > 0 guard everywhere (invertibility)
  └─[7] UNROLL OUTPUT
        • per-winding meshes in our own format (optional one-off export to diff
          vs the GP mesh / Paul — validation only, not a compat requirement)
        • per-patch SLIM for final flat 2D rendering
```

Building blocks already present: `src/segmentation/snic.{c,h}` (supervoxel coarsening,
re-indexed to taberna layout). Everything else is new first-party code.

## 3. Human annotation & correction strategy

We will invest ~**a week of hand processing** up front; the design goal is to make
that effort *targeted* (the pipeline says where) and *reusable* (every annotation is
both a constraint and held-out validation).

### What to annotate (ranked by leverage / cost)

| # | Annotation | Cost | Leverage | Consumed by |
|---|---|---|---|---|
| 1 | **Umbilicus** `z→(y,x)` polyline | ~1 h | very high | θ-origin + radial dir for the winding field; `r` datum for the spiral fit |
| 2 | **Absolute-winding seeds** (point = "wrap N") | low | high (global) | Dirichlet seeds + integration-constant calibration of the winding field; partition seeds |
| 3 | **Relative-winding links** at breaks/touching (pair = "same wrap" / "Δ=1") | medium | **highest at hard spots** | must-link / cannot-link in the partition; Δ-winding constraints in the field; loop-closure |
| 4 | **Correction points** (placed post-fit, at errors, with winding) | low/point | high (targeted) | anchor/correction loss (collection-coupled); re-seed next solve |
| 5 | **Damage / lost-material masks** | low | medium | down-weight/extrapolate region in the field; exclude from metrics |
| — | dense sheet traces | high | — | **skip** — use existing VC/GP/instance-label data for validation instead |

### How corrections are used (four mechanisms)

1. **Partition constraints** — must-link (same wrap) and cannot-link (touching but
   distinct wraps) become hard edges in GASP/MWS, converting the silent
   wrap-merge failure into a satisfied constraint.
2. **Field boundary conditions** — absolute seeds are Dirichlet values; relative
   links are Δ-winding equality constraints in the Poisson solve.
3. **Fit anchors** — collection-coupled correction loss (Paul's "superconductor"
   averaging: points in a collection agree on winding up to their relative offsets).
4. **Validation** — hold out ~20 % of annotations; score TRE + `satisfied_*` on them.

### The targeted loop (the actual workflow)

The pipeline tells you *where* to annotate via three automatic "where's it broken"
signals:
- **loop-closure holonomy** — walk the winding graph; any cycle accumulating ≠ 0
  winding is inconsistent (Paul's `find_inconsistent_windings`);
- **merge hotspots** — regions of high `VOI_merge` vs current seeds;
- **det-J < 0 folds** in the deformation field.

**Week-1 budget:** Day 1 — umbilicus + an absolute-seed backbone (a few seeds per
~10 wraps along 2–3 radial spokes). Days 2–5 — run fit → tool surfaces flagged hard
spots → annotate only those with relative links / cannot-link → refit → repeat. Hold
out 20 % for validation.

### Annotation format & tooling

**Our own** simple format — we take the *concept* (umbilicus + winding-annotated
point collections + must/cannot-link) from prior art, not the schema:
- `umbilicus.txt` — one `z y x` per line.
- **point collections** — plain JSON: `{collection_id, winding_is_absolute:bool,
  points:[{id, zyx, wind}]}`; absolute collections calibrate, relative collections
  link. Plus an explicit `must_link` / `cannot_link` pair list for the partition.
- Designed so it's trivial to write/parse in C and to render as overlays.

Tooling is our choice, not tied to any villa tool. Week-1: place points on slices in
**whatever viewer is convenient** (or a 50-line slice-clicker we write) and export to
the JSON above; we render inconsistency-overlay images to direct attention. A polished
annotation UI is a later nicety, never a blocker. (If diffing against the GP mesh ever
wants their format, that's a one-off export, not a coupling.)

## 4. Data & ground truth

(Full inventory from research; see `segmentation-sota.md` references.)

- **Feasibility / fast iteration (no gating):** `instance-labels-harmonized.zip`
  (~1.4 GB) — Scroll 1, **256³ cubes, per-voxel sheet-INSTANCE labels**, `.nrrd`.
  Ideal for the de-risking experiment (validates sheet detection *and* separation).
- **Quantitative surface benchmark:** Kaggle "Vesuvius Challenge – Surface
  Detection" (2026) 256³/512³ chunks + binary sheet masks; official metric =
  `0.30·TopoScore + 0.35·SurfaceDice@τ + 0.35·VOI` (verify weights on the live page).
- **Unwrapping gold standard:** the Grand-Prize banner region + segment
  `20231231235900_GP` (`.obj` + `.ppm`); compare our flattening directly.
- **Spiral benchmark:** Paul's `spiral-v2` per-winding meshes on the GP region.
- **Umbilicus** reference polylines (ThaumatoAnakalyptor layout) for validating our
  own center estimate.

## 5. Volume ingest, compression & compute budget (vesuvius-c dropped as obsolete)

**The working volume is LOCAL, not streamed.** It is stored compressed with
matter-compressor (~10–25× at the codec quality setting τ≈8–16 — *distinct from the
SurfaceDice tolerance in §6*) and **fits on NVMe (PCIe Gen5)**. So random-access cube
reads are fast and there is **no remote-streaming problem to solve** — we deleted that
from scope. taberna consumes decompressed **u8 cubes** `(z,y,x,dims)`; how a cube is
produced (matter-compressor's external decode API, or a pre-extracted cube) is opaque
to us — we do not modify matter-compressor (it is being rewritten elsewhere).

**Lossy-compression budget (must design for it):** decompressed cubes carry
**~1 % MAE ≈ 2.5 graylevels on u8**. The sheet detector and structure tensor run on
*this* data, not pristine CT. That is non-trivial at the ~1-voxel inter-wrap gap,
where 2.5 levels of error can blur the very boundary we must not merge across — it
strengthens the case for the contrast-invariant phase-symmetry channel and a tight
structure-tensor integration scale, and it makes "raw vs compressed" a measured
robustness check (E2b).

**Public ground-truth data** (instance-labels, GP cubes) is separate and read as
`.nrrd` (hand-rolled ~30-line reader: text header + raw/gzip blob) or TIFF via the
`third-party/tiff` submodule. No zarr/S3 reader is needed for the planned work; if a
direct remote reader is ever wanted, revisit then (minimal zarr-v2 + `c-blosc2`, or
TensorStore).

**Compute budget (target workstation):** Ryzen 9 7950X (16C/32T), 192 GB RAM, 4 TB
NVMe Gen5, 2× RTX 3090 (48 GB total). Implications: a 256³ u8 cube is ~16 MB
(~67 MB as f32) — hundreds fit in RAM; block-wise stages parallelize across 32
threads with the full volume's working set comfortably resident. GPUs are **available
but not required**: the CPU pipeline is the baseline; the eventual global spiral fit
(stage 6) or heavy field solves may optionally use a 3090 if it pays off. Local NVMe
makes the block-wise I/O effectively free, so the scale test (E6) is compute-bound,
not I/O-bound.

## 6. Verification metrics

Tiered — cheap GT-free checks every iteration, heavier GT-based checks at gates.
Prefer permissive libs (trimesh/potpourri3d MIT; libigl MPL-2.0 — **avoid its CGAL
self-intersection module, GPL**); reimplement the small ones in C.

**Tier 0 — inner-loop, GT-free, every iteration:**
- **det-J > 0** everywhere (% foldings) — the field-level wrap-merge / fold detector;
  the single strongest invertibility guard.
- **self-intersection count** on the output mesh (0 = good).
- **radial winding monotonicity** — fraction of steps with dθ ≤ 0 along umbilicus
  rays ("MVF"; ideal 0).
- **loop-closure holonomy** — winding-graph cycles must sum to 0.

**Tier 1 — segmentation quality (vs instance labels):**
- **asymmetric VOI**: weight **VOI_merge ≫ VOI_split** (a cross-wrap merge is
  catastrophic; an intra-wrap split is recoverable). Report split/merge separately.
- **ARAND** precision (merges) vs recall (splits), per-wrap then aggregated.

**Tier 2 — surface detection (vs sheet masks):**
- **SurfaceDice@τ** (τ≈2 voxels), **TopoScore** (Betti-0/1 matching), precision/
  recall/F1, ASSD/HD95.

**Tier 3 — unwrapping topology (vs GP mesh / Paul):**
- **WJF** (winding-jump fraction), **single-manifold** checks (components=1,
  manifold edges/vertices, Euler χ), and Paul's **`satisfied_patches / satisfied_area
  / satisfied_unattached_pcls`** (radius band `|Δshifted_r| ≤ 0.5·dr` ∧ scan-space
  distance ≤ 4 voxels; patch satisfied at ≥95 % quads).

**Tier 4 — flattening distortion (intrinsic, GT-free):**
- **symmetric Dirichlet** (σ₁²+σ₂²+1/σ₁²+1/σ₂², floor 4), **quasi-conformal**
  k=σ_max/σ_min, **area** det J, geodesic-stress (libigl/SLIM).

**Tier 5 — registration (vs held-out landmarks):**
- **TRE** on held-out annotated landmarks; **% foldings / SDlogJ** on the field.

**End gate:** text legibility / ink signal on the unrolled result — the only true
end metric, but a final human gate, not an inner-loop signal.

## 7. Experiment plan (ordered, gated)

**E1 — Signed-affinity feasibility (highest priority; de-risks the whole bet).**
On `instance-labels-harmonized` 256³ cubes: build sheetness + normals → signed
structure-tensor affinities → GASP avg-linkage + cannot-link → score **asymmetric
VOI + ARAND**, focusing on *do adjacent wraps stay distinct in crushed regions*.
**Gate:** if VOI_merge is unacceptable even with a few cannot-link hints, the no-ML
signed-graph approach needs rethink before further investment.

**E2 — Sheet-detector bake-off.** OOF vs Descoteaux vs phase-symmetry vs
structure-tensor-only on faint / strong / near-touching wraps. Pick the scalar and
tune radius/σ to resolve one wrap without bridging (ρ < gap).

**E2b — Compression-robustness check.** Run E2 (and E1) on matter-compressor
round-tripped cubes (~1 % MAE / ~2.5 graylevels) vs the pristine GT cubes. Quantify
how much the lossy compression degrades sheet detection and wrap separation at the
inter-wrap gap. Informs whether we need the phase-symmetry channel and/or a tighter
codec quality τ for the regions we actually unwrap.

**E3 — Eulerian winding-field prototype.** On a multi-wrap region: solve the
winding-phase Poisson field, check monotonicity across wraps and **halo-stitch
consistency** of accumulated winding across a brick boundary (Tier-0 metrics).

**E4 — Delayed-merge stitch prototype.** Two adjacent bricks: verify global labels
are **chunking-invariant** (Lu/Zlateski/Seung). (Caveat: invariance is proven for
single-linkage/mean-affinity, *not* cleanly for globally-optimal MWS — open problem.)

**E5 — Spiral fit + deformation field on the GP region.** Fit `r=a+b·θ`, build the
invertible field (det-J > 0), output per-winding meshes; score **WJF + `satisfied_*`
+ Chamfer** vs the GP mesh and vs Paul's `spiral-v2`.

**E6 — Scale test.** Run the block-wise pipeline across a larger sub-scroll; record
wall-clock + peak RAM per brick; confirm bounded memory.

Human annotation (§3) is woven through E1, E3, E5 (umbilicus + seeds bootstrap;
targeted corrections refine).

## 8. Risks & open problems

1. **The central bet:** are structure-tensor *signed* affinities clean enough in
   crushed/delaminated/low-contrast regions? Unvalidated — E1 exists to test it.
2. **Sub-resolution gap, worsened by lossy compression:** adjacent wraps are
   ~1 voxel apart, at/below stable derivative-filter scale — and we read u8 cubes
   carrying ~1 % MAE (~2.5 graylevels) from matter-compressor, which can blur exactly
   that gap. OOF + phase-symmetry mitigate, human relative-links are the backstop, and
   E2b measures the actual damage.
3. **Single-object topology mismatch:** connectomics segments many closed objects;
   a scroll is one open sheet wound N times. Labeling each *wrap* distinctly via
   signed-graph agglomeration is genuinely novel work.
4. **Chunking-invariant stitching for signed-graph partitioners** is unsolved
   (proven only for single-linkage/mean-affinity).
5. **Joint deformation-field + graph maintenance** across re-merges — no precedent.
6. **Validation GT** for per-wrap spiral segmentation is limited; lean on
   instance-labels + GP mesh + held-out annotations.

## 9. Repo layout & milestones

```
src/segmentation/   snic.{c,h} [done]; sheet_tensor.{c,h} (structure tensor+OOF);
                    affinity.{c,h}; partition.{c,h} (GASP/MWS); stitch.{c,h}
src/unwrap/         winding_field.{c,h}; spiral_fit.{c,h}; deform.{c,h}
src/io/             nrrd.{c,h} (now); zarr.{c,h} (later)
src/eval/           metrics.{c,h} (VOI, surface-dice, det-J, WJF, satisfied_*)
src/annotate/       point_collection.{c,h}, umbilicus.{c,h} (our own simple JSON)
tools/              extract_cube, run_pipeline, eval, show_inconsistencies
docs/               segmentation-sota.md, unwrapping-plan.md
```

- **M0** ingest (NRRD) + a working feasibility cube.
- **M1** sheet detection (E2).
- **M2** signed-graph segmentation feasibility (E1) — **the gate**.
- **M3** winding field (E3) + stitch (E4).
- **M4** spiral fit + deformation + metrics vs GP / Paul (E5).
- **M5** scale (E6).

## 10. Decisions (resolved 2026-06-18)

- **Annotation UIs:** build good, purpose-built UIs as each task needs them — assume
  good tooling, not a constraint. (First likely tool: a slice viewer that places
  winding-annotated points + renders the pipeline's flagged-inconsistency overlays.)
- **GASP/MWS:** implement directly in C (small algorithms — sorted-edge union-find for
  MWS, priority-queue agglomeration for GASP); no C++ Higra dependency.
- **Compute:** 7950X / 192 GB / 4 TB NVMe Gen5 / 2× 3090. CPU pipeline is the
  baseline; GPU optional for stage 6 / heavy field solves. Working volume is local
  and fits compressed, so I/O is not a bottleneck.
- **Volume access:** local matter-compressor archive (external decode), consumed as
  u8 cubes; design for ~1 % MAE.
