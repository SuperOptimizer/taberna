# Touching-sheets improvement plan

Adapting the 2026-06-23 classical-SOTA research (see memory `touching-sheets-sota`) into taberna.
Goal: separate physically TOUCHING / FUSED scroll wraps — the one place the working winding pipeline
(`sheet_sep3d` + `tile_winding_mr` + `unroll_wind`) still merges two wraps into one (undercount +
winding leak), because the local signal `across·sheetness` is identical for a within-sheet radial
contact and an inter-wrap touch.

## The governing principle (from the research)

No classical operator invents a boundary the data didn't record. Every method that survives touches
injects a **non-local** signal. Three independent production systems (ThaumatoAnakalyptor, FASP/VC3D,
Henderson's diffeomorphic spiral) all use the **same** discriminator: attach each surface element a
**winding number** (unwrapped angle about the umbilicus) and gate merges on Δwinding. A within-sheet
contact has Δw≈0; an inter-wrap touch has Δw≈±1 (a full turn) — same local geometry, opposite global
gap. **We already have a coarse winding field and a SNIC supervoxel graph**, so the actionable fix is
cheap; the rigorous end-states (ordered max-flow, spiral fit) build on the same prior.

Plan is phased by **risk×yield**: instrument first, the near-free winding-gate next, then harden its
one assumption, then the heavyweight global optimizers, with detector refiners and the spiral north-star
as parallel tracks.

---

## Phase 0 — Instrument: measure touches before fixing them (foundation)

We currently have NO touch-specific metric; `backward-switch`/`climb` are whole-crop. Can't tune a
touch fix we can't see.

**0a. Multiscale sheetness + argmax-σ fused-wrap detector** (cheapest win in the whole report).
- `st_sheet_detect` runs at one `sigma_tensor≈1.5`. The two-peak (Sparrow) limit is σ_crit≈gap/2;
  if 1.5 > gap/2 at the working LOD, the tensor fuses two wraps *by itself*.
- Build `tools/sheetscale.c` (or a flag in `sheet_sep3d`): compute sheetness at a σ-set
  {0.7,1.0,1.5,2.0}, take per-voxel **max**, and emit the **argmax-σ map**. Voxels whose selected σ
  ≈ 2× the lone-sheet σ are fused-wrap candidates → that map *is* a touch detector and a ground-truth-
  free touch-density metric.
- Deliverable: a per-region "touch fraction" number + a PPM overlay. Feeds every later phase's eval.

**0b. Touch-region benchmark.** Extend `scripts/robustness_bench.py` with a `touch_frac` column (from
0a) and a "winding-monotonicity-across-touch" probe: along rays through high-touch-fraction voxels,
does the winding field still climb by ~1 per wrap, or does it flatten (leak)? This is the target metric
for Phases 1–3. Still ground-truth-free (we have no labelled wraps).

Exit: we can quantify, per region, how many touches exist and how badly the current winding leaks
through them.

---

## Phase 1 — Winding-gated supervoxel merge (Family A; the actionable fix)

Reuse `tools/svaff_seg.c` (SNIC supervoxels + per-supervoxel orientation + signed RAG + MWS). The wall
was that `across·sheetness` can't tell within-sheet from inter-wrap. Fix: add a per-edge **Δwinding**
and gate on it.

**1a. Geometric (non-circular) winding per supervoxel.** Do NOT only sample the coarse field (it may
have leaked through the very touch we're judging — circular). Compute, per supervoxel centroid, polar
coords about the umbilicus (we already auto-detect the center): radius r, angle θ. Then
`w_geom ≈ r/pitch + θ/2π` (pitch is already data-calibrated in `sheet_sep3d`). For two **radially
adjacent** supervoxels at ~equal θ, `Δw ≈ Δr/pitch`: one pitch apart in radius ⇒ different wraps ⇒
Δw≈1 — derived from radius+pitch (robust), **not** from the merged field. This is exactly Thaumato's
"angle-to-umbilicus" winding, made non-circular by the radius term.

**1b. Cross-check against the coarse field.** Also sample `cw_trilin`-style Δw from the coarse winding
volume (add `priorvol`/`priorlod` args to `svaff_seg`, mirroring `sheet_sep3d`). Agreement between
`Δw_geom` and `Δw_field` = trustworthy edge; disagreement flags a leaked/ambiguous region.

**1c. New merge rule.** Per RAG edge: if `|Δw| < 0.5` → attractive (`w = +k_attract·(1−across)`,
keep wraps' tangential continuity); if `|Δw| ≥ 0.5` → hard **repulsive/mutex** regardless of how
high `across·sheetness` is. This is the line that breaks the wall: a touch with identical local
geometry is forbidden because its Δw is a full turn. Feed to existing `mws_partition`.

**1d. Validate** on Phase-0 touch regions: do touching wraps now land in separate segments without the
collapse↔shatter failure (the sweep that killed plain affinity)? Render mid-z segments + check that
each concentric arc is one color and adjacent arcs differ.

Exit: touching wraps separated on real CT at supervoxel resolution, gated by a non-circular winding
number. Risk: if the coarse winding is *globally* wrong in a region (bad pitch/center), Δw is wrong —
mitigated by 1b's cross-check and Phase 2.

---

## Phase 2 — Harden the precondition: winding that survives touches

Phase 1's gate is only as good as the winding's turn-count *through* a touch. Two reinforcing tracks.

**2a. Preserve faint inter-wrap valleys in the sheetness field.** The research's strategic corollary:
a field that keeps even sub-noise valleys between fused crests is higher-leverage than any downstream
separator. Concretely: use the small-σ multiscale sheetness from 0a as the detection field (small σ
preserves the dip), and only use large σ for the *normal* estimate. Re-measure touch leak (0b).

**2b. Morse-Smale / persistence valley check.** Where 2a leaves a (possibly tiny) saddle between two
crests, persistence finds it at arbitrarily low contrast and *quantifies degree-of-fusion*. Prototype
with **DisPerSE** (pure C, drop-in spirit with our C submodules) on the sheetness field of a touch
region: does a 1-saddle / 2-separatrix exist between touching wraps? If yes → its separatrix surface is
a real inter-wrap cut we can inject as a `bar[]` barrier in `sheet_sep3d` (this time data-grounded, not
the failed `touchprom` guess). If the field is monotone across the touch → confirms truly-fused, defer
to Phase 1/3's order prior. Either outcome is informative.

Exit: the winding field climbs ~1/wrap through touches (0b metric green), OR we've proven the touch is
truly fused and must be carried by global order (Phase 3).

---

## Phase 3 — Ordered global optimization (Family B; the rigorous end-state)

Turn the soft winding-gate into a hard guarantee: K nested wraps, ordered, min-separation δ, forced
apart even at zero local contrast. K and δ come from the winding field; local sheetness only biases
placement.

**3a. Ishikawa / continuous-max-flow ordered-label regularization (lower-friction entry).** We already
have a winding field — reframe touch-fixing as: regularize the (possibly leaked) winding into the
nearest **monotone-outward integer wrap-index** label field, snapping its jumps onto sheetness ridges.
This is 3D phase-unwrapping-with-edges. Build as a new pass consuming `merged.f32`:
data term ties label to winding; pairwise term rewards label-increment at high-sheetness. Continuous
max-flow is matrix-free, tileable. Output: a wrap-index volume where touches are forced to increment.

**3b. LOGISMOS optimal-surface (heavier, explicit surfaces).** Columns = rays from umbilicus (or
winding-gradient integral curves); node cost = −sheetness; inter-surface ∞-arcs enforce
`δ_l ≤ f_i−f_{i+1} ≤ δ_u`. Solve with single-file Boykov–Kolmogorov max-flow. Use only if 3a's labels
aren't enough or we want sub-voxel ordered surfaces.

**Scale (both).** Global cut on full LOD0 is infeasible. Mandatory: (i) thin per-column/per-voxel
**bands** around the coarse winding init (shrinks the graph by ~100×), (ii) overlapping **tiles** +
stitch (already in `tile_winding_mr.py`), (iii) hard adjacency constraints instead of literal ∞
capacities.

Exit: wrap order is a structural guarantee; fused wraps provably separated within the banded region.

---

## Phase 4 — Detector refiners (parallel track; optional extra edge signals)

These don't replace Phase 1/3 — they add cheap independent evidence to the RAG edge weight or seed it.

- **4a. Fiber recto/verso orientation edge (novel, unexplored).** A second finer structure tensor →
  in-plane fiber orientation; a ~90° orientation jump across a contact localizes a boundary density
  can't. Inject as a soft orientation-discontinuity term into the `svaff_seg` edge weight. Proven analog:
  paperboard-ply separation under low density contrast. Blind where wraps are co-oriented; degrades in
  compressed core — budget as partial assist.
- **4b. OOF targeted refiner.** Optimally Oriented Flux localizes closely-spaced parallel sheets the
  structure tensor merges (needs gap≳2 vox; sub-PSF hopeless). FFT-based → NOT tileable → run only on
  Phase-0-flagged slab ROIs as a second opinion, not whole-volume.
- **4c. Power-watershed seeded carving.** Carves a unique contrast-invariant divide through a signal-
  free zone when seeded both sides — seed from where wraps ARE separated (low touch-fraction) and let
  it propagate the cut into the touch. Reference C/C++ exists.

Exit: extra orthogonal cues available for the hardest regions; each kept only if it moves the 0b metric.

---

## Phase 5 — North star: diffeomorphic Archimedean spiral fit

Henderson WACV 2026. Fit a topology-preserving spiral (affine∘flow∘gap-scale diffeomorphism) globally
to the sheetness + winding field; fused wraps **cannot** merge because the map stays injective. Makes
the local touch decision moot entirely and *is* the scroll-un-deformation goal already in memory
`taberna-segmentation`. Largest build; do after Phases 1–3 prove the order prior works at small scale.
Losses portable to C: L_normal (fit sheetness normals), L_radius (constant canonical radius/wrap),
L_windings (K-apart in data ⇒ K-apart in canonical). Reuses our umbilicus, pitch, winding field.

---

## Cross-cutting: validation harness (ground-truth-free — we have no labelled wraps)

Every phase reports against Phase-0 metrics:
1. **Touch fraction** (argmax-σ map, 0a) — how much touch exists in a region.
2. **Winding-monotonicity through touches** (0b) — climb ≈1/wrap along rays crossing touches; flat = leak.
3. **Segment-per-wrap** (Phase 1) — each concentric arc one segment, adjacent arcs differ; no
   collapse↔shatter.
4. **Unroll quality** (`unroll_wind`) — touches that previously smeared two wraps into a half-wrap band
   now render as two distinct flattened sheets; per-wrap mean flat across all wraps.
5. **Cross-method agreement** — `Δw_geom` vs `Δw_field` (1b); persistence saddle existence (2b);
   independent Poisson winding (`wind_poisson`).

Beware the known traps from `sheet-sep-spiral-fit`: metrics blind to *flattening* (a collapsed field
has no backward-switches) — always pair any separation claim with the fill-climb / monotonicity guard.

---

## Sequencing & dependencies

```
Phase 0 (instrument)  ──► Phase 1 (winding-gate merge)  ──► Phase 3 (ordered max-flow)
   │                          │                                  ▲
   └─► 2a (small-σ field) ────┴──► 2b (persistence valley) ──────┘
                                                                  
Phase 4 (refiners)  : parallel, feeds 1 & 3 edge weights
Phase 5 (spiral fit): after 1–3 validate the order prior
```

**Build first:** Phase 0 (no fix is tunable without it) → Phase 1 (near-free, reuses `svaff_seg` +
existing winding field, exactly what the production systems ship). That two-step is the minimum path to
a real result; everything after is hardening and rigor.

## Risks & kill-criteria

- **Phase 1 circularity** — coarse winding leaked through the touch ⇒ Δw_field wrong. Mitigation: the
  geometric Δw (1a) is the primary; field is cross-check only. Kill if Δw_geom is itself unreliable
  (bad center/pitch) over touch regions even after Phase 2.
- **Phase 3 scale** — if banding+tiling can't bring the max-flow graph into memory at the needed LOD,
  fall back to Phase 1's soft gate as the production path.
- **Truly sub-PSF fused touches** — a sampling hard limit no method escapes (all four streams agree).
  Managed by the order prior (carry the trace), not solved; document where this floor bites.
- **Detector refiners (Phase 4)** — each must move the 0b metric or it's dropped; no speculative knobs
  (cf. the tested-negative `touchprom`).
