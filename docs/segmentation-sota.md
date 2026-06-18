# Taberna: SOTA for Undeforming a Scroll to an Archimedean Spiral (No-ML)

*Synthesis for the taberna maintainer. Goal: take a deformed CT volume of a spiral-wound papyrus sheet and undeform it into an ~ideal Archimedean spiral, recording the spiral + an invertible deformation field, at tens-of-TB scale, in chunks/bricks with halos, using classical (non-ML) methods. ML is reported where it is the bar to beat or where its non-ML core is borrowable.*

---

## 1. Executive summary

**What the SOTA actually is.** Two facts dominate everything below:

1. **The published method that literally does taberna's job already exists.** Paul Henderson's *Diffeomorphic Spiral Fitting* (WACV 2026, [arXiv 2512.04927](https://arxiv.org/abs/2512.04927), code [pmh47/spiral-fitting](https://github.com/pmh47/spiral-fitting)) globally fits an idealized **extruded Archimedean spiral** to a scroll via an explicit, invertible, topology-guaranteeing **diffeomorphism** (per-slice affine + ODE-integrated flow field + radial gap-scaling), optimized against surface-normal / constant-radius / winding-spacing / fiber-orientation losses. **Its entire geometric core is classical**; only the per-voxel surface/fiber inputs come from three nnU-Nets — exactly the part taberna replaces with its structure-tensor sheetness detector. This is the single most important reference in this report. (Caveat, *high confidence*: it is a single **global** ~19 h RTX 3090 optimization on a ~30-winding subvolume, **not** block-wise — the out-of-core re-engineering is genuinely new work.)

2. **For the front-end (sheet detection -> graph -> stitching), the connectomics field has spent a decade building precisely taberna's pipeline at TB/PB scale, and its segmentation core is almost entirely classical** — only the affinity prediction is learned. You over-segment into supervoxels, build a region-adjacency graph (RAG) with geometry on nodes and affinity on edges, agglomerate greedily, and stitch blocks with a delayed-merge rule that makes the global result independent of the chunking. Swap learned affinities for structure-tensor affinities and the no-ML constraint is fully compatible with this SOTA core.

**The 3–5 highest-leverage takeaways:**

- **The lever is the affinity/sheetness field, not the clustering algorithm.** Every compact superpixel/supervoxel method (SLIC, SNIC, SEEDS, VCCS) will leak across adjacent wraps or fragment the sheet unless the distance/edge weight is computed on a **sheetness + orientation** field, not raw CT intensity. Spend your effort there.

- **Use a two-channel classical primitive for sheet detection:** the **structure tensor** for the robust orientation/normal field (first derivatives + smoothing, noise-robust, you already have it), and the **Hessian-based Descoteaux sheetness** or **Optimally Oriented Flux (OOF)** for the mid-plane-peaking sheetness scalar + width. OOF specifically resolves *closely-located adjacent structures* the Hessian merges — directly relevant to ~1-voxel inter-wrap gaps ([Law & Chung ECCV 2008](https://cse.hkust.edu.hk/~achung/eccv08_law_chung.pdf)).

- **The scroll's defining constraint is that adjacent wraps touch but must NOT merge.** This makes **signed-graph partitioning** (attractive in-plane edges, *repulsive* across-sheet edges) the right tool — Mutex Watershed / GASP — and it is what prevents the catastrophic failure mode of the whole scroll collapsing into one segment.

- **Copy the connectomics block-wise architecture wholesale** (per-brick supervoxel + local agglomeration -> persist a small supervoxel RAG out-of-core -> global partition on that reduced graph -> emit a relabel LUT). The two-level local-then-global pattern is proven: E11 Bio's **Volara** (May 2025) does exactly this with Mutex Watershed; **Lu/Zlateski/Seung 2021** did it at >1.5 trillion edges with a chunking-invariant delayed-merge rule.

- **For the undeform itself, prefer Eulerian field methods over a global mesh.** Solve a Laplace/Poisson harmonic field per brick for the across-winding coordinate, integrate the structure-tensor orientation into a winding-phase scalar, accumulate a **winding-number field** (the natural cross-brick-stitchable coordinate), then fit `r = a + b·θ` by trivial linear least squares. This is block-wise-friendly in a way Henderson's global optimizer is not.

---

## 2. Findings by axis

Each axis below has *high* confidence per adversarial fact-check unless flagged otherwise.

### 2.1 3D supervoxel / superpixel oversegmentation

**SOTA:** Classical oversegmentation is mature; the favored front-end at TB scale is **seeded watershed on an affinity/sheetness map + RAG agglomeration** (the connectomics workhorse), not compact superpixels. Among compact generators, two 2024 surveys reframe the field: the peer-reviewed **ACM Computing Surveys taxonomy** ([arXiv 2409.19179](https://arxiv.org/abs/2409.19179), benchmarks 20 methods on 9 criteria incl. connectivity/running-time/delineation/stability) and the *"ill-posed problem"* survey ([arXiv 2411.06478](https://arxiv.org/pdf/2411.06478), regularity-vs-adherence is a fundamental conflict; irregular methods can game boundary-recall).

**Methods that matter:**

| Method | C/ML | Verdict for taberna |
|---|---|---|
| **SNIC** (CVPR 2017, [paper](https://openaccess.thecvf.com/content_cvpr_2017/papers/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.pdf)) | classical | Best fit for your *existing* direction: non-iterative single pass, **connectivity guaranteed by construction**, lower memory than SLIC, trivial 3D + intensity-only. **But** the global priority queue is inherently sequential (run per-brick + stitch), and it will under-segment across a thin sheet unless the distance is computed on the sheetness field. |
| **3D SEEDS** (2025, [arXiv 2502.02409](https://arxiv.org/abs/2502.02409), [code](https://github.com/Zch0414/3d_seeds)) | classical | Strongest recent **volumetric no-ML** result: 10× faster than 3D SLIC, +6.5% Dice, open source, designed for grayscale medical CT. Local boundary updates parallelize better than SNIC. **Downside:** blocky/regular supervoxels fight thin-sheet adherence. Worth benchmarking as a drop-in. |
| **SLIC / SSLIC** (PAMI 2012; [SSLIC arXiv 1806.08741](https://arxiv.org/abs/1806.08741)) | classical | The baseline; SSLIC is a verified n-D/3D parallel implementation tested on the Visible Human volume. Now dominated by SNIC (memory, connectivity) and 3D SEEDS (speed). Connectivity is a post-hoc fix that can orphan thin fragments. |
| **ERS** (CVPR 2011, [TR](https://www.merl.com/publications/docs/TR2011-035.pdf)) | classical | Highest classical boundary adherence (1/2-approx guarantee) — but a global graph optimization with **no streaming variant**. Use only as a small-ROI quality ceiling. |
| **ISF / DISF / SICLE** ([DISF arXiv 2007.04257](https://arxiv.org/pdf/2007.04257)) | classical | Path-cost/seeded-forest with an adherence-maximizing connectivity function is conceptually the best classical match for thin high-contrast sheets and integrates naturally with affinity weights — but reference code is 2D/in-core; *research-prototype*. |
| **Seeded watershed + RAG agglomeration** (GALA, waterz, mean-affinity) | **hybrid** | The actual SOTA paradigm at TB/PB scale. Fully no-ML if you keep the affinity hand-crafted (structure tensor) and use mean-affinity agglomeration. **This is what your sheetness + SNIC RAG already approximates.** |
| VCCS (CVPR 2013), StreamGBH (ECCV 2012), gSLICr (GPU) | classical | Borrow ideas, not implementations: VCCS's normal-aware distance; StreamGBH's *process-once/freeze/stitch-at-boundary* streaming discipline; GPU SLIC as a throughput optimization only. |
| Deep methods (SSN, AINet, SEAL, LNS-Net, SAM+maskSLIC) | **ml** | Out of scope; borrowable core is the learned affinity map, which you replace with structure-tensor sheetness. Report-only. |

**Scalability verdict:** Compact methods are all per-brick feasible with halos; **none** ships out-of-core cross-chunk stitching — you build it. Watershed+RAG is the canonical proven out-of-core design.

### 2.2 Affinity-graph / agglomerative segmentation of thin sheets at scale (connectomics)

**SOTA:** The connectomics recipe — affinity graph -> watershed fragments -> hierarchical RAG agglomeration — is classical at its core (the graph algorithms don't care whether weights came from a CNN or a structure tensor). Parameter-free signed-graph partitioning via **Mutex Watershed / GASP** is the segmentation SOTA; **block-wise agglomeration with cross-block RAG stitching** (daisy+waterz, Volara) is the engineering SOTA.

**Methods that matter:**

| Method | C/ML | Verdict |
|---|---|---|
| **Mutex Watershed** ([arXiv 1904.12654](https://arxiv.org/abs/1904.12654), [code](https://github.com/sciai-lab/mutex-watershed)) | classical | Arguably the single best fit. **Attractive** (in-plane coherence) + **repulsive** (across-sheet, normal-direction) edges, **no seeds, no threshold**, near-linearithmic union-find. The repulsion is exactly the mechanism that keeps touching wraps apart. |
| **GASP** ([CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Bailoni_GASP_a_Generalized_Framework_for_Agglomerative_Clustering_of_Signed_Graphs_CVPR_2022_paper.pdf), [code](https://github.com/abailoni/GASP)) | classical | The design-space map. Proves AbsMax+constraints = Mutex Watershed; average-linkage = mean-affinity. **Average-linkage-with-cannot-link was most accurate AND most noise-robust** on CREMI — important because structure-tensor affinities are noisier than learned ones. Use this as your experimentation substrate. |
| **waterz** ([code](https://github.com/funkey/waterz)); **Zlateski/Seung quasilinear watershed** ([arXiv 1505.00249](https://arxiv.org/abs/1505.00249)) | classical | Battle-tested unsigned watershed+mean-affinity. Fast fragment generator/baseline, **but unsigned + threshold-dependent — cannot separate touching wraps by itself.** |
| **Volara** (E11 Bio, May 2025, [blog](https://www.e11.bio/blog/volara), [code](https://github.com/e11bio/volara)) | **hybrid** (classical core) | **The current reference design for exactly taberna's need:** per-block MWS supervoxels -> supervoxel-center nodes in PostgreSQL -> aggregated within/cross-block edge costs -> **global** MWS on the supervoxel graph -> relabel LUT. Linear scaling 1->120 EC2 workers. Study this repo. |
| **daisy + LSD blockwise pipeline** ([Nature Methods 2022](https://www.nature.com/articles/s41592-022-01711-z), [daisy](https://github.com/funkelab/daisy)) | classical engineering | Checkerboard block scheduling; persist a supervoxel RAG (not voxel labels) across blocks; resolve identity via relabel LUT. Adopt wholesale. |
| **Multicut / lifted multicut** ([Beier et al. 2017](https://hci.iwr.uni-heidelberg.de/sites/default/files/publications/files/217205318/beier_17_multicut.pdf), [lifted arXiv 1505.06973](https://arxiv.org/pdf/1505.06973)) | classical | Most globally consistent; natural place for "these fragments are different wraps" priors. NP-hard; reserve for a **final cleanup on the small inter-brick supervoxel graph**, never the voxel grid. |
| **GALA** ([arXiv 1303.6163](https://arxiv.org/pdf/1303.6163)) | **ml** | Learns the agglomeration *policy*. Skip the classifier; **borrow its hand-crafted RAG edge features** (boundary intensity stats, region size/shape) as deterministic merge scores. |

**Scalability verdict:** Excellent and proven at TB/PB scale. MWS/GASP have no native distributed global optimum; the Volara two-level (local-then-global-on-supervoxel-graph) trick is the empirically-proven workaround.

### 2.3 Structure tensor & classical sheet/surface detection

**SOTA (caveat, honest):** The *actual* SOTA on Vesuvius scrolls is **ML** (nnU-Net surface/fiber prediction; this is the bar). Among classical options, no single operator gives both a clean scalar and a polarity-aware normal — use **two channels**.

**Methods that matter:**

| Method | C/ML | Role |
|---|---|---|
| **Structure tensor** (Förstner/Bigün/Weickert) | classical | **Normal/orientation field.** Most noise-robust (first derivatives + averaging). Polarity-blind and double-responds at the two surfaces, so it is an orientation estimator, not a localizer. **Critical tuning: integration scale ρ must be strictly below the inter-wrap gap (~1–2 voxels)** or you average two wraps' normals. You already have this. |
| **Descoteaux sheetness** (Hessian, [PubMed 17127650](https://www.ncbi.nlm.nih.gov/pubmed/17127650); [ITK impl](https://github.com/InsightSoftwareConsortium/LesionSizingToolkit/blob/main/include/itkDescoteauxSheetnessImageFilter.h)) | classical | **Sheetness scalar + mid-plane localization + width + normal** in one eigen-decomposition. Production ITK filter. Noise-fragile (2nd derivatives); keep σ small. |
| **Optimally Oriented Flux (OOF)** ([Law & Chung ECCV 2008](https://cse.hkust.edu.hk/~achung/eccv08_law_chung.pdf)) | classical | **Best lever for inter-wrap bleed:** flux over a sphere of radius r; the paper explicitly shows it resolves closely-located structures the Hessian merges, at **no extra compute** over the Hessian. r maps directly to "stay within one wrap." Original framing is tubular; sheet is the complementary eigenvalue signature. |
| **Beyond Frangi (Jerman)** / **strain-energy filter** ([SPIE 2015](https://ui.adsabs.harvard.edu/abs/2015SPIE.9413E..2AJ/abstract); [MedIA](https://pubmed.ncbi.nlm.nih.gov/20879421/)) | classical | Parameter-light Frangi successors. **Avoid vanilla Frangi** — its β/c constants are contrast-dependent and won't transfer across TB-scale bricks. |
| **Phase congruency / symmetry (Kovesi)** ([2000](https://pubmed.ncbi.nlm.nih.gov/11195306/); [3D C++ impl](https://github.com/chvillap/phase-congruency-features)) | classical | **Contrast-invariant** sheetness scalar — fires equally on faint and strong wraps, addressing taberna's varying-thickness requirement (the Achilles heel of all amplitude-based filters). 3D is fiddly/expensive; *medium confidence on real-scroll behavior — prototype it.* |
| **Coherence-Enhancing Diffusion (CED)** ([Weickert IJCV 1999](https://link.springer.com/article/10.1023/A:1008009714131)) | classical | Structure-tensor-driven pre-filter: denoise *along* the wrap, preserve the cross-wrap edge. **Cap iterations tightly** or it bridges adjacent wraps. (Not in the public Vesuvius pipeline; a recommended addition.) |
| **nnU-Net surface/fiber + Sobel normals + spiral fit** (Henderson, Schilling/Johnson) | **hybrid/ml** | The ML bar. Reusable without ML: (a) the global diffeomorphic spiral fit (operates on paths+normals, not voxels); (b) 3D-Sobel-on-probability is just a crude gradient the structure tensor does better. ML wins on the raw scalar in low-contrast regions — the gap a classical pipeline must close. |

**Scalability verdict:** All classical primitives are separable convolution + per-voxel eigensolve — embarrassingly parallel per brick, halo ≈ max(3σ, OOF radius, CED diffusion length).

### 2.4 Vesuvius scroll segmentation & virtual unwrapping (the application SOTA)

**SOTA:** The community explicitly factors the problem into **Representation** (predict sheet surface + fiber orientation per voxel — *"via ML or classic filters"*) and **Fitting** (geometry turning noisy predictions into a continuous manifold). The fitting half is almost entirely classical. **No automated method has yet matched the 2023 manual text recovery** — the 2024 $200K Grand Prize went *unclaimed*, and both autosegmentation submissions received $30K Gold Aureus++ awards rather than the full First Automated Segmentation Prize ([prize writeup](https://scrollprize.substack.com/p/awarding-the-amazing-autosegmentation)).

**Methods that matter:**

| Method | C/ML | Verdict |
|---|---|---|
| **Diffeomorphic Spiral Fitting** (Henderson, [arXiv 2512.04927](https://arxiv.org/abs/2512.04927), [code](https://github.com/pmh47/spiral-fitting)) | **hybrid** (classical core) | **This is taberna's stated end-goal, published.** Outputs an ideal Archimedean spiral + an invertible deformation field + topology guarantee. **Steal:** (1) the sparse 1D-path representation (threshold@50% -> connected components -> skeletonize -> longest-path) that lets the optimizer avoid the TB volume; (2) the flow-integrated diffeomorphism guaranteeing one contiguous sheet through gaps; (3) constant-radius + winding-spacing priors. Inputs replaceable by structure-tensor. *Global, not block-wise (~19h/RTX3090).* |
| **FASP / VC3D** (Schilling & Johnson, [FASP repo](https://github.com/hendrikschilling/FASP), [VC3D](https://github.com/ScrollPrize/villa/tree/main/volume-cartographer)) | **hybrid** | Production tracer: greedy **quad-patch growth + windowed least-squares + Ceres** in a **moving window** (~linear in surface area) → the closest existing analog to block-wise/halo processing. ML confined to surface prediction + ink detection. |
| **ThaumatoAnakalyptor** ([repo](https://github.com/schillij95/ThaumatoAnakalyptor)) | **hybrid** | Validates several taberna primitives: derivative-based surface-point extraction (≈ structure-tensor response), **umbilicus + winding-number-per-patch** bookkeeping (the canonical spiral-ordering mechanism), Poisson meshing, **chunked C++/GPU graph solver across subvolumes** (battle-tested at-scale stitching). Patch instance-seg is ML; stitching is heuristic. |
| **Surface Tracer** ([unwrapping](https://scrollprize.org/unwrapping)) | **hybrid** | Local mesh growth with a **fidelity-vs-bending objective** — a clean ML-optional template: structure-tensor planarity = data term, thin-sheet bending penalty = smoothness term. Block-friendly. |
| **Volume Cartographer + Optical Flow Segmentation** (Seales/EduceLab, [En-Gedi Sci Adv 2016](https://www.science.org/doi/10.1126/sciadv.1601247)) | classical | Foundational 3-stage pipeline (segment/flatten/texture). Use the **texturing/rendering backend**; the 2D slice-propagation tracer is too manual/brittle. |
| **nnU-Net Representation** (Kaggle Surface Detection 2026, [comp](https://www.kaggle.com/competitions/vesuvius-challenge-surface-detection)) | **ml** | The component taberna replaces with structure-tensor. **Steal the topology-aware metric** (SurfaceDice@τ + VOI split/merge + TopoScore) for your own validation. |

**Scalability verdict:** Henderson is global (sparse-path trick decouples memory from volume size but the optimizer needs domain-decomposition rework); FASP moving-window and Thaumato chunk-graph already operate at scroll scale but stitch heuristically.

### 2.5 Surface flattening / parametrization / coordinate-field & winding-field unrolling

**SOTA:** Mesh parametrization is mature (LSCM, ARAP, SLIM, BFF, CETM/CEPS) but **none scales globally to a whole TB scroll** — confine to per-patch final layout. For taberna's constraint the promising direction is **Eulerian field methods** (no global mesh): harmonic coordinate fields + orientation-field-to-scalar Poisson integration + a global winding-number field, then a closed-form Archimedean fit.

**Methods that matter:**

| Method | C/ML | Role |
|---|---|---|
| **Eulerian harmonic/Laplace coordinate field** ([Yezzi & Prince TMI 2003](https://lccv.ece.gatech.edu/docs/tmi_EulerianA_Tissue_Thickness.pdf)) | classical | **The TB-scale-friendly primitive.** Sparse linear solve on the regular grid (multigrid + domain decomposition + halo exchange); level sets = sheet layers, streamline length = local sheet spacing (≈ thickness). No mesh, no global correspondence. *A single Poisson coordinate is monotone over only a few wraps — multi-revolution still needs winding accumulation.* |
| **Orientation-field → scalar via Poisson (Heat-Method integration)** ([Crane et al.](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/); [parallel/scalable Tao 2019](https://arxiv.org/abs/1812.06060)) | classical | Integrate the structure-tensor azimuthal director to a winding-phase scalar with one Poisson solve. Director sign-ambiguity must be propagated across halos. |
| **Generalized / relative winding number** ([Jacobson SIGGRAPH 2013](https://igl.ethz.ch/projects/winding-number/robust-inside-outside-segmentation-using-generalized-winding-numbers-siggraph-2013-compressed-jacobson-et-al.pdf); [one-shot 2024](https://arxiv.org/abs/2408.04466)) | classical | **The natural global unrolling coordinate** — accumulated winding / 2π indexes the wrap and is the cross-brick-stitchable scalar. Robust-to-holes matches missing/crushed material. Fast-summation evaluable per brick. |
| **Direct Archimedean-spiral least-squares** ([Mishra 2004](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=531542)) | classical | The final cheap step: fit `r = a + b·θ` (tiny LSQ per cross-section); per-voxel residual = in-plane deformation. |
| **SLIM / slim-flatboi** ([SLIM TOG 2017](https://cims.nyu.edu/gcl/papers/SLIM2017.pdf); [slim-flatboi](https://github.com/giorgioangel/slim-flatboi)) | classical | **What Vesuvius already uses** to flatten individual segments (flip-free, millions of faces). Use for **final per-patch 2D rendering only**, not the global engine. |
| LSCM / ARAP / BFF / CETM-CEPS / Ricci flow | classical | Per-patch flattening alternatives (ARAP isometric for near-developable papyrus; BFF for boundary control; CEPS/Ricci for topological robustness after damage). Per-mesh, not out-of-core. |
| **Field-aligned integer-grid maps** (Bommes/Campen) | classical | Borrow the *concept* (seamless integer transition = +1 winding increment per wrap), **not** the NP-hard MIQP solver. |
| **Conformal colon/cortical unrolling (holomorphic 1-form)** ([Hong/Gu 2006](https://dl.acm.org/doi/10.1145/1128888.1128901)) | classical | Methodological analog: "cut tube, map to rectangle" is the mesh sibling of integrating two orthogonal coordinate fields. Inspiration, not a TB tool. |

**Scalability verdict:** Eulerian PDE fields + winding-number + closed-form spiral fit are the only combination that is genuinely block-wise. Mesh parametrization is per-patch only.

### 2.6 Scalable graph representations (RAG, merge trees, out-of-core stores, multiresolution)

**SOTA:** Over-segment per block -> RAG (geometry on nodes, affinity on edges) -> priority-queue HAC -> chunking-invariant block stitching, with the RAG stored out-of-core (KV/doc store) and a merge-tree kept for cheap re-thresholding.

**Methods that matter:**

| Method | C/ML | Role |
|---|---|---|
| **Block-wise distributed agglomeration with delayed merges** ([Lu/Zlateski/Seung 2021, arXiv 2106.10795](https://arxiv.org/abs/2106.10795)) | classical | **The single most relevant method here.** Defer any boundary-crossing merge until both regions live in one (recursively doubled) chunk → final result **provably identical to non-chunked**. Proven at >1.5T edges / 135B supervoxels. This is your seam-free stitching rule. |
| **PyChunkedGraph / CAVE** ([repo](https://github.com/CAVEconnectome/PyChunkedGraph/blob/pcgv1/README.md)) | classical | **Octree-of-chunks dynamic RAG** for O(local) edits/queries — the model for a persistent, editable global spiral-graph aligned to spatial chunks + LODs. |
| **Size-dependent single-linkage on watershed basin graph** ([Zlateski & Seung 2015](https://arxiv.org/abs/1505.00249)) | classical | Canonical quasilinear agglomeration kernel; parameter-light size knob; produces a re-cuttable hierarchy. |
| **Mutex Watershed / GASP** (see 2.2) | classical | Signed-graph partitioning — encode "same sheet" vs "different wrap" as ± edges. |
| **Higra** ([SoftwareX 2019](https://hal.science/hal-02309938/document), [repo](https://github.com/higra/Higra)) | classical | **Most directly reusable codebase given taberna is in C:** RAG, MST, watershed cuts, binary partition trees. In-memory → the in-block engine, not the global store. |
| **Out-of-core BPT — Join/Select/Insert** ([DGMM 2022 arXiv 2210.02218](https://arxiv.org/abs/2210.02218); [JMIV 2025](https://link.springer.com/article/10.1007/s10851-025-01234-0)) | classical | For when the merge tree itself exceeds RAM: tile-and-join per-brick BPTs out-of-core. |
| **TeraHAC** ([SIGMOD 2024, arXiv 2308.03578](https://arxiv.org/abs/2308.03578)); RAC/SCC | classical | If you want exact/near-exact **global** agglomeration instead of the delayed-merge approximation: (1+ε)-approx HAC to 8 trillion edges. Constrains you to reducible (mean/Ward-like) linkage. |
| Irregular pyramids / **Binary Partition Tree** ([Salembier & Garrido TIP 2000](https://ieeexplore.ieee.org/document/841934/)) | classical | The classical multiresolution-graph-pyramid; modern practice favors a single BPT/merge-tree as the pyramid. |
| **funlib/LSD MongoDB RAG**; **GALA** | hybrid | Out-of-core RAG store design (KV/doc, keyed by chunk); GALA's RAG-feature-update loop. |

**Scalability verdict:** This is the most thoroughly solved axis. Delayed-merge stitching + out-of-core RAG store + persistent merge-tree LUT is the proven workstation-to-cluster path.

---

## 3. Recommended architecture for taberna

A four-stage, block-wise, no-ML pipeline. Rationale follows each stage.

```
                 ┌─────────────────────────────────────────────────────────────┐
                 │  Per 256³ brick (+ halo), embarrassingly parallel             │
                 │                                                               │
  CT brick  ──►  │  [1] SHEET DETECTION (two-channel, classical)                 │
                 │      • CED pre-filter (capped iters)                          │
                 │      • Structure tensor → normal/orientation field (ρ < gap)  │
                 │      • OOF (or Descoteaux) → sheetness scalar + width         │
                 │      • (prototype) phase-symmetry → contrast-invariant scalar │
                 │            ↓ surface-prob + normals + fiber dirs              │
                 │  [2] LOCAL GRAPH/FIELD                                        │
                 │      • SNIC supervoxels (sheetness-driven distance) = RAG     │
                 │        nodes: centroid, size, normal, sheetness               │
                 │      • signed edges: +in-plane coherence, −across-sheet gap   │
                 │      • GASP avg-linkage+cannot-link (or MWS) local partition  │
                 │      • Eulerian Laplace field → across-winding coord + spacing │
                 └─────────────────────────────────────────────────────────────┘
                            ↓ persist SMALL supervoxel RAG (PostgreSQL/CSR-on-disk)
                 ┌─────────────────────────────────────────────────────────────┐
                 │  [3] GLOBAL STITCH (on reduced supervoxel graph, out-of-core) │
                 │      • Lu/Zlateski/Seung delayed-merge → chunking-invariant   │
                 │      • global GASP/MWS partition (Volara two-level pattern)   │
                 │      • integrate orientation→winding-phase (Poisson), accum.  │
                 │        generalized WINDING-NUMBER FIELD across bricks         │
                 │      • merge-tree / BPT kept for cheap re-thresholding        │
                 └─────────────────────────────────────────────────────────────┘
                            ↓
                 ┌─────────────────────────────────────────────────────────────┐
                 │  [4] ARCHIMEDEAN FIT + DEFORMATION FIELD                      │
                 │      • LSQ fit r = a + b·θ  (winding-coord vs accum. angle)   │
                 │      • Henderson-style composable diffeomorphism, tiled on    │
                 │        coarse grids (per-slice affine + Euler-flow + gap-scale)│
                 │      • output: spiral params + invertible deformation field   │
                 │      • SLIM/slim-flatboi for final per-patch 2D rendering     │
                 └─────────────────────────────────────────────────────────────┘
```

**Adopt:**

- **The connectomics block-wise RAG + Volara two-level + Lu/Zlateski/Seung delayed-merge stack.** Rationale: it is the only proven path to TB scale with seam-free global consistency, and its core is classical. Study [e11bio/volara](https://github.com/e11bio/volara) and [funkelab/daisy](https://github.com/funkelab/daisy) as concrete references.
- **Signed graphs (MWS/GASP) with structure-tensor affinities.** Rationale: the scroll's defining constraint (touching wraps must not merge) *is* a repulsive-edge problem. Prefer **GASP average-linkage + cannot-link** over raw MWS because structure-tensor affinities are noisier than learned ones (CVPR 2022 robustness result). Higra gives you C/C++ primitives to build on.
- **Henderson's diffeomorphic Archimedean spiral fit as the undeform model**, fed classical fields instead of nnU-Net. Rationale: it is literally taberna's goal, published, with a topology guarantee and an invertible field — exactly what you must record.
- **Eulerian field + winding-number coordinate recovery** as the dense backbone (not a global mesh). Rationale: block-wise-friendly, gives thickness for free (streamline length), and the winding-number field is the natural cross-brick stitch scalar.

**Where SNIC fits (and doesn't):** SNIC is a *reasonable fragment/representation layer* — keep it to coarsen bricks ~2–3 orders of magnitude before agglomeration, **but its value is only realized if (a) its distance is computed on the sheetness/orientation field, not raw intensity, and (b) its RAG carries signed structure-tensor edge weights plus per-node normal/sheetness, not just adjacency.** Critically, **no leading Vesuvius method uses supervoxels** — they stitch at the level of surface patches / point-clouds / sparse skeleton paths carrying winding numbers. So either attach winding+normal+fiber attributes to your SNIC nodes, or prefer quad-patch / sparse-path nodes. SNIC's global single queue is sequential — run per-brick. **Do not over-invest in the supervoxel algorithm choice; the affinity field dominates outcome.**

**Drop:** Single global passes (SNIC or otherwise) over the whole scroll; vanilla Frangi; mesh-based global flattening (SLIM/BFF/Ricci) as the *engine*; MIQP integer-grid solvers; per-slice 2D optical-flow tracing.

---

## 4. What to avoid / known dead-ends

- **A single global pass at TB scale (any method).** Even Henderson's is global and ~19h on a subvolume — it will not survive 80k³ without domain decomposition. Always brick + halo + stitch.
- **Vanilla Frangi vesselness-as-plateness.** β/c constants are contrast-dependent and won't transfer across bricks; use Beyond-Frangi/strain-energy or Descoteaux/OOF instead.
- **Mesh-based global parametrization (SLIM, BFF, ARAP, CETM/CEPS, Ricci flow, multicut MIQP, field-aligned integer-grid) as the whole-scroll engine.** All per-mesh/in-core; reserve for bounded per-patch rendering only.
- **Unsigned agglomeration (waterz, plain SLIC compactness) as the final partitioner.** No repulsive edges → the entire scroll collapses into one segment once a boundary weakens in a crushed region. This is the catastrophic failure mode.
- **ERS / DISF / SICLE at scale.** No streaming variant; small-ROI quality ceilings only.
- **Integration scale ρ ≥ inter-wrap gap.** Averages two wraps' normals together — silently corrupts the orientation field. Keep ρ below the gap (~1–2 voxels).
- **Over-diffusing with CED.** Bridges adjacent wraps; cap iterations and cross-conductance.
- **Trusting boundary-recall metrics.** Irregular methods game BR ([ill-posed survey](https://arxiv.org/pdf/2411.06478)); use VOI / split-merge (connectomics) and topology-aware scores instead.

---

## 5. Consolidated key-references table

| # | Method / Reference | Class | Link | Year |
|---|---|---|---|---|
| 1 | **Diffeomorphic Spiral Fitting (Henderson)** — *the* end-goal method | hybrid (classical core) | https://arxiv.org/abs/2512.04927 · [code](https://github.com/pmh47/spiral-fitting) | 2026 |
| 2 | Mutex Watershed | classical | https://arxiv.org/abs/1904.12654 | 2019 |
| 3 | GASP (signed-graph agglomeration) | classical | https://arxiv.org/abs/1906.11713 · [CVPR'22](https://openaccess.thecvf.com/content/CVPR2022/papers/Bailoni_GASP_a_Generalized_Framework_for_Agglomerative_Clustering_of_Signed_Graphs_CVPR_2022_paper.pdf) | 2022 |
| 4 | Volara (block-wise MWS at scale) | hybrid (classical core) | https://www.e11.bio/blog/volara · [code](https://github.com/e11bio/volara) | 2025 |
| 5 | Block-wise delayed-merge (Lu/Zlateski/Seung) | classical | https://arxiv.org/abs/2106.10795 | 2021 |
| 6 | Zlateski & Seung quasilinear watershed | classical | https://arxiv.org/abs/1505.00249 | 2015 |
| 7 | SNIC | classical | https://openaccess.thecvf.com/content_cvpr_2017/papers/Achanta_Superpixels_and_Polygons_CVPR_2017_paper.pdf | 2017 |
| 8 | 3D SEEDS (volumetric, no-ML) | classical | https://arxiv.org/abs/2502.02409 · [code](https://github.com/Zch0414/3d_seeds) | 2025 |
| 9 | SSLIC (n-D parallel SLIC) | classical | https://arxiv.org/abs/1806.08741 | 2018 |
| 10 | Descoteaux sheetness (Hessian) | classical | https://www.ncbi.nlm.nih.gov/pubmed/17127650 · [ITK](https://github.com/InsightSoftwareConsortium/LesionSizingToolkit/blob/main/include/itkDescoteauxSheetnessImageFilter.h) | 2006 |
| 11 | Optimally Oriented Flux (OOF) | classical | https://cse.hkust.edu.hk/~achung/eccv08_law_chung.pdf | 2008 |
| 12 | Beyond Frangi (Jerman) | classical | https://ui.adsabs.harvard.edu/abs/2015SPIE.9413E..2AJ/abstract | 2015 |
| 13 | Phase congruency (Kovesi) | classical | https://pubmed.ncbi.nlm.nih.gov/11195306/ · [3D impl](https://github.com/chvillap/phase-congruency-features) | 2000 |
| 14 | Coherence-Enhancing Diffusion (Weickert) | classical | https://link.springer.com/article/10.1023/A:1008009714131 | 1999 |
| 15 | Eulerian harmonic thickness field (Yezzi & Prince) | classical | https://lccv.ece.gatech.edu/docs/tmi_EulerianA_Tissue_Thickness.pdf | 2003 |
| 16 | Heat Method (orientation→scalar Poisson) | classical | https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/ · [scalable](https://arxiv.org/abs/1812.06060) | 2017 |
| 17 | Generalized winding number (Jacobson) | classical | https://igl.ethz.ch/projects/winding-number/robust-inside-outside-segmentation-using-generalized-winding-numbers-siggraph-2013-compressed-jacobson-et-al.pdf | 2013 |
| 18 | Archimedean spiral LSQ fit (Mishra) | classical | https://papers.ssrn.com/sol3/papers.cfm?abstract_id=531542 | 2004 |
| 19 | SLIM / slim-flatboi (per-patch flatten) | classical | https://cims.nyu.edu/gcl/papers/SLIM2017.pdf · [Vesuvius](https://github.com/giorgioangel/slim-flatboi) | 2017 |
| 20 | Higra (C++ RAG/MST/BPT) | classical | https://github.com/higra/Higra | 2024 |
| 21 | Out-of-core BPT (Join/Select/Insert) | classical | https://arxiv.org/abs/2210.02218 | 2022 |
| 22 | PyChunkedGraph (octree RAG) | classical | https://github.com/CAVEconnectome/PyChunkedGraph | 2021 |
| 23 | TeraHAC (trillion-edge HAC) | classical | https://arxiv.org/abs/2308.03578 | 2024 |
| 24 | daisy + LSD blockwise pipeline | classical eng. (hybrid) | https://www.nature.com/articles/s41592-022-01711-z · [daisy](https://github.com/funkelab/daisy) | 2022 |
| 25 | FASP / VC3D (production tracer) | hybrid | https://github.com/hendrikschilling/FASP | 2024 |
| 26 | ThaumatoAnakalyptor (winding + chunk-graph) | hybrid | https://github.com/schillij95/ThaumatoAnakalyptor | 2024 |
| 27 | En-Gedi virtual unwrapping (Seales) | classical | https://www.science.org/doi/10.1126/sciadv.1601247 | 2016 |
| 28 | ACM CSUR superpixel taxonomy (survey) | survey | https://arxiv.org/abs/2409.19179 | 2024 |
| 29 | GALA (borrow features, skip classifier) | ml | https://arxiv.org/pdf/1303.6163 | 2014 |
| 30 | Kaggle Surface Detection (ML bar + metric) | ml | https://www.kaggle.com/competitions/vesuvius-challenge-surface-detection | 2026 |

---

## 6. Open problems & suggested next experiments

**Open problems (the genuine unknowns):**

1. **Single-object topology mismatch.** Every connectomics method segments *many closed objects*; a scroll is *one open sheet wound thousands of times*. There is no published recipe for labeling every *wrap* of a single spiral as a distinct segment via signed-graph agglomeration. The repulsive across-sheet stand-off distance tied to local sheet thickness is the obvious lever but is unproven and thickness varies. **This is taberna's core novel work.**
2. **Affinity quality without learning.** The entire connectomics SOTA assumes high-quality *learned* affinities. Whether structure-tensor planarity/orientation gives a clean enough signed graph in crushed/delaminated/low-contrast regions is **the central unvalidated bet.** No benchmark exists for structure-tensor affinities feeding MWS/GASP.
3. **Inter-wrap sub-resolution separation.** Adjacent wraps are ~1 voxel apart, at/below the scale of any stable derivative filter. OOF mitigates but does not eliminate this — a key reason ML currently wins on the raw scalar.
4. **Chunking-invariant stitching for signed-graph partitioners.** Solved for single-linkage/mean-affinity (Lu/Zlateski/Seung), **not** cleanly for globally-optimal MWS/multicut — their optimality is per-graph, so block decomposition can change the answer.
5. **Joint deformation-field + graph-hierarchy maintenance.** No precedent in connectomics tooling for keeping an undeforming map in sync with re-merges. Novel.
6. **Validation.** No ground truth or accepted metric for per-wrap spiral-segmentation or undeformation accuracy; VOI/split-merge and the Kaggle TopoScore are the closest borrowable metrics, but internal validation must be designed.

**Suggested next experiments (concrete, ordered):**

1. **Signed-affinity feasibility study (highest priority, attacks open problem #2).** On a small annotated ROI from a public scroll, derive signed structure-tensor affinities (in-plane coherence = attractive; normal-direction stand-off = repulsive) and run GASP avg-linkage+cannot-link via the [abailoni/GASP](https://github.com/abailoni/GASP) reference. Measure VOI split/merge vs whether adjacent wraps stay distinct in crushed regions. **This single experiment de-risks the whole bet.**
2. **OOF vs Descoteaux vs phase-symmetry bake-off** for the sheetness scalar on faint vs strong wraps and near-touching wraps. Pick the radius/σ that resolves one wrap without bridging.
3. **Henderson fed classical fields.** Run [pmh47/spiral-fitting](https://github.com/pmh47/spiral-fitting) on a subvolume but replace nnU-Net surface/fiber inputs with your structure-tensor outputs; quantify the gap to the published ML-fed result. This directly tests whether the classical front-end is good enough for the spiral fit.
4. **Delayed-merge stitching prototype** on two adjacent bricks: verify the global label is chunking-invariant per Lu/Zlateski/Seung, using Higra for the in-block RAG/agglomeration.
5. **Eulerian winding-coordinate prototype.** Solve the Laplace across-winding field + Poisson orientation integration on a small multi-wrap region; check monotonicity across wraps and whether the accumulated winding-number field stitches cleanly across a halo (open problem #4 in miniature).
6. **3D SEEDS drop-in benchmark** vs your SNIC as the fragment layer (10× speed claim; blocky-regularity downside on thin sheets).

**Low-confidence flags:** phase-congruency behavior on real scroll CT (2.3) is *medium confidence* — prototype before committing. The Henderson method's technical details and the connectomics graph stack are *high confidence* (adversarially fact-checked). The "structure-tensor affinities are clean enough" assumption underpinning the whole no-ML architecture is **unvalidated** — experiment #1 exists precisely to test it before large engineering investment.
