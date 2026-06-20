# matter-compressor archive format v2 (`.mca` v2)

Design locked 2026-06-20. Implementation target: `third-party/matter-compressor`
(mc is in-scope to edit). Supersedes the v1 static-fixed-slot layout whose logical
size equalled the *uncompressed* volume (slammed the ext4 16 TiB / 48-bit-mmap walls).

## Goals
1. **Compact & dense**: logical = physical = compressed bytes. **No sparse-file dependency** —
   v1 needed holes because static slots made the file logically uncompressed-sized; v2 packs the
   data heap, so the file is naturally compact and is a normal dense file (portable, predictable
   size, no `cp --sparse` / filesystem-hole reliance).
2. **Streamable**: lock-free random-position shard writes (many uncoordinated workers fetching from a remote source), *and* a single-writer full export — same file format.
3. **Variable-rate, noise-floor quantization** per shard.
4. **Uniform 16-bit addressing** in all three axes.

---

## Implementation status (2026-06-20): RESTORED & VERIFIED
v2 is **the appendable archive that existed before the static-slot rewrite** (`mc_archive @ 519b5b1`),
surgically restored: `mc_archive.{c,h}`, `mc_stream.c`, `mc_volume.c`, `mc_codec.{c,h}` restored to
519b5b1 + the NaN-guard re-applied; the `295b73b` fixes to cache/solve/trace/zarr kept. **Builds clean**
(lib + the viewer consumer); **round-trips** a real 256³ (q=32 → MAE 8.1, air preserved, ~41× on a
1-shard archive). The restored code already realizes this design; sections below are reconciled to it.
A couple of real deltas from the paper design — both improvements over what I'd sketched:
- **Index is a 3-level dense node tree** (root→inner→shard, each a flat 4096-u64 child-offset array
  indexed by the chunk-coord nibble; slot 0 = absent), updatable in place + crash-safe — NOT a flat
  dense `{offset,len}` table. It's appendable, sparse (only materialized nodes exist), and covers
  **2^12 chunks/axis = 2^20 voxels/axis** (the 2^16 cube fits with room to spare).
- **The 2^16 rotation goes in the user-metadata region** (`mc_archive_set_metadata`, 128 KB carve-out),
  not a new header field — no format change. The export pipeline does the resample + writes the matrix.

## 1. Coordinate frame — uniform 2^16 cube
Some scrolls are long in one axis (real content 68k–76k voxels — exceeds 2^16). **Rotation is
conditional**: only when the longest axis > 2^16, rotate so the umbilicus (scroll axis) lies on
the **(1,1,1) body diagonal** — the empty cylinder corners then let all three axes fit a
**uniform 2^16 cube** (measured: PHercParis4 74272→~58k, PHercParis3 68128→~62k). Volumes that
already fit are stored axis-aligned (e.g. **PHerc0332 is 15761×15761×33592 — z=33592 < 65536,
no rotation needed**). The header flags whether a rotation was applied.

- Per-axis coordinate: **16 bits**. Morton key = 3×16 = 48 bits + 3 LOD bits = 51 bits (u64).
- The rotation is folded into the existing preprocessing resample (one interpolation; the
  data is already heavily interpolated from CT reconstruction, so this is free in context).
- **Header stores the 3×3 rotation matrix + origin (f64)** mapping archive→world/source
  coords, so coordinates stay recoverable and the unwrap pipeline knows the umbilicus dir.
- Cost accepted: ~3× lower occupancy (cheap — extra ALL_ZERO shards at 2 bits each).

## 2. Hierarchy (per LOD; **LODs are fully independent** — own grid/index/data, separate fetch/encode/decode)
| level | size | role |
|---|---|---|
| **DCT block** | 16³ voxels | transform + quant unit; self-contained (own air mask) |
| **shard** | 16³ blocks = **256³ voxels** (4096 blocks) | network-IO / fetch / write-atomicity unit; a Z-order rope of blocks |
| **shard grid** | up to 2^8 = 256 shards/axis at LOD0 (halves per LOD) | 2^24 shards max for the full 2^16 cube |
| **LODs** | 8, independent | LOD k coarser ratio (halve target per level, ~+17% pyramid) |

## 3. File layout (per LOD region) — dense, no holes
```
[header]
[shard occupancy : 2 bits/shard, dense, allocated up front]
[shard offset table : u64/shard, dense, allocated up front, computed slot]
[shard data heap : packed REAL shards, append-grown]
```
Header: magic, version, per-LOD dims, rotation flag + matrix + origin, codec params, region offsets.
The index regions are **sized to the volume's actual shard dims** (not the full 2^16 cube) and
allocated outright; the heap grows by append as shards are written. No sparse file.

**Access model.** The runtime **mmaps the entire archive** and randomly accesses compressed data:
read a shard's `{offset,len}` from the computed offset-table slot, then the shard's bytes are a
**single contiguous region in the mapping** — decoded directly from the live mmap, no read()/gather.
Because each shard is contiguous, it is also the **base unit of network IO**: a `DONT_KNOW` shard is
fetched from the remote source in **one contiguous range request**, appended to the heap, and its
offset published — after which it's live in the mmap. This is v1's mmap model, but it now scales:
v2 is compact (~21 GB dense), so the whole-file mapping is trivial (vs v1's 9.75 TB mapping that hit
the 48-bit VA ceiling).

## 4. Shard index
- **Occupancy — 2 bits/shard, 3-state**, dense, addressed by shard Morton index:
  `DONT_KNOW` (not fetched → ask remote source) · `ALL_ZERO` (fetched, all air → skip, no data) · `REAL` (has data).
  This map *is the streaming manifest*: a client pulls it up front and knows what to skip vs fetch.
  Sized to actual shard dims (e.g. PHerc0332 = 62×62×132 shards → ~127 KB; ~4 MiB only at the full 2^16 cube).
- **Offset table — dense `{offset:u64, len:u32}` per shard at a *computed* slot**
  (`offset_table_base + morton(shard)*16`), pointing into the data heap. The **length** is stored so a
  reader fetches the whole shard in **one contiguous I/O** `[offset, offset+len)` (one HTTP range request) —
  no header round-trip, no scatter-gather. This computed slot is also what makes **lock-free random writes**
  work (slot is computed, not rank-derived). Absent/non-REAL → 0. Allocated outright, sized to actual dims
  (PHerc0332 ≈ 507 K shards × 16 ≈ **8 MB**; ≤256 MB only at the full 2^16 cube) — ~0.04% of a ~21 GB archive.
  `DONT_KNOW` vs `ALL_ZERO` is disambiguated by the occupancy 2-bit (offset 0 alone can't).

## 5. Shard body (each REAL shard = a contiguous chunk-blob in the heap)
Actual restored blob layout (`encode_chunk_blob`), all contiguous:
```
[f32 q][u64 xxh64][u16 fmaplen][fracmap (fmaplen B)][bitmask 512 B][u16 lengths × npresent][block payloads]
```
- **q is per-blob (f32)** — the format supports per-chunk q; we just write the *global* q to every blob.
- **xxh64** over (fracmap|bitmask|lengths|payloads) → crash-safe verify-on-decode.
- **fracmap**: context-coded per-block material fractions (for the ML sampler / air detection).
- **bitmask** (512 B, 4096 bits): block present (REAL) vs air. **Lengths are u16, present-blocks only**;
  **offsets are implicit = prefix-sum** of the lengths (no offset table, no exceptions, no checkpoints).
- **Locate block k**: `bit k` of bitmask (air→zeros); `r = rank(bitmask, k)`; `offset = Σ lengths[0..r-1]`.
- u16 always fits (a 16³ block ≤ ~4 KB). Overhead ≈ 512 B + fracmap + 2·npresent ≈ **~1–3%** at q≈32.
  (A later 1-byte-length + checkpoint refinement could trim this; u16 + prefix-sum is the proven form.)
- (No 64 KB page layer — the chunk-blob *is* the variable-size IO unit, ~200 KB–1.3 MB at q≈32.)

## 6. Block codec (16³ DCT) — unchanged from mc
Mask-aware: own context-coded 16³ air mask; **harmonic air-fill** (red-black SOR, the
energy-minimizing fill — confirmed optimal) before integer DCT-16 + dead-zone quant + CABAC;
air force-zeroed on decode. Per-block self-contained / independently decodable.

## 7. Quantization — **global q**, `tau = 2q`
- **One global q for the archive** (stored in the header, like v1's `MCAH_QUALITY`). The caller (e.g.
  the fysics export pipeline) chooses it — picking quality is a user/export decision, not the format's.
  **No per-shard q, no per-chunk noise-floor estimation** in the codec.
- `tau = 2q` (tau is the binding cost — archive size ∝ tau; tau=q is wasteful).
- Recommended default (informed by the noise-floor analysis on this data): **q≈32 / tau=64 (~60×)** —
  added error is noise-level throughout (p50 0.45σ, p90 1.2σ, p99 2.2σ, max 4.8σ). Use q=16 (43×) or
  q=8 (27×, edge-safe) if more edge fidelity is wanted. Whether q=32 costs sheet-edge detail is the
  export pipeline's call via its choice of q — not something the format decides.

## 8. Write modes
- **Full export** (single writer): iterate shards in Morton order, append each to the heap, fill
  occupancy + offset. Optimal read locality, no fragmentation. This is also the **compaction** pass.
- **Streaming** (many uncoordinated writers): each fetches a `DONT_KNOW` shard from source, encodes,
  **atomic-bump-appends** its bytes to the heap, then publishes `{occupancy=REAL, offset}` (offset
  written last = crash-safe; a torn write reads back as absent → re-fetch). Lock-free, random-position.
  Re-fetch/overwrite orphans old heap bytes → reclaimed by periodic compaction.

## 9. Sizing (measured on PHerc0332)
- q=32/tau=64 ≈ 60× ⇒ full archive **~21 GB** (vs ~103 GB at v1 q=3 — **~5× smaller**; the savings
  are sub-noise-floor precision v1 was wasting on noise).
- Index allocated outright, sized to actual dims (PHerc0332: offset table ~4 MB + occupancy ~127 KB ≈ 0.02%); intra-shard metadata ≈ 1%. No sparse file.
