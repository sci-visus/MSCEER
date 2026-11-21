# Out-of-Core Path to 4K³ MicroCT: Instructions (Track C)

Stand-alone direction-setting guide. Written against branch `cuda-gradient`, HEAD
`598edd2`. This track is staged farther out than Tracks A/B
(docs/forest_1skeleton_instructions.md, docs/gpu_setdim_instructions.md); C1 is
implementable any time, C2 is the master plan's 3D kernel stage unchanged, C4 is the one
research-grade design in the whole roadmap. Re-verify line references before editing.

---

## Goal

A credible, staged path from today's fits-in-VRAM 2D pipeline to discrete-gradient +
extremum-network analysis of ~4096³ microCT volumes — beyond both GPU memory (8 GB) and
host RAM. Not implemented today; each stage below is a committable brief with a gate, so
the track can be picked up piecemeal.

## The design law (derived from measured numbers — internalize before designing anything)

At 4096³: 68.7 G voxels; doubled lattice 8191³ = **549.4 G cells**.

| Quantity | Size / rate | Consequence |
|---|---|---|
| Gradient bytes (1 B/cell) | 549 GB | cannot be device- OR host-resident |
| Function (4 B/voxel) | 275 GB | same |
| Today's CPU 3D path (grad + max/min label bytes) | ~1.92 TB | dead; label fusion (proven in 2D Stage 2) drops it to ~825 GB — still dead |
| int32 flow-terminal labels | 68.7 G > 2³¹ | **overflow ×32** — int64 or (brick, local) ids everywhere |
| int64 dense terminal map | 550 GB/direction | never materialize a dense global terminal map |
| Fused kernel throughput (measured 338 GB/s DRAM, ~60-70 ms/Gvoxel) | 4K³ ≈ **~4 s** | compute is free |
| Full-gradient D2H at pinned ~20 GB/s | **~7.6 hours** | dense outputs are the enemy; the architecture is I/O shaping |
| `msc_3d_lib` bootstrap (float per topological cell) | 2.2 TB | the 3D facade must be bypassed at scale, not adapted |

**Law: dense global anything is dead. Outputs must be compacted** — extremum network +
merge forest + optional 1-skeleton polylines are MB-scale at any volume size. The merge
forest (`include/gi_extremum_merge_forest.h`) is exactly the right output shape for this
scale, which is why Track A completes it first.

Host-side law too: the master plan's Sub-plan 6 gate (1280³ at constant *device*
footprint) assumes the volume fits in host RAM. 4K³ does not ⇒ **disk-brick out-of-core**,
not just VRAM streaming.

## Existing foundations (verified in-repo; build on these, don't reinvent)

- **`steepest_lstar_tiled/localized_morse_complex.h`** — the only per-block gradient with
  halos in the repo. 64³-vertex tiles; 1-low ghost / 2-high (bridge+ghost) ownership
  (`:121-125`) that partitions the global cell set **exactly** (proven by the stitched
  gradient-file writer `save_discrete_gradient_to_file1` `:1063`); halos are re-read from
  the input function — **no gradient exchange needed**. The memory-lean variant
  `compute_localized_morse_complex_from_grid1` (`:733`) keeps only
  `tile1{position, size, gradient_size, GradBitfield*}` (~2.1 MB per tile) and frees all
  scratch. The `:3` TODO already states the goal: "build directly from a file so the
  global data does not need to be in memory." Known bug to fix in passing: cross-tile
  neighbor indexing at `:661-663` (the dead `query_descending_manifold*` traversals are
  not salvageable; ignore them).
- **`TopologicalSparseBlockedMesh`** (`include/gi_topological_regular_grid.h:3072`) —
  blocked doubled-lattice topological mesh whose blocks match base-grid blocks 1:1, with
  8 submesh archetypes; built on `RegularGrid3DDecompositionImplicit/Fixed`
  (`include/gi_regular_grid.h:1248/:1816`) with packed `(block << 32 | local)` ids and
  block-local 6/27-neighborhood iterators.
- **Sparse / lazy storage**: `SparseBlockedLabeling` (`include/gi_labeling.h:252`, per-
  block lazy allocation + `has_block`), `IndirectLabeling` (`:699`, compacted-domain
  gradient storage already used in production by `sparsesteepest` — the masked-domain
  option for microCT with large air fraction).
- **Disk layer**: `GetSubBlockRaw` / `GetSubBlockIdx`
  (`sparseblockedcompute/sparseblockedcompute.cxx:183/:227`) — seek-read a sub-box from a
  raw file / OpenVisus IDX dataset.
- **Partitioned merge machinery**: `ExtremumMergeForest::Build` partitioned mode
  (frozen/taint deferral + serial residue) and the partitioned-MSC freeze/lineage
  pattern — the in-repo precedents for cross-block topology resolution.
- **2D kernel structure**: the fused kernel's halo-staging tile
  (`gpu_dgrad/cls/cls2d_kernel.cu:285-303`) generalizes to per-brick 3D launches;
  `LabelCtx2D` is the persistent-device-resource pattern.

**Missing and load-bearing (each is a stage below):** any CUDA stream / `cudaMemcpyAsync`
/ pinned buffer in a compute path (zero hits repo-wide; pinned appears only in the
bandwidth smoke test); any 3D device kernel (`ComputeDiscreteGradient3D` is a stub); a
cross-brick flow-terminal convergence design (master plan Sub-plan 7b names it, no
design); any notion that host RAM may not fit.

---

## Stage C1 — Pinned + streams transfer infrastructure (2D now, 3D-ready)

Small utility in `gpu_dgrad/`: pinned buffer pool, 2-3 `cudaStream_t`,
`cudaMemcpyAsync`, event-based buffer recycling. Convert the batched fused path and the
labeling/ctx entry points. **No algorithm change.** Measured rates: pageable ~11/9 GB/s
vs pinned ~20.9/20.1 (the link trains Gen5 ×8) — the plan's ~21 GB/s number is currently
realized by no real path.

**Gate**: byte-identical outputs on the full validation set; measured in-pipeline H2D/D2H
at ~21/20 GB/s recorded in `gpu_dgrad/docs/README.md`; batched-slice throughput improvement
reported.
**Files**: new `gpu_dgrad/common/gpu_transfer.h/.cu`, `cls2d_kernel.cu`, `labeling2d.cu`.

## Stage C2 — 3D fused kernel

Master plan Sub-plan 5 (`~/.claude/plans/look-at-the-discrete-iterative-map.md`),
**unchanged**: 27-slot state, `(B+2)³` shared staging, per-boundary-class sub-expansions,
tables generated from `Explicit3x3x3SmallRegularGrid`. Gate = full-array memcmp vs
`SlidingWindowRobinsNoalloc::ComputePairing_sliding` on ≥3 volumes (256³/512³ CPU
baselines) + Euler/critical-count checks. Nothing in this doc alters it; C3+ consumes its
kernel. Note for later: the 2D offset-doubling successor generalizes (vertex→edge→vertex
and hex→quad→hex are ±2 steps; middle-dimension flows branch and need the general form).

## Stage C3 — Brick container unification (Sub-plan 6 + host-exceeding extension)

**Spine**: *tile = GPU streaming unit = disk brick = MSC/forest partition.* Build a
`BlockedGradientField` from the three foundations: tile1 ownership (64³ vertices,
1-low/2-high halo, halos re-read from input), `TopologicalSparseBlockedMesh` addressing,
`GetSubBlockRaw`/`GetSubBlockIdx` I/O.

**Budget per brick**: function tile ~66³ floats ≈ 1.1 MB; gradient 129³ bytes ≈ 2.1 MB;
triple-buffered device slots ≈ 10 MB total ⇒ device footprint is never the constraint.
**I/O model**: input bricks read ≥ once (halo re-reads amortized by brick-row scheduling);
gradient bricks written exactly once (or kept only transiently if downstream consumes
them immediately); nothing else dense.

**Extension beyond Sub-plan 6**: bricks stream disk → pinned → device → pinned → disk
with a bounded LRU host brick cache; **no allocation proportional to volume anywhere**.
Fix `localized_morse_complex.h:661-663` in passing.

**Gate**: streamed == monolithic == CPU byte-identical at ≤ 768³; 1280³ at constant
DEVICE footprint (Sub-plan 6's gate); NEW: **2048³ synthetic with the host cache capped
far below volume size completes at constant HOST footprint**.

## Stage C4 — Cross-brick flow terminals + extremum forest, compacted-only (research stage)

The one genuinely novel design. Goal: the extremum network + merge forest of a volume no
single memory can hold, with **no dense global terminal map ever existing**.

**Per-brick labeling**: run 3D pointer doubling *inside* the brick with successors
clamped at brick faces. The dense int32 label array is transient and device-side
(129³ × 4 B ≈ 8.6 MB), discarded after compaction. Persistent outputs per brick:
- *local terminal records*: critical cells owned by the brick (int64 global ids, or
  (brick, local) pairs);
- *frontier records*: one per face-crossing V-path — `(entry cell or interior origin) →
  exit successor global id`. An entry→exit map, NOT per-lattice-point.

**Global stitch**: the frontier records form a successor function on a small set (only
face-crossing paths). **Pointer-double THAT** (host or device) — log passes resolve every
frontier to a true critical terminal even when flat regions chain through many bricks.
Downstream only ever needs terminal lookups *at critical-cell neighbors* (arc
extraction), which are brick-local + frontier-resolved.

**Extremum network + forest**: per-brick arc lists (2-cell/1-cell saddle → terminal
lookups, the 3D analogue of `ExtractArcsBanded2D`) with cross-brick arcs frozen — feed
`ExtremumMergeForest::Build` partitioned mode **as-is** (it is dimension-agnostic; the
frozen/taint machinery is the precedent). Mitigations if the global residue is too
large: hierarchical (brick → octree-level) merges; brick-level pre-collapse of sub-noise
extrema (local persistence pre-merge) before globalizing — also the mitigation for forest
node-count growth on noisy microCT.

**Gate**: on a 512³ real volume, the brick pipeline's terminals-at-critical-cells, arc
set, and forest `(parent, pers, mergeSaddle)` records match the in-core reference (the
C5 milestone pipeline) exactly modulo documented tie order; a peak-memory audit shows
nothing scaling with global voxel count except brick I/O buffers.

## Stage C5 — 3D milestones and the 4K³ gate

- **First concrete milestone when 3D work starts** (build BEFORE C4 — it is C4's
  reference): a 512³-768³ single-GPU **resident** forest pipeline — C2 fused gradient on
  device, device 3D pointer doubling, compacted arc extraction D2H, host forest build,
  and the Track-A-style skeleton trace (this is where the deferred GPU two-pass trace
  from A3 lands, per brick with exit records). Proves the 3D MSC-free path end-to-end at
  in-core scale.
- **4K³ gate**: 4096³ (procedurally generated bricks acceptable) streamed end-to-end at
  constant host+device footprint (target ≤ 6 GB device, ≤ 32 GB host); outputs = extremum
  network + merge forest + optional 1-skeleton polylines, all compacted; optional
  gradient bricks on disk; wall clock reported against the calibrated model (kernel
  ~4 s; total dominated by brick I/O).

---

## Risks

- **Frontier-record spikes**: flat/plateau regions crossing many bricks lengthen chains
  and multiply face crossings. Entry→exit compaction bounds records per brick face;
  frontier-graph pointer doubling bounds resolution passes at log; plateau pre-resolution
  inside bricks if measurements demand.
- **id discipline**: int32 overflows ×32 at 4K³ — int64 or packed (brick, local) ids in
  every frontier/terminal/arc record from day one. Forest compact node ids may also
  approach 2³¹ on noisy microCT; brick-level pre-collapse bounds the node count.
- **Sub-plan 6's host-RAM assumption**: its device-only gate silently permits
  volume-proportional host allocations; the C3 constant-host-footprint gate is what
  actually points at 4K³ — do not skip it.
- **3D middle-dimension flows branch** (unlike 2D): the quad/edge saddle traces and any
  3D skeleton work need the general successor form, not just offset doubling. The
  extremum network only needs the extreme-dimension flows (vertex and hex), which do NOT
  branch — that is why the forest path scales and the full MSC path does not.
- **Masked domains**: microCT with large air fractions should route through
  `IndirectLabeling`-style compaction (precedent in `sparsesteepest`) to cut the 549 GB
  proportionally to occupancy — worth a design note in C3, not a separate stage.
