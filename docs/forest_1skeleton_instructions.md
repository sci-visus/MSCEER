# Forest-Native 1-Skeleton: Instructions (Track A)

Stand-alone implementation guide. Written against branch `cuda-gradient`, HEAD `598edd2`
("Merge-forest simplification mode: MSC-free 2D segmentation"). Re-verify cited line
numbers against your HEAD before editing; the anchors are function names, not lines.

**Recommended execution order across tracks: A1 → B1 (docs/gpu_setdim_instructions.md)
→ A2 → A3 memo → C1+ (docs/outofcore_4k_microct_instructions.md).** A1 is hours and
gates everything else in this doc; A2 is the headline capability and is independent of
Track B (parallelizable with it).

---

## Goal

`Simplification::MergeForest` mode produces persistence-simplified region segmentations
without ever building the Morse-Smale complex — but it currently discards every saddle,
every arc, and all geometry. Any request for `criticalPoints()`, `arcGeometry()`, or
`computePolylineGraph()` falls back to lazily building the full MSC on the CPU
(`libs/msc_2d_lib/msc_2d_lib.cxx:461-469`), which costs `EnsureAscendingManifoldDimensions()`
(~580 ms at 3232²) plus MSC construction — the exact CPU path the forest was built to avoid.

After this track, forest mode serves geometric ridge/valley lines (the 1-skeleton),
critical points, and the polyline graph **natively from the forest + gradient bytes**,
with the MSC fallback demoted to a diagnostic switch. A3 additionally shapes the tracing
code so a GPU trace kernel can be dropped in later without API change.

## Verified context

### What the forest stores and what it destroys

`include/gi_extremum_merge_forest.h`, class `ExtremumMergeForest`:

- Per node (`:325-337`): `mNodeCell` (INDEX_TYPE, the extremum's cell id — the output
  identity), `mNodeValue` (float), `mNodeTieKey` (INDEX_TYPE, `Cell2HighestVertex`),
  `mParent` (int32, survivor at death, self for roots), `mMergePers` (float, +inf for
  roots), plus build/memo scratch.
- `Arc { e1, e2, saddleValue, saddleCell, pers }` — **consumed by `Build()`**: moved into
  a per-partition `priority_queue` (`:132-134`) and the input vectors cleared. Nothing
  arc-shaped survives.
- `applyMerge(r1, r2, const Arc& a)` (`:310-323`) records `mParent[rd] = rs` and
  `mMergePers[rd] = a.pers` — **and drops `a.saddleCell`, which it has in hand.** That is
  the single missing datum for geometry.
- CRITICAL NUANCE (verified): by pop time, a stale arc's `e1/e2` have been lazily
  rewritten to *current representatives* (`:144-151`), so the arc's original **base**
  endpoints are unrecoverable from the Arc. Do NOT try to store them in A1. They are
  recoverable later from the retained terminal maps: the saddle cell's two facet
  vertices (min side) / two cofacet quads (max side) looked up in `forestTermLat` give
  the base endpoints. The saddle cell alone is sufficient.
- `BuildRemap(t, outRep)` (`:221-240`): memoized parent-chain climb while
  `mMergePers <= t`; O(n) amortized; NOT thread-safe (min/max forests are separate
  objects, so directions still run concurrently).

### The V-path facts that make MSC-free geometry possible

- In 2D the descending V-path from a saddle side is **deterministic and non-branching**:
  from critical edge `s`, each of its 2 facet vertices starts a walk where
  `succ(c) = 2*getPair(c) − c` (offset doubling); symmetric for the 2 cofacet quads
  ascending to maxima. Needs only `getCritical` + `getPair` — **no ldir**, so it works in
  forest mode where `setAscendingManifoldDimensions` never ran.
- Typical path length ~10-30 cells (measured 232 arcs / 6738 points ≈ 29 cells/arc,
  2 cell ids per gradient step); scales with feature size, not image size.
- `ExtNet2D::CpuFlowTerminals2D` (`include/gi_extremum_arc_extract_2d.h:77-118`) already
  walks exactly these chains with memoization — but discards the ordering (stamps the
  terminal only). `ExtractArcsBanded2D` (`:262-379`) already enumerates critical edges
  and does the two terminal lookups per side, skipping boundary edges on the max side
  (single-cofacet edges), exactly as `trace_up_2saddle` did.

### The MSC machinery you must match (for gating)

`include/gi_morse_smale_complex_basic.h`:
- Base arc geometry = the DFS stack of `rec_tdcr` (`:523-549`): alternating cell sequence
  `[upper crit cell, facet, paired cell, facet, …, lower crit cell]`, ordered
  upper → lower, ~2 cells per gradient step, endpoints included.
- Merged arc geometry = recursive triples `arc_merged_geometry{lower-arc, cancelled-arc,
  upper-arc}` resolved by `fillArcGeometry` (`:1584` → `_fillArcGeometry:1477`) with a
  2-deep backtrack canceller at joints. Orientation invariant: always upper → lower.
- `arcGeometry()` (`msc_2d_lib.cxx:1326-1363`) emits world coords = `cellid2Coords * 0.5`
  (half-integer coordinates on edge midpoints); `computePolylineGraph()` (`:1364-1376`)
  dumps arc cells as a soup into `MeshCellsGraphBuilder` (`include/gi_graphs.h:359-492`),
  which re-derives junctions by degree and smooths with
  `ConstrainedLine::IrreversableSmooth(10)` (clamped ±0.5/axis).

---

## Stage A1 — Record the merge saddle

**Objective**: stop destroying the one datum that lets the forest answer geometric
questions. ~15 lines + a validator. Do this first: every forest built afterwards
(including by unrelated work) then carries the data.

**Changes** (`include/gi_extremum_merge_forest.h`):
1. Member `std::vector<INDEX_TYPE> mMergeSaddle;` next to `mParent`/`mMergePers`
   (`:330-331`).
2. In `Build()` init block (`:97-104`): `mMergeSaddle.assign(mN, (INDEX_TYPE)-1);`.
3. In `applyMerge` (`:310-323`): `mMergeSaddle[(size_t)rd] = a.saddleCell;`.
4. Accessor: `INDEX_TYPE NodeMergeSaddle(int32_t id) const { return mMergeSaddle[(size_t)id]; }`
   (-1 for roots), next to `NodeCell/NodePers/NodeParent` (`:242-245`).

Cost: 8 B/node — noise (nodes are extrema, ~10⁵ at 3232²).

**Gate** (add to `libs/msc_2d_lib/msc_2d_lib_smoke.cxx`, run on the standard case matrix):
1. *Structural*: every non-root's `mMergeSaddle` is a critical edge
   (`grad->getCritical(s)` and `mesh->dimension(s) == 1`); every root is -1; the
   saddle's two terminal extrema (facet vertices → asc terminal map for the min forest;
   cofacet quads → dsc terminal map for the max forest) resolve to valid compact node ids.
2. *Tie-safe replay*: sort non-root records by `(mMergePers, mMergeSaddle)` (the heap
   order), replay a **fresh** union-find; at each record `(d, s = mParent[d], saddle)`,
   require `{find(nodeOf(termA)), find(nodeOf(termB))} == {d, find(s)}` **before**
   applying the union. This states "representatives just below the merge" without an
   epsilon (quantized data makes epsilon-below ill-defined; replay order is ground truth).
3. *Parity*: serial (`numPartitions = 1`) and partitioned/banded builds produce identical
   `mMergeSaddle` arrays (extends the existing monotone + frozen-closure gates).

**Files**: `include/gi_extremum_merge_forest.h`, `libs/msc_2d_lib/msc_2d_lib_smoke.cxx`.

---

## Stage A2 — Geometric 1-skeleton from the forest (no MSC)

**Objective**: forest mode serves `arcGeometry()`, `criticalPoints()`, and
`computePolylineGraph()` natively. Keep the lazy-MSC path behind a diagnostic switch
(needed for the gate; also the escape hatch if a consumer needs true MSC semantics).

### Definitions (per direction; ridge = max/descending forest, valley = min forest)

- *Living extremum at t*: node `i` with `mMergePers[i] > t` (equivalently
  `BuildRemap(t)` maps it to itself).
- *Living saddle at t*: a critical edge `s` **not** in
  `{ mMergeSaddle[d] : mMergePers[d] <= t }` for this direction's forest. Its arc
  endpoints are `Rep(nodeOf(termA), t)` and `Rep(nodeOf(termB), t)`. Equal endpoints =
  loop arc — keep it (real ring geometry; the MSC keeps loop arcs too).
- Optional stricter mode (closer to the MSC's interleaved min/max cancellation): also
  drop saddles consumed by the *opposite* forest at ≤ t. Ship both; default =
  per-direction. Document the difference.

### Base V-path trace

`BasePath(s, side)`: from critical edge `s`, take the side-k cofacet quad (ridge) or
facet vertex (valley); walk `cur → pair = getPair(cur) → succ = 2*pair − cur`, recording
**every cell** (the d-cell and the pair edge, matching the MSC's alternating sequence),
until `getCritical`; append the terminal. Lock-free `omp parallel for` over saddle-sides,
each writing its own vector — no critical section (unlike the MSC's `createArc` lock).
Deliberately ldir-free.

Cache per (saddle, side): computed once on first request, reused across thresholds and
shared by many stitched lines.

### Stitching — the link tree (the one genuinely new algorithm)

Naive parent-chain concatenation is **geometrically wrong**: when node `d` dies via
saddle `s₂`, the base endpoints of `s₂` are two base extrema `(bA, bB)` — one in `d`'s
subtree, one in the survivor's — and `bA` is generally NOT the base extremum your
polyline currently ends at. The connected path from your endpoint to `bA` runs through
*earlier* merge saddles of the subtree.

Build once per direction the **link tree** `L`:
- Nodes = base extrema (compact ids).
- For every non-root `d`: one edge `(bA_d, bB_d)` labeled
  `(saddle = mMergeSaddle[d], pers = mMergePers[d])`, where `bA/bB` are the saddle's two
  base terminal endpoints from the retained terminal maps.
- Each merge joins two previously-distinct components ⇒ `L` restricted to `pers <= t` is
  a forest whose components are exactly the t-classes, and the unique tree path between
  two same-class extrema uses only edges with `pers <= t` (the subforest path IS the tree
  path).
- Root each component; store `depth`, `linkParent`, per-edge `(saddle, pers)`. Path query
  = LCA walk, O(path length).

`RidgeLine(s, side, t)` = `BasePath(s, side)` ending at base extremum `b0`, then for each
link edge on the tree path `b0 → Rep(b0, t)`: append that edge's saddle's two `BasePath`s
(one reversed) oriented to continue the walk. Where forest merges coincide with MSC
cancellations this reproduces the MSC's `{lower-arc, cancelled-arc, upper-arc}` recursion
exactly, including the joint duplication the 2-deep backtrack canceller removes — apply
the same canceller at joints.

### API mapping (`msc_2d_lib`)

- `arcGeometry()`: fill `ArcGeometry{id, lowerNodeId, upperNodeId, dim, line}`; world
  coords = cell coords × 0.5 (existing convention). Forest-native id space: saddle ids =
  enumeration of living critical edges; extremum ids = compact node ids. Document that
  ids differ from MscHierarchy-mode MSC node ids.
- `criticalPoints()`: living extrema from the forest + living saddles from the scan;
  values via `meshfunc->cellValue`. No MSC.
- `computePolylineGraph(useValleys)`: feed stitched polyline cells into the existing
  `MeshCellsGraphBuilder` → `BuildGeometricGraphFromMeshGraph(…, 10)` path unchanged.

**Gate** (case matrix; MSC comparisons use the diagnostic switch):
1. **t = 0 exact**: every saddle's two-sided polyline matches the MSC base
   `fillArcGeometry` cell sequence memcmp-exactly (same V-paths; exact match expected).
2. **t > 0 self-consistency (exact)**: every emitted polyline starts at a living saddle,
   ends at the cell of the living representative per `BuildRemap(t)`; consecutive cells
   are facet/cofacet-adjacent; every merge saddle traversed has `pers <= t`.
3. **t > 0 vs MSC (statistical, frozen)**: forest merges are not MSC cancellations
   (documented ~99.7-99.9% region agreement), so simplified arc sets may legitimately
   differ. Metric = % of living MSC arcs at t whose saddle is forest-living and whose
   polyline matches; measure once, FREEZE as a regression floor. Not 100%, by design.

**Files**: new `include/gi_extremum_skeleton_2d.h` (trace + link tree + stitch; header-only
like its siblings), `libs/msc_2d_lib/msc_2d_lib.cxx` + `.h` (forest-native dispatch +
diagnostic switch), `libs/msc_2d_lib/msc_2d_lib_smoke.cxx` (gates).

---

## Stage A3 — Tracing framework: keep CPU and GPU paths open (NO implementation)

Per project decision: do not implement GPU tracing now. Instead, shape A2's code so a GPU
tracer drops in later, and record the decision data.

**Interface contract** (put in `gi_extremum_skeleton_2d.h`):
- Input: gradient bytes + a list of (saddle cell, side) requests.
- Output: CSR-style — `offsets[i]` into one flat `cells[]` array of INDEX_TYPE polyline
  cell ids, order = request order, each path upper→lower with endpoints included.
- The CPU walker (A2) and any future GPU tracer must both produce exactly this; the
  stitcher consumes only this contract. Keep the trace behind a single free function
  (e.g. `TraceBasePaths2D(grad, requests, offsets, cells)`) so the backend is swappable.

**GPU design to document (not build)** — two-pass, on a device-resident gradient
(`GradCtx2D` from Track B is the landing zone):
- Pass 1: one thread per (saddle, side) walks `succ = 2*pair − c` counting steps.
- Exclusive scan of counts → offsets.
- Pass 2: re-walk, emitting cell ids into the compacted buffer; D2H downloads only the
  compacted skeleton. Warp divergence on unequal path lengths is acceptable
  (bandwidth-bound gather).

**Why deferred**: after simplification, living saddles are few (hundreds-thousands) with
~10-30-cell paths — the parallel CPU walk is sub-10 ms and kernel launch + D2H overhead
loses. In 2D the gradient is host-resident anyway for MSC/file consumers, so GPU tracing
saves nothing end-to-end. The GPU case is (a) the t=0 base skeleton of ALL saddles at
large sizes and (b) 3D, where compacted-output-only is mandatory (see Track C, stages
C4/C5 — the two-pass design ships there).

**Deliverable**: a short section in `gpu_dgrad/README.md` after A2 lands: measured
base-skeleton trace time at 1024²/3232², and the revisit trigger — *base-skeleton trace
exceeding ~20% of total forest build time, or the move to 3D*.

---

## Risks

- **Forest-vs-MSC divergence surfaces in geometry gates**: only t=0 can be exact.
  Mitigation: the exact self-consistency gate + the frozen statistical baseline; offer
  the stricter both-directions living-saddle mode.
- **Stitched polylines can double back** through merge saddles where forest merge order
  diverges from MSC cancellation order. Correct by construction, aesthetically different;
  `IrreversableSmooth` hides most of it. Document, don't fight it.
- **Terminal-map retention**: A2 needs `forestTermLat` kept alive after `Build()` (it
  already is, for painting); if a future memory optimization drops it, the link tree
  loses its base endpoints. Note the dependency where the maps are stored.
- **Loop arcs and boundary saddles**: max-side boundary saddles have one cofacet (skip,
  as `ExtractArcsBanded2D` does); equal-representative saddles are loops — keep, and make
  the gate count them explicitly so they aren't silently dropped.
