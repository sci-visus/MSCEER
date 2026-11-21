# GPU Dim-of-Ascending-Manifold Pass: Instructions (Track B)

Stand-alone implementation guide. Written against branch `cuda-gradient`, HEAD `598edd2`.
Re-verify cited line numbers against your HEAD before editing; anchors are function names.

**Position in the roadmap**: B1 comes right after A1 (see
docs/forest_1skeleton_instructions.md). It is self-contained, has a hard memcmp gate, and
its residency refactor (`GradCtx2D`) is the landing zone for the future GPU trace (A3)
and the pinned/streams work (C1 in docs/outofcore_4k_microct_instructions.md).

---

## Goal

Replace the ~580 ms (at 3232²) CPU `setAscendingManifoldDimensions` pass with ~10-20 ms
of kernels running on the gradient bytes that are *already on the device*, and fix two
transfer wastes in the same files. After this stage, `MscHierarchy` mode's GPU path does:
fused gradient kernel → GPU set-dim → ONE flow-terminal call (both directions) → one
gradient D2H — instead of today's: kernel → D2H → CPU set-dim → two label calls that each
re-upload the gradient.

## Verified context

### What the CPU pass does

`TopologicalGradientUsingAlgorithms::setAscendingManifoldDimensions()` —
`include/gi_topological_gradient_using_algorithms.h:430-478`, worker
`rec_man_trace_up_marking` `:395-425`.

It writes the **ldir bits (5-7)** of the 1-byte `GradBitfield`
(`include/gi_discrete_gradient_labeling.h:11-19`; accessors `getDimAscMan/setDimAscMan`
`:114-126`). Three passes, in order:
1. Every cell ← `ldir = 3` (full-lattice read-modify-write) while collecting critical
   cells per dimension under **`omp critical`** (the serialization hotspot).
2. For every critical **edge** (2D saddle): BFS upward through cofacets following pairs,
   stamping `ldir = 2` on each cofacet and its pair — marks the saddle's ascending
   V-path cells (the ridge tree).
3. Every critical **quad** (2D maximum): self-mark `ldir = 1`.

The tags are the **shifted 3D values** (the function is hardcoded 3D and 2D just
inherits): vertices stay 3, saddle-ascending-path cells 2, maxima 1. This is NOT
dimension-truthful in 2D — **replicate it exactly, do not "fix" it**, because:
- every consumer uses the relative predicate `getDimAscMan(c) + 1` (the `rec_tdcr`
  tracer family, `gi_morse_smale_complex_basic.h:538/564/590`; the partitioned tracer;
  `fillGeometry`), and
- the `.grad` **file format carries ldir**: `outputToFile`
  (`gi_discrete_gradient_labeling.h:200-226`) writes the whole byte, and `extractmsc`,
  `extractgraph`, `extractsimplifiedmorse`, `extractfilteredlinesoup` load `.grad` and
  build the MSC **without recomputing** the pass. A GPU pass must write the identical
  bits back into the host array before any `WriteGrad`.

Consumers are ONLY the MSC arc tracers. The forest path never reads ldir (verified: zero
`getDimAscMan` hits in `gi_extrema_region_builder.h`, `gi_extremum_merge_forest.h`,
`gi_extremum_arc_extract_2d.h`, `labeling2d.cu`), which is why `598edd2` could skip
it in forest mode (`SetNeededAscManDims(false)`, `gi_discrete_gradient_computer.h:107-108`,
deferred re-entry `EnsureAscendingManifoldDimensions` `:554-560`).

**Known failure mode (documented twice, keep it impossible):** building the MSC with
ldir never set produces the correct nodes and **essentially no arcs — silently** (the
tracer predicate becomes `0 == 1`). See `gi_discrete_gradient_computer.h:534-541` and
`libs/msc_2d_lib/msc_2d_lib.cxx:463-469`; the gate `forest_lazy_msc` measured arcs
288438 → 48244 when the ensure call was removed.

### Current GPU residency and the transfer wastes

- The fused kernel writes bits 0-4 and leaves **ldir = 0**
  (`gpu_dgrad/cls/cls2d_kernel.cu:184`: `grad[cid] = 1u | (code << 2)`), then the launcher
  frees `d_grad` (`cls2d_kernel.cu:513-514`) and the CPU pass runs on the downloaded
  bytes (`gi_discrete_gradient_computer.h:760-765`, task "Set Dim of Asc Manifolds").
- Waste 1: `MscHierarchy` mode calls `Label2DFlowTerminals` **twice** — asc at
  `msc_2d_lib.cxx:541`, dsc at `:576` — and each call re-uploads the full gradient
  (`labeling2d.cu:134`). Forest mode already does it right: one call, both outputs
  (`:697-699`).
- Waste 2: `run_doubling` (`labeling2d.cu:79-91`) does a **blocking 4-byte D2H of
  the changed flag every pass** (`:87-88`) — a full sync + PCIe round trip × ~15-20
  passes.
- The only persistent device object today is `LabelCtx2D`
  (`gpu_dgrad/labeling2d/label_ctx2d.cu`); it is the pattern to copy.

### Why pointer doubling is the wrong tool here

Pointer doubling answers "what is my root" (descendant → ancestor info flow). Set-dim
needs the opposite aggregation: a cell must learn whether any **descendant** in the flow
forest is a seed (a saddle cofacet). Rather than build reverse edges, exploit two facts:
- **Seed detection is purely local**: a quad is on a ridge iff reachable from a critical
  edge's cofacet; the seeds are just the critical edges — found by a parallel scan, no
  `omp critical` (kills the CPU pass-1 bottleneck outright).
- The marked set is the union of seed→root chains, and chains share tails. One thread
  per seed walking its chain, with early exit on already-stamped cells, does O(total
  ridge length) work — the same total work as the CPU BFS, massively parallel over seeds,
  with **zero extra memory** (stamps go into the resident gradient bytes in place).

---

## Stage B1 — GradCtx2D + set-dim kernels + transfer fixes

### 1. Residency: `GradCtx2D` (new, mirrors `LabelCtx2D`)

`gpu_dgrad/api/dgrad_gpu_api.h` + new `gpu_dgrad/dgrad2d_grad_ctx.cu` (or extend
`cls2d_kernel.cu`):

```cpp
GradCtx2D* CreateGradCtx2D(int64_t X, int64_t Y);            // nullptr => CPU fallback
bool CtxComputeFused(GradCtx2D*, const float* values, Dgrad2DTimings*);
        // fused kernel into ctx-owned d_grad; do NOT free (today :513-514 frees)
bool CtxSetDimAscMan(GradCtx2D*, Dgrad2DTimings*);           // the new pass, below
bool CtxLabelFlowTerminals(GradCtx2D*, int32_t* out_v, int32_t* out_q, Dgrad2DTimings*);
        // reads resident d_grad; BOTH directions, one call
bool CtxDownloadGrad(GradCtx2D*, uint8_t* host);             // the single D2H
void DestroyGradCtx2D(GradCtx2D*);
```

Existing one-shot entry points (`ComputeDiscreteGradient2DFused`,
`Label2DFlowTerminals`) remain as thin wrappers (create → run → download → destroy) so
`cls2d_validate` and other callers are untouched.

### 2. The set-dim kernels (three ordered launches, replicating CPU pass order)

1. **Init**: elementwise `ldir = 3` over all `(2X−1)(2Y−1)` cells. Whole-byte compose is
   safe elementwise (each thread owns its cell): `g = (g & 0x1F) | (3u << 5)`.
2. **Saddle flood**: one thread per critical edge. Generate the exact cell set from the
   CPU `rec_man_trace_up_marking(cid, 2)` semantics — port the CPU function's traversal
   (cofacets of the current cell, skip critical/unassigned, follow `getPair`, stamp both
   the cofacet and its pair), do not re-derive from intuition. Stamping 3 → 2 **clears a
   bit, so plain `atomicOr` is insufficient**: use a 32-bit `atomicCAS` loop composing
   the target byte within its word (`set_ldir(cell, v)` helper). Races are same-value
   idempotent (shared ridge tails all write 2) — the CPU has the same benign race,
   here made explicit and safe. Early-exit: stop the walk on reading a cell already
   stamped 2 (converged tails) — keeps total work O(ridge length).
3. **Maxima**: one thread per critical quad, `ldir = 1`. Launch after pass 2 completes
   (stream order gives the CPU's "maxima stamped last" for free).

Vertices stay 3 (the CPU never touches `criticals[0]`); the 3D `criticals[3]` pass is
empty in 2D. Critical-cell lists: either reuse pass-1 to compact critical edges/quads
into two device arrays (small), or re-scan per pass — either is fine; scan is simpler.

### 3. Pipeline integration (`msc_2d_lib.cxx` + `gi_discrete_gradient_computer.h`)

- The gradient override lambda (installed at `msc_2d_lib.cxx:392-406`) moves to the ctx
  flow: `CtxComputeFused` → (`MscHierarchy` mode only) `CtxSetDimAscMan` →
  `CtxDownloadGrad` into `grad->m_dgrad->LabelArray()`.
- Collapse the two `Label2DFlowTerminals` calls (`:541`, `:576`) into one
  `CtxLabelFlowTerminals` with both outputs (mirror the forest-mode call at `:697-699`).
- **Ensure-flag interaction**: when the GPU stamped ldir, the CPU pass must be skipped
  WITHOUT reintroducing the silent no-arc failure. Add an explicit
  `MarkAscManDimsDone()` on `DiscreteGradientBuilder`, set **only after the gated
  download succeeded** — never inferred from "GPU mode was on".
  `EnsureAscendingManifoldDimensions` (`gi_discrete_gradient_computer.h:554-560`) then
  no-ops on the flag.

### 4. Bundled cheap win in `labeling2d.cu`

Replace the every-pass blocking flag D2H (`:82-89`) with either a lagged check (read the
flag every 4th pass) or a fixed pass count (`ceil(log2(X+Y)) + 1` — the chain length
bound). Keep the flag write for the lagged variant.

Optional preview of C1 (fine to defer): pinned staging for the gradient D2H and
terminal-map D2H in the ctx paths (pageable measured ~11/9 GB/s vs pinned ~21/20).

### Gate

1. **Full-lattice memcmp of ALL 8 bits** of every gradient byte: GPU ctx path vs the CPU
   `setAscendingManifoldDimensions` run on a copy of the same post-kernel bytes. Case
   matrix: constant, ramp, uniform noise, quantized-tie noise, gaussian bumps, odd+even
   dims, tiny (3×3, 5×3) + 1024² + 4096² (the `cls2d_validate` case set is the
   template).
2. **End-to-end MSC parity**: `MscHierarchy` mode with GPU set-dim produces identical
   `criticalPoints()`, `arcGeometry()` (or arc counts if geometry off), and manifold
   labelings to the CPU-set-dim run. **Arc count is the sentinel** for the silent no-arc
   failure.
3. **Timings + transfer audit** recorded in `gpu_dgrad/docs/README.md`: set-dim CPU vs GPU
   (~580 ms → target ~10-20 ms at 3232²); upload count before/after (2 gradient
   re-uploads + 1 download → 0 re-uploads + 1 download); doubling-pass wall time
   before/after the flag fix.

**Files**: `gpu_dgrad/api/dgrad_gpu_api.h`, `gpu_dgrad/cls/cls2d_kernel.cu`,
`gpu_dgrad/labeling2d/labeling2d.cu`, new `gpu_dgrad/dgrad2d_setdim.cu` (or `dgrad2d_grad_ctx.cu`),
`include/gi_discrete_gradient_computer.h`, `libs/msc_2d_lib/msc_2d_lib.cxx`,
`gpu_dgrad/tests/msc2d_gpu_parity.cxx` (gate harness).

---

## Risks

- **ldir semantics drift**: the 2D values are the shifted 3D tags (3/2/1), baked into the
  `.grad` format and every `+1` predicate. Replicate exactly; the all-8-bits memcmp gate
  is the guard. Any temptation to make 2D "dimension-truthful" breaks every downstream
  `.grad` consumer.
- **Silent no-arc MSC**: an ensure-flag set without a real GPU stamp reproduces the
  known nodes-but-no-arcs failure. Explicit `MarkAscManDimsDone()` only after the gated
  download; end-to-end arc-count sentinel in the gate.
- **atomicCAS byte contention**: shared ridge tails and plateaus create same-word CAS
  traffic. Expected negligible (ridge cells are sparse); fallback = a mark-flood variant
  (iterate `mark[succ] |= mark[cur]` with a changed flag, O(chain length) passes) if
  Nsight shows contention. Note: Nsight Compute is not installed on the dev machine —
  use `ptxas -v` + wall timings, as Stages 1-3 did.
- **CPU BFS fidelity**: `rec_man_trace_up_marking` has subtle skip conditions (critical,
  unassigned, pair re-entry). Port from the source, and let the memcmp gate catch any
  divergence — that is what it is for.
- **Ordering between passes**: maxima must land after the saddle flood (a maximum quad
  can sit on a ridge path); default-stream launch order suffices, but if streams are
  introduced later (C1), add an event dependency.
