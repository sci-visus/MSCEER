# gpu_dgrad — GPU discrete gradient path

GPU (CUDA) implementation of the Robins lower-star discrete gradient, developed
on the `cuda-gradient` branch. Stage 0 provided the build scaffolding, the
host-facing API boundary, and a smoke test. Stage 1 adds the 2D
interior-vertex kernel, validated bit-exact against the CPU 2D reference. The
staged plan and its technical appendix live in the branch's master plan
document.

## Layout

| File | Role |
|---|---|
| `dgrad_gpu_api.h` | Host-facing interface. Plain C++ — any MSCEER tool can include it and link `gpu_dgrad` without a CUDA compiler. Documents the `GradBitfield` byte-exact output contract and the generated-table types. |
| `dgrad_gpu_api.cu` | Host-side launcher implementations (`QueryDevice`; `ComputeDiscreteGradient3D` stub until the 3D stage). |
| `dgrad2d_kernel.cu` | Stage 1: 2D Robins lower-star kernel (interior vertices, thread-per-vertex, 9-slot state) + `ComputeDiscreteGradient2D` launcher. |
| `dgrad2d_validate.cxx` | Stage 1 harness: generates the 9-slot topology tables from the real `TopologicalRegularGrid2D` iterators, runs the CPU reference (`MyRobinsNoalloc<...,4,6>::ComputePairing`), and byte-compares every interior-owned cell. |
| `dgrad2d_label_ctx.cu` | Stage 6: persistent `LabelCtx2D` handle - device-resident base labeling plus the remap gather, for interactive relabeling without re-uploading the labeling. |
| `smoke_test.cu` | `gpu_dgrad_smoke` executable: device query, kernel round-trip check, bandwidth measurement. |

Design rule (master plan, Appendix A.5): device code never includes the `gi_*`
template stack. The CPU and GPU implementations share only generated topology
tables and the requirement that outputs `memcmp`-match.

## Building

Requires: CUDA toolkit **>= 12.8** (for native sm_120 / Blackwell; this machine
has v13.0 installed alongside a v12.1 that is too old), Visual Studio 2022 with
the CUDA build customizations installed.

```bash
cmake -S . -B build_cuda -G "Visual Studio 17 2022" -A x64 -T cuda=13.0 -DGPU_DGRAD_ENABLED=ON
cmake --build build_cuda --config Release --target gpu_dgrad_smoke
build_cuda/gpu_dgrad/Release/gpu_dgrad_smoke.exe
```

Notes:
- `-T cuda=13.0` pins the toolkit for the Visual Studio generator (the PATH
  `nvcc` is 12.1 and must not be used). With Ninja/Makefiles instead pass
  `-DCMAKE_CUDA_COMPILER=".../CUDA/v13.0/bin/nvcc.exe"`.
- `CMAKE_CUDA_ARCHITECTURES` defaults to `native` (resolves to `120` on this
  machine); override for other targets.
- `GPU_DGRAD_ENABLED` defaults to `OFF`; a default configure/build is untouched.
- The `LNK4098 LIBCMT/MSVCRT` link warning is the usual benign static-cudart
  vs `/MD` mix; static cudart is kept deliberately so the executables run
  without the CUDA v13 `bin` directory on `PATH`.

## Measured baseline (Stage 0, 2026-08-27)

RTX 5070 Laptop GPU, driver 596.08 (CUDA driver 13.2), toolkit 13.0.88,
MSVC 19.44:

| Quantity | Measured |
|---|---|
| Compute capability / SMs / memory | 12.0 / 36 / 8.55 GB |
| Pinned H2D | **20.9 GB/s** |
| Pinned D2H | **20.1 GB/s** |
| Device DRAM effective (D2D, read+write) | **338 GB/s** |
| PCIe link | Gen5 **x8** active (x16 physical max) |

Calibration vs the master-plan performance model: device DRAM lands inside the
modeled 300–350 GB/s window (kernel-time estimates stand). PCIe is ~21 GB/s,
not the 45–55 GB/s modeled for Gen5 x16 — the link trains at x8 on this laptop.
Streamed end-to-end estimates roughly double their transfer legs (1024^3:
H2D 4.3 GB ≈ 0.2 s, D2H 8.6 GB ≈ 0.43 s, kernel ≈ 0.07 s → ~0.5–0.7 s wall),
which strengthens the case for keeping downstream consumption on-device and
compacting outputs before download.

## Stage 1 results (2026-08-27): 2D interior kernel, bit-exact

`dgrad2d_validate` — **0 mismatches across all cases**: constant image
(pure index tie-breaking), x-ramp, uniform noise, 5-level quantized noise
(odd dims 251x247), gaussian bumps, plus the two timing sizes below. Every
cell owned by an interior vertex byte-matches the CPU
`MyRobinsNoalloc<TopologicalRegularGrid2D,...,4,6>::ComputePairing` output,
and no other cell is touched. (Real-image datasets are not checked into the
repo; the quantized-noise and bumps cases stand in until one is pointed at.)

Measured baselines (CPU = this laptop, 16 OpenMP threads; GPU = RTX 5070
Laptop, kernel is L2-only in Stage 1 — no shared staging yet, labels
CPU-precomputed):

| Case | CPU labeling | CPU pairing | GPU h2d / kernel / d2h |
|---|---|---|---|
| noise 1024x1024 | 30.7 ms | 56.9 ms | 2.4 / **0.61** / 0.5 ms |
| qnoise 4096x4096 | 498 ms | 731 ms | 17.1 / **8.8** / 7.0 ms |

Kernel vs CPU (labeling+pairing) compute: **~140x** at 4096^2; end-to-end
including transfers: **~37x**. The Stage 1 CPU pairing baseline includes the
known dead per-lower-star heap allocation in
`MyRobinsNoalloc::HomotopyExpand` (`gi_modified_robins.h:1005-1013`, 2D-only
code path) — flagged for a separate CPU-side fix; Stage 2 (label fusion)
removes the label arrays from the GPU path entirely.

Run it:

```bash
build_cuda/gpu_dgrad/Release/dgrad2d_validate.exe
```

## Stage 2 results (2026-08-27): label fusion, still bit-exact

`dgrad2d_kernel.cu` now carries two kernels sharing one expansion core:
the Stage 1 kernel (CPU labels, L2-only reads — kept for A/B via the
`max_labels`/`min_labels` parameters) and the fused kernel
(`ComputeDiscreteGradient2DFused`): a 16x16-vertex block stages an 18x18
tile of function values in 1.4 KB of shared memory and computes lower-star
membership and lowest-vertex queries on the fly under the (value, index)
total order. Verified against `gi_max_vertex_labeling.h`: the CPU max
labeling breaks value ties toward the higher index and the min labeling
toward the lower index — exactly argmax/argmin under `IsGreater` — so
fusion reproduces the label semantics without the arrays.

**Both kernels: 0 mismatches on all 7 cases.**

| 4096x4096 | CPU (16 thr) | GPU labels | GPU fused |
|---|---|---|---|
| labeling prerequisite | 474 ms | 474 ms (CPU) | **none** |
| pairing / kernel | 832 ms | 9.3 ms | 10.2 ms |
| h2d / d2h | — | 17.1 / 6.3 ms | 6.1 / 7.3 ms |
| end-to-end incl. prerequisites | 1306 ms | ~507 ms | **23.5 ms (~56x)** |

Fusion's win in 2D is structural, not kernel-time: the fused kernel is
~10% slower than the labels kernel (extra on-the-fly comparisons) but
deletes the 474 ms CPU labeling pass, 2 B/cell of label storage, and
134 MB of label upload. Kernel-only vs CPU labeling+pairing: **~129x**.

Traffic model: ~9.5 B/vertex fused -> measured ~15 GB/s effective — the 2D
kernel is **ALU/latency-bound, not DRAM-bound** (the plan's >60%-of-peak
DRAM gate is a 3D-stage criterion; at 2D sizes the expansion loop dominates).
ptxas (Nsight Compute is not installed on this machine): fused kernel
**50 registers, 0 spill bytes**, 112 B local stack, 1368 B smem/block;
labels kernel 58 registers, 0 spills. No occupancy concern.

## Stage 3 results (2026-08-28): full domain + batched slices, bit-exact

The fused kernel now covers **every vertex** — no CPU fallback anywhere. The
CPU (`MyRobinsNoalloc::ComputeLowerStar`) partitions each lower star into
independent subsets by `boundaryValue(cell)` and expands each separately
(which is why domain corners are always critical minima of their stratum);
the kernel replicates this by computing a per-slot validity + boundary-class
mask and running the shared expansion core once per non-empty class. Batched
slices run via `blockIdx.z` (`ComputeDiscreteGradient2DBatched`); the kernel
is templated on block shape for tuning.

**Gate (strongest test so far): full-array memcmp vs CPU incl. all edges and
corners — 0 mismatches on 9 cases** (previous 7 + tiny 3x3 and 5x3 where
boundary classes dominate), plus per-dimension critical counts equal and
Euler characteristic == 1 on every case. No CPU interior-vs-boundary
inconsistency surfaced.

| Measurement | Result |
|---|---|
| Batched 64 x 512x512 noise slices | **0 mismatches**; kernel 11.3 ms vs CPU serial loop 1654 ms (**~146x**); scaling 0.89x of linear (launch overhead amortized) |
| Block-shape matrix, 4096^2 kernel | 16x16: 11.13 ms, **32x8: 10.43 ms**, 8x8: 11.26 ms — within ~7%, consistent with ALU-bound; 16x16 kept as default |
| 4096^2 full-domain fused | h2d/kern/d2h 5.7 / 11.0 / 7.0 ms vs CPU 538 + 1346 ms |

Found in passing (pre-existing CPU bug, flagged for separate fix):
`RegularGridMaxMinVertexLabeling2D::ComputeOutput` heap-corrupts when the
image height is smaller than the OpenMP thread count (its guarding assert is
compiled out in Release); the harness clamps threads for tiny images as a
workaround. (Since fixed: the labeling now clamps its thread count in 2D and
3D, and the harness workaround was removed.)

## Stage 4 results (2026-08-28): msc_2d_lib integration + flow-terminal labeling

**Integration.** `Accurate2D::DiscreteGradientBuilder` gained a
`SetGradientOverride` hook (called after `ClearAllGradient`; returning false
falls back to CPU pairing), keeping the `gi_*` stack free of CUDA
dependencies. `Msc2D::ComputeOptions::useGpuGradient` wires the fused GPU
kernel through it when the library is built with `GPU_DGRAD_ENABLED`
(`MSC2D_HAS_GPU_DGRAD`); slot tables come from the shared
`dgrad2d_tables.h` builder run on the library's own mesh. A non-CUDA build
compiles unchanged and the option degrades to the CPU path with a notice.

**Parity gate (`msc2d_gpu_parity`): ALL PASS** — the full pipeline with CPU
vs GPU gradient produces identical critical points, arc geometry, and
ascending/descending manifold label images, in BOTH builder modes
(serial and partitioned, 4 partitions), at two persistence values, on noise
and bumps images. Per-partition node/arc/delayed-arc counts are identical.
One test-design note: arc comparison is order-independent with arc ids
excluded, because the CPU serial geometry path creates arcs under
`#pragma omp critical` in a parallel loop — arc ordering is nondeterministic
run-to-run on the CPU itself (`gi_morse_smale_complex_basic.h:512`).

**Flow-terminal labeling (`Label2DFlowTerminals`).** Implements the
offset-doubling observation: like-dimensional flow continuation never needs
the intermediate edge — descending vertex flow steps v -> v + 2*pairdir and
ascending quad flow steps q -> q + 2*pairdir, one gradient byte read and one
add per step on the vertex/quad lattices. Terminal (base-manifold) labels
are computed by in-place pointer doubling. Validated cell-for-cell against
the serial `fillGeometry` base-manifold paint (the reverse traversal of the
same successor):

| Case | CPU fillGeometry paint | GPU flow labeling (h2d/kern/d2h) | mismatches |
|---|---|---|---|
| qnoise 512^2 | 56 ms | 0.1 / 1.1 / 0.3 ms | 0 |
| noise 1024^2 | 363 ms | 0.5 / 1.0 / 1.8 ms | 0 |
| qnoise 4096^2 | 5461 ms | 7.6 / **7.9** / 16.7 ms (**~690x kernel, ~170x end-to-end**) | 0 |

Both labelings (per-pixel basin of minimum, per-quad descending region of
maximum) come back fully labeled — no unlabeled cells — matching the
deterministic single-successor flow partition.

## Stages 6-8 results (2026-08-31): persistent label context, region API, parallel remap

Three changes aimed at the *interactive* cost - what happens when the user
moves the persistence slider, with the gradient and the base labeling already
computed.

**Stage 6 - `LabelCtx2D` (`dgrad2d_label_ctx.cu`).** The first persistent
handle in this library. Everything before it is
alloc -> upload -> compute -> download -> free, which would re-upload the
per-pixel base labeling (the largest buffer in the problem) on every repaint.
A context keeps that labeling resident per direction; a relabel uploads `m`
int32 and runs one coalesced gather,
`out[i] = base[i] < 0 ? -1 : remap[base[i]]`, with an out-of-range compact id
also yielding -1 so the contract is total. `CtxRelabelDevice` skips the
download for callers that consume the labels on the GPU.

Gate: random compact images x random remaps over six shapes, including `m = 0`,
an all-unlabeled image, an all -1 remap, out-of-range ids, 1x1, a re-upload
onto a live context, and a relabel-before-upload that must fail rather than
scribble - all `memcmp`-equal to a CPU reference. Plus a 200x
create/upload/relabel/destroy loop with `cudaMemGetInfo` before and after:
**delta = 0 bytes**.

**Stage 7 - facade integration and region-scale API.** The manifold calls keep
their `LabelImage`-by-value contract; only the final per-pixel projection
changes, becoming one `CtxRelabel` when a context is available and an OpenMP
gather otherwise. The context is populated from the *host* compact base
labeling, deliberately independent of whether that labeling was GPU- or
CPU-painted. New public API (`baseRegionCount` / `baseLabeling` /
`baseToLiving` / `paintLabels` / `deviceBaseLabels` / `paintLabelsDevice`)
lets a caller work at region granularity instead of per pixel.

**Stage 8 - parallel walk-down remap build.** The per-select remap was
allocation-bound, not traversal-bound: one `GatherNodes` into a `std::set` per
living extremum, then an `unordered_map` keyed by raw base identity, then a
projection pass - roughly one set-node allocation plus two hash operations per
base region, per threshold, per direction. The walk-down performs the same
traversal `recGatherNodes` does, but from all living extrema concurrently,
stamping straight into the compact-indexed array through a node -> compact CSR
built once with the base labeling. No `std::set`, no `unordered_map`, no
projection pass, no allocation in the inner loop.

Correctness rests on each base region having exactly one living owner, which
holds for the ascending manifolds of minima and the descending manifolds of
maxima but *not* in general - a cancellation can hand the same manifold to
several surviving neighbours as `merged[1]`, which is what makes the
1-manifolds (2D and 3D) and the 2-manifolds in 3D multi-membership. Those are
out of scope. Rather than trust the argument, the walk counts its own
double-writes and falls back to the reference build if it ever sees one, and
the gate proves map equality directly.

Remap build, 3200x3200 (10.24 Mpix), 2049302 ascending / 1138183 descending
base regions, averaged over four persistence thresholds:

| Builder | Direction | GatherNodes | Parallel walk-down | Speedup |
|---|---|---|---|---|
| serial | ascending | 1602.0 ms | 74.9 ms | **21.4x** |
| serial | descending | 956.3 ms | 75.9 ms | **12.6x** |
| partitioned | ascending | 2290.5 ms | 62.5 ms | **36.6x** |
| partitioned | descending | 1235.1 ms | 62.2 ms | **19.9x** |

The partitioned speedup is larger because the CSR also hoists the
per-constituent lineage vector copy (`GetAscendingLineageCells` returns by
value) out of the per-threshold work. A full repaint at 10 Mpix - remap plus
GPU paint - is now ~75 ms + ~14 ms per direction, down from ~1.6 s + ~19 ms.

`MSC2D_PARALLEL_REMAP=0` restores the original build; that switch is what the
equivalence gate compares against. The gate is CPU-only by construction, so it
runs (and can fail) on a machine with no GPU: two identically computed `Msc2D`
objects, one forced onto each path, compared for exact elementwise equality of
`baseToLiving()` - stronger than comparing label images, since it also covers
compact ids that never label a pixel - across 21 persistence values x both
directions x serial and partitioned x two images.

## Stage A0 results (2026-09-04): tile spine + pinned/streams ring, bit-exact

First stage of the accurate-gradient plan (see the plan file referenced in the
project memory): the fused 2D kernel now runs on **tile views**
(`dgrad2d_internal.cuh: DeviceTileView2D`) - a box of vertices to compute, the
staged values covering their 1-ring, and the owned cell box that receives the
bytes, selected per block by `blockIdx.z`. The in-core batched path is the view
"computed = owned = whole slice", so the tiled and in-core paths share every
line of kernel code; "tiled == in-core" is a property of the code, checked by
the gate rather than hoped for.

Ownership (`tile2d.h`): tile (tx,ty) owns vertices `[ox,ox+tw) x [oy,oy+th)`;
a doubled-lattice cell is owned by the tile owning vertex `floor(c/2)` - the
1-high bridge convention of `steepest_lstar_tiled/localized_morse_complex.h`,
so a tile's cell box is `2tw` wide (`2tw-1` on the last tile) and every cell is
owned exactly once. The kernel computes owned vertices plus a 1-ring (bridge
cells lie in the neighbour's lower stars), stages one more ring of values, and
masks the scatter to owned cells. Halos are re-read from the source; nothing is
exchanged between tiles.

`ComputeDiscreteGradient2DTiled(source, grid, tables, sink, options, stats)`
streams tiles through a `TransferRing` (`gpu_transfer.h`: N slots of pinned-in
/ device-in / device-out / pinned-out buffers, three streams, events) with
`tilesPerLaunch` tiles packed per launch. The image and the gradient are only
ever touched through callbacks (`TileSource2D` / `ByteSink2D`), so neither has
to exist as a whole in host or device memory; the footprint is
`slots x tilesPerLaunch x ((T+4)^2 floats + (2T)^2 bytes)`.

**Gate (`dgrad2d_tiled_validate`): ALL PASS.**

| Gate | Result |
|---|---|
| tiled vs in-core `ComputeDiscreteGradient2DFused`, T in {64,256,512,1000} x (slots,tpl) in {(1,1),(3,1),(3,4)}, 10 cases incl. 3x3, 5x3, 1001x777, 4096^2 | **0 mismatches** everywhere |
| tiled vs CPU `MyRobinsNoalloc::ComputePairing` (same cases) | **0 mismatches** |
| in-core fused vs CPU after the view refactor (+ `dgrad2d_validate` re-run) | **0 mismatches** |
| procedural 8192^2 and **32768^2** (1.07 Gpix; 4.3 GB values + 4.3 GB gradient, never materialized) at T=512, 3 slots, 4 tiles/launch | every byte assigned, Euler chi = 1 on both; device allocated **25.4 MB**, device free-min **identical** (0 B difference) between the two sizes; free before/after delta **0 B**; host pinned 25.4 MB, working set flat |

Measured (RTX 5070 Laptop, driver 596.08, toolkit 13.0):

| Run | H2D / kernel / D2H (summed per-batch windows) | wall |
|---|---|---|
| 32768^2 procedural, T=512, 3 slots, 4 tiles/launch (1024 launches) | 244 / 767 / 241 ms -> **17.9 / 17.8 GB/s** pinned, async | **2.04 s** (4.36 GB up, 4.29 GB down) |
| 8192^2 procedural, same options | 16 / 47 / 15 ms -> 16.6 / 18.4 GB/s | 129 ms |
| 4096^2 qnoise in-core fused (pageable, sync) | 6.3 / 10.4 / 7.4 ms | ~24 ms |
| 4096^2 qnoise tiled T=512, 3 slots, 4 tiles/launch (16 launches) | 19 / 12 / 5.6 ms | 33-37 ms |
| 4096^2 qnoise tiled T=64, 3 slots, 1 tile/launch (4096 launches) | 321 / 280 / 233 ms | 545 ms |

Reading the numbers: the pinned+async ring reaches ~18 GB/s each way on the
streaming run (the pageable in-core path measured ~11/9 GB/s), i.e. the
transfer legs are within ~15% of the Stage 0 bandwidth ceiling once the batch
is large enough. At in-core sizes the tiled path is slower in wall time than
the single-shot path (33 vs 24 ms at 4096^2) because the callback I/O copies
the image into pinned staging and the gradient out of it on the host thread;
that cost is the price of never holding the whole image, and it disappears
when the source/sink are the real out-of-core endpoints. Tiny tiles are launch-
bound (T=64: 4096 launches, 545 ms) - keep T >= 256 and pack tiles per launch;
T=512 with 4 tiles per launch is the default choice. Per-batch event windows
include launch latency, so the GB/s printed for small batches understate the
link; the 32K^2 run is the representative measurement.

Run it (about 25 s; `--quick` skips the timing and footprint cases):

```bash
build_cuda/gpu_dgrad/Release/dgrad2d_tiled_validate.exe
```

## Stage A1 results (2026-09-04): tiled successor resolver, bit-exact

The engine every later stage of the accurate-gradient plan runs on: a
successor function on a lattice (here the discrete flow), resolved per tile
and stitched across tiles, with nothing proportional to the image ever held.
First client: `Label2DFlowTerminalsTiled` (`successor2d.cu`, `frontier2d.h`).

Per tile (lattice box `w x h`) one int32 word per lattice point: tag `00` =
local successor index, `10` = `TERMINAL(v)`, `11` = `FRONTIER(r)` where `r` is
the **ring position** (perimeter index) of the cell whose successor leaves the
box; its exit target (a global lattice index, at most one step outside) is
recorded per ring slot. Ring position as the frontier slot id is
deterministic, needs no hashing, and bounds records per tile by `2(w+h)`.
In-tile pointer doubling `W[i] = W[W[i]]` runs a **fixed `ceil(log2(n))+1`
passes** (a chain in a box of `n` cells is at most `n` long), so the tiled
path has no per-pass host flag at all - the in-core path's blocking flag
read per pass is gone here, not grouped. Terminal cells get a scan-order rank
(CUB exclusive scan) so compact ids are `tileBase + rank`, deterministic in
tile order.

Two streaming phases over the gradient tiles: phase 1 downloads only each
tile's ring words, exit targets and terminal count (perimeter-sized); the
host resolves the frontier graph (`ResolveFrontier2D`: memoized walk over
slots, rooted at terminals - the exit target of a slot is a ring cell of the
neighbour, whose root word is known); phase 2 recomputes the in-tile doubling
(cheap) and a device gather kernel substitutes the resolved payloads straight
into the output tile, which the host copies to the sink. The gradient source
is read twice; that is the price of holding no intermediate words.

**Gate (`succ2d_validate`): ALL PASS.**

| Gate | Result |
|---|---|
| lattice-index mode vs in-core `Label2DFlowTerminals`, both directions, T in {64,256,512} x (slots,tpl) in {(1,1),(3,4)}, 11 cases incl. 3x3, 5x3, 1001x777, 4096^2, and a 2048x512 staircase whose chains cross every tile | **0 mismatches** everywhere |
| compact-id mode vs `ExtNet2D::BuildBaseRegions2D` (consistent bijection + equal region count) | **0 failures** |
| frontier slots per tile <= 4T+4 | max 512 at T=512 (staircase), 81 at T=64 |
| run twice, memcmp | **0 differences** |

Measured (RTX 5070 Laptop), per direction, 3 slots, 4 tiles/launch:

| Case | T | phase 1 | stitch | phase 2 | wall | frontier slots | in-core (both dirs, h2d/kern/d2h) |
|---|---|---|---|---|---|---|---|
| noise 1024^2 | 512 | 6.7 ms | 0.06 ms | 1.7 ms | **8.5 ms** | 798 | 0.5 / 1.1 / 1.2 ms |
| qnoise 4096^2 | 512 | 16.6-22.9 ms | 1.2 ms | 14.7-17.8 ms | **33-42 ms** | 23-26k | 5.2 / 8.1 / 13.3 ms |
| qnoise 4096^2 | 256 | 32-35 ms | 2.4-4.0 ms | 34-43 ms | 70-80 ms | 210-234k | |
| qnoise 4096^2 | 64 | 352 ms | 11 ms | 347-397 ms | 710-765 ms | 4096 tiles, launch-bound | |

Host stitch memory is perimeter-sized: 2.6 MB at 4096^2 / T=512 (20.6 MB at
T=64). Device footprint 63 MB at T=512 (words, flags, ranks, scan scratch for
3 slots x 4 tiles), independent of the image. The tiled path costs ~1.5-2x
the in-core kernel+transfer time at in-core sizes (two passes over the
gradient, 20 fixed doubling passes instead of ~8 adaptive ones); at
out-of-core sizes there is no in-core path to compare with.

Run it (about 40 s; `--quick` skips the 1024^2/4096^2 timing cases):

```bash
build_cuda/gpu_dgrad/Release/succ2d_validate.exe
```

## Stage A2 results (2026-09-05): numeric region decide + hop, on the tiled resolver

`LabelNumericRegions2DTiled` (`successor2d.cu`, `numeric2d.cuh`) labels every
vertex with the living region (basin of a minimum at a persistence) it flows
into NUMERICALLY, using the A1 engine with a different successor function:

- **Decide** first, per vertex: a base minimum -> `remap[base]`; a vertex with
  `f(x) < s*(living)` (lowest saddle joining the living region to any other,
  from the same arcs the forest was built from) is CERTAIN -> its living id.
  A static per-vertex test: no traversal order, no neighbour state.
- **Hop** otherwise: `AdvectThroughVoxelNoCheck` ported line for line (central-
  difference vertex gradients with the zero-on-edge stencil, bilinear
  interpolation with a fixed voxel base, 0.5-cell step cap, squared error
  1e-4, x1.8/x0.5, descent), computed on the fly from the staged values. It
  ends on entering a voxel whose 4 corners are certain with one living id
  (-> that id), on leaving the tile's owned box, on touching the domain
  boundary, at `hopCells` voxels, or at `iterCap` iterations; then it SNAPS to
  the lowest corner of the voxel it was in, **accepted only if strictly
  below the start vertex under the (value, index) order**, else the discrete
  V-path step. That guard is what makes the field acyclic for any vector
  field - the gate checks it edge by edge rather than trusting the argument.
- Frontier slots are now indexed by the exit **target's** position on the
  tile's outer ring (a numeric hop can leave from an interior vertex, but its
  snap target is always within one cell of the box); the A1 flow kind uses
  the same scheme, targets are implied by geometry, and nothing but inner-
  ring words crosses the bus in phase 1.

The tile size is an ALGORITHM parameter (hops clamp at the tile box), so
runs at different T legitimately differ by a few 0.01-0.2% of vertices; what
must be identical is the same T at different residencies.

**Gate (`numeric2d_validate`): ALL PASS.**

| Gate | Result |
|---|---|
| every successor chain ends at a terminal; no cycles; every vertex labeled | 0 failures, 6 cases up to 2048^2 (1.9M chain edges walked) |
| strict (value, index) decrease along EVERY successor edge | 0 violations |
| same T, residency (1 slot,1 tile/launch) vs (3,4) memcmp; run twice | 0 differences |
| numeric labels vs host-integrated basins (`AdaptiveEulerAdvector2D` from every vertex to its minimum, double precision, no snapping) | bumps 512x384: **99.991%** (discrete forest labels: 98.92%); multiscale 512^2: **99.377%** (discrete: 99.19%); noise: 86-87% (discrete: 79%) - inside the certain core agreement is 99.997-100% |

Termination histogram (multiscale 2048^2, hopCells 32, iterCap 500): certain
24.8%, minimum 6.9%, corner-certain exit 29.4%, cell cap 20.7%, iteration
cap 17.5%, boundary 0.4%, guard fallback 0.02%. Only ~38% of vertices ever
integrate; the certain test and the corner-certain exit decide the rest.

Measured (RTX 5070 Laptop), in-core (one tile) vs tiled T=512, 3 slots, 4
tiles/launch:

| Case | in-core phase1 + phase2 | tiled T=512 wall | cross-T diff vs in-core |
|---|---|---|---|
| bumps 512x384 | 1.4 + 1.1 ms | 3.6 ms | T=64: 0.057%, T=256: 0.034% |
| multiscale 512^2 | 9.2 + 9.3 ms | 21 ms | T=64: 0.22%, T=256: 0.03% |
| multiscale 2048^2 | 138 + 142 ms | 276 ms | T=64: 0.21% (993 ms, launch-bound), T=256: 0.035%, T=512: 0.017% |

Phase 2 recomputes the hops (the price of holding no intermediate words);
the hop kernel is ~35 ms per Mvertex on this field, dominated by the
iteration-capped chains in noisy areas. The boundary-rework pass (P5b of the
plan) is not implemented yet: the agreement numbers above are without it.

Run it (about 40 s; `--quick` skips the 2048^2 case):

```bash
build_cuda/gpu_dgrad/Release/numeric2d_validate.exe
```

## Stage A3 results (2026-09-05): consistency pass

`ConsistencyPass2DTiled` (`successor2d.cu`, kind `kKindConsistency`) is the
GPU replacement for the CPU's `IsolatedCCRegionRemoverNEW`: it makes a
vertex region labeling *discretely consistent*, meaning every non-minimum
vertex has a lower AXIS neighbour carrying its own label - the property the
constrained lower-star pairing needs so it never manufactures a boundary
minimum. Successor of a vertex = its lowest axis neighbour (under the
(value, index) order) with the same label, else its lowest axis neighbour
overall; a vertex with no lower axis neighbour is a Robins minimum and a
terminal carrying its own label. The output label is the chain terminal's,
so the property holds **by construction** (a vertex and its successor share
a terminal). Successors are +-1 steps, so the result is independent of the
tile size, and the pass is idempotent. It runs on the same engine as A1/A2
(inputs: values + region words with a 1-vertex halo).

**Gate (`numeric2d_validate`, A3 section): ALL PASS.**

| Case | violations before -> after | relabeled | idempotent | tiled T=64 == in-core | ms |
|---|---|---|---|---|---|
| bumps 512x384 | 3 -> **0** | 5 (0.003%) | yes | yes | 2.2 |
| multiscale 512^2 | 494 -> **0** | 1855 (0.71%) | yes | yes | 1.7 |
| quantized noise 251x247 | 4009 -> **0** | 5728 (9.2%) | yes | yes | 1.2 |
| noise 256^2 | 4727 -> **0** | 6227 (9.5%) | yes | yes | 1.9 |

Agreement with the host-integrated basins after the pass: bumps 99.989%
(from 99.991%), multiscale **99.418%** (from 99.377%), quantized noise 87.2%
(from 85.8%), noise 89.2% (from 86.7%) - i.e. enforcing discrete consistency
does not cost numeric accuracy on structured fields and helps on noise.

The quad-lattice (maximum / descending-region) direction of A2/A3 is the
symmetric follow-up and is not implemented yet.

## Stage A4 results (2026-09-05): numeric separatrix lines, deterministic and thin

`ComputeLineClasses2DTiled` (`lines2d.cu`, `lines2d.cuh`) produces the per-cell
line classes the CPU's accurate 2D path paints (1 = ascending quad-edge-quad
chain toward a maximum, 2 = descending vertex-edge-vertex chain toward a
minimum), the input of the constrained lower-star pairing (Stage A5).

- **Seeds** replicate the CPU's depth-10 BFS from every critical edge along
  both V-paths (`gi_discrete_gradient_computer.h:237-271, :382-416`), on the
  device from the pair codes, clamped to the tile's staged gradient box.
- **Walks** are chains of numeric hops (the A2 integrator) of
  `hopCellsLine` = 4 voxels, digitized with the CPU's plane-crossing rules,
  snapping to the nearest vertex (descending, guarded by the (value, index)
  order) or the containing quad (ascending, guarded by the quad order), with
  the CPU's **500-iteration budget per streamline** (`walkIter`) - without
  it a walk in a noisy field costs thousands of times the CPU's work, since
  every hop would burn its own iteration cap before snapping one cell.
- **Thinness without races.** The CPU keeps lines thin by stopping a line
  the instant it touches a painted cell - earlier in *thread timing*. The
  GPU reproduces that rule deterministically in two phases: phase 1 walks
  every seed and records at each lattice point the most senior walk passing
  it (`atomicMin` of a 64-bit key: saddle value in the CPU's processing
  order, then seed index, then a saddle hash); phase 2 walks again and marks
  only until the first point a more senior walk owns, marking that point so
  the junior line joins the senior one. The union of marks is a function of
  the gradient alone.
- **Cross-tile**: a walk leaving its owned box records its key on the tile's
  outer ring for the exit target (again an `atomicMin`: the most senior walk
  resumes); the host re-seeds the owning neighbour tile in later rounds.
  Owner keys (8 B per lattice point, two lattices) are the one array
  proportional to the image in this stage: a tile-streamed state, kept
  in-core here.

**Gate (`lines2d_validate`): ALL PASS.**

| Gate | Result |
|---|---|
| run twice; same T at residency (1,1) vs (3,4) and the (1,4)/(3,1) factor test | **0 differences** on all 5 cases (up to 2048^2, 525k saddles) |
| every critical edge carries both classes | 0 missing |
| thickness (line vertices / quads with all 4 incident edges marked), GPU vs the same metric on a SERIAL CPU run of the line integrators | bumps 0.00% vs 0.00%; multiscale **0.53% vs 0.67%**; noise 1.94% vs 2.20% (asc) - the GPU lines are as thin as or thinner than the CPU's |
| line cells vs CPU (serial), Jaccard | bumps asc 2916 vs 2907 cells (J 0.38), multiscale 373k vs 419k (J 0.75), noise 111k vs 110k (J 0.87) - equal budgets, different but equally valid line placements where lines merge |

Measured (RTX 5070 Laptop; ascending + descending, both phases):

| Case | in-core | tiled T=64 (3 slots, 4 tiles/launch) | serial CPU line integrators |
|---|---|---|---|
| bumps 512x384 | 21 ms | 30 ms (10 rounds, 67 resumes) | 8 ms |
| multiscale 512^2 (33k saddles) | 37 ms | 123 ms (10 rounds, 14.7k resumes) | 474 ms |
| noise 256^2 (20k saddles) | 12 ms | 34 ms | 552 ms |
| multiscale 2048^2 (525k saddles) | **614 ms** | 1.78 s (12 rounds, 267k resumes) | (not run) |

Walk-end histogram at 2048^2: merged into a senior line 7.7M, terminal
494k, boundary 18k, budget 29k - i.e. almost every seed's walk is a short
stub joining an existing line, exactly the CPU's structure. Cross-T
differences (T=64 vs in-core: 0.25% of cells at 2048^2) are the clamped-hop
effect documented in A2. The residency mismatch seen during development
(2 bridge-edge cells with 4 tiles per launch) did not survive the rewrite of
the seed staging; the factor test stays in the gate.

Run it (about 60 s; `--quick` skips the 2048^2 case):

```bash
build_cuda/gpu_dgrad/Release/lines2d_validate.exe
```

## Stage A5 results (2026-09-05): constrained fused kernel, bit-exact vs the CPU

`ComputeDiscreteGradient2DConstrained` / `...ConstrainedTiled`: the fused
lower-star kernel with a per-cell restriction class byte (the A4 line
classes, or any labeling). Each lower star is split into subsets keyed by
`class * maxDim() + boundaryValue` and every subset is homotopy-expanded
independently in ascending key order - exactly the 4-argument
`MyRobinsNoalloc(mesh, maxv, restriction, grad)` (`gi_modified_robins.h:
1281-1340`), including that key's collisions, so the output is a byte-exact
match of the CPU's conforming pairing rather than a re-interpretation of it.
The kernel stages a `(2B+3)^2` class tile next to the value tile; a null
class box is the Stage 3 boundary-class grouping, so the unconstrained path
is unchanged (re-gated). Every vertex is recomputed - on the device that is
cheaper than finding the touched lower stars, and identical for single-
class stars.

**Gate (`dgrad2d_constrained_validate`): ALL PASS.**

| Gate | Result |
|---|---|
| memcmp vs CPU constrained `MyRobinsNoalloc::ComputePairing` on all-zero, uniform random {0,1,2} per cell, random thin chains, and the REAL A4 line classes; 7 cases incl. 5x3, 251x247, 2048^2 | **0 mismatches** everywhere |
| all-zero classes == unconstrained `ComputeDiscreteGradient2DFused` | **0 mismatches** |
| tiled (T = 64, 256; residency (1,1) and (3,4)) == in-core | **0 mismatches** |
| every cell assigned; Euler characteristic | chi = 1 on every case and labeling |

Kernel time is essentially the unconstrained kernel's (0.11-0.55 ms at
512^2, the CPU's constrained pass 9-18 ms there). Note the critical-cell
counts under the A4 classes at persistence 0: multiscale 512^2 goes from
26k/33k/7k (min/saddle/max) unconstrained to 72k/140k/68k constrained - the
CPU produces the identical numbers (it is a memcmp match). Line constraints
manufacture low-persistence critical cells along the lines by design; the
`lineSaddles = livingOnly` knob of the plan (trace only saddles alive at the
target persistence) is the lever to reduce them, and belongs to A6.

Run it (about 40 s; `--quick` skips the 2048^2 case):

```bash
build_cuda/gpu_dgrad/Release/dgrad2d_constrained_validate.exe
```

## Stage A6 results (2026-09-05): the GPU accurate gradient through msc_2d_lib

`Msc2D::ComputeOptions::useGpuAccurate` (implies `useGpuGradient`) routes the
accurate step through a new builder hook,
`Accurate2D::DiscreteGradientBuilder::SetAccurateOverride` - called in place
of the CPU's numeric re-integration + constrained re-pairing block with the
base pairing already in the gradient; returning false falls back to the CPU
block, so `gi_*` stays CUDA-free exactly as with `SetGradientOverride`. The
override runs A4 (`ComputeLineClasses2DTiled`, class bits from the base
gradient) and A5 (`ComputeDiscreteGradient2DConstrained`) and writes the
conforming pairing back in place; everything downstream (MSC, hierarchy,
forests, label context) is untouched. Knobs: `gpuTile`, `gpuLineHopCells`,
`gpuLineWalkIter`.

**Gate (`msc2d_gpu_parity`, new section): ALL PASS.**

| Gate | noise 512^2 | bumps 512x384 |
|---|---|---|
| `useGpuAccurate` run twice: critical points, arcs, both manifold images at 5 thresholds | **0 differences** | **0 differences** |
| the CPU accurate path run twice (its own baseline) | 2.4M differing label entries | 1.7M differing label entries |
| critical points GPU vs CPU accurate, @0.02 range and at base | 155805 vs 155805; 125035 vs 125035 | 17 vs 17; 13 vs 13 |
| mutual segmentation agreement GPU vs CPU accurate, minimum over the sweep (incl. persistence 0) | asc 88.3%, dsc 97.2% | asc 78.8%, dsc 84.2% |
| accuracy off + `useGpuAccurate` == `useGpuGradient` | 0 differences | 0 differences |

The CPU accurate path is not reproducible - its parallel saddle loop races
on `HIT_PREASSIGNED` - so agreement with it is a statistic, not a memcmp;
its own run-to-run figure is the ceiling. The segmentation agreement
minima are taken at persistence 0, where the two line placements (equal
budgets, different but equally valid choices where lines merge) shuffle
many tiny base regions; at the thresholds a user selects the agreement is
higher, and the critical-point counts are identical.

Measured accurate stage only (lines + constrained kernel, both directions):

| Case | GPU accurate stage | CPU accurate stage (SaddleIntASC + DSC + ConformingGrad) |
|---|---|---|
| noise 512^2 | **41 ms** (lines 41 ms, kernel 0.6 ms) | 611 ms (338 + 246 + 27) |
| bumps 512x384 | **29 ms** | 202 ms (103 + 77 + 22) |

`compute()` wall (serial builder, MSC hierarchy, arc geometry on): noise
954 -> 482 ms, bumps 504 -> 474 ms - the remaining time is the serial MSC
build, which is outside this plan's scope. Small images no longer pay for
worst-case residency: every driver clamps slots and tiles-per-launch to the
tile count, and the line driver sizes its seed lists to an exact per-tile
bound (2(bfsDepth+1)+2 per critical edge) instead of a fixed capacity -
a fixed capacity that could overflow was found to make the result depend
on atomic order.

Not implemented in this pass (knobs named in the plan): `lineSaddles =
livingOnly`, the boundary-rework pass (P5b), the quad-lattice direction of
the region stages, and the region-conforming key (`regionConform`); the 2D
product path does not need them.
