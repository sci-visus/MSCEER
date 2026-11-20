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
| `dgrad_gpu_api.cu` | Host-side helpers (`QueryDevice`). |
| `dgrad2d_kernel.cu` | Stage 1: 2D Robins lower-star kernel (interior vertices, thread-per-vertex, 9-slot state) + `ComputeDiscreteGradient2D` launcher. |
| `dgrad2d_validate.cxx` | Stage 1 harness: generates the 9-slot topology tables from the real `TopologicalRegularGrid2D` iterators, runs the CPU reference (`MyRobinsNoalloc<...,4,6>::ComputePairing`), and byte-compares every interior-owned cell. |
| `dgrad2d_label_ctx.cu` | Stage 6: persistent `LabelCtx2D` handle - device-resident base labeling plus the remap gather, for interactive relabeling without re-uploading the labeling. |
| `dmt3d_gradient.cu` | 3D discrete Morse theory rules: `ComputeLowerStarGradient3DGPU` (PLS, PLS-BR) and `ComputeDMFGradient3DGPU` (Shivashankar-Natarajan), running the device-safe core of `include/dmt_core` unchanged. |
| `dmt3d_validate.cxx` | Gate for the above: byte-compares every output (gradient, step/part journals, DMF stage and order, stats) against the CPU raw drivers of `include/gi_dmt_*.h`, plus a timing table. |
| `tile3d.h` | 3D tile geometry (ownership by vertex floor(c/2); computed, staged and class boxes), callback I/O types and in-core adapters, tiled options and stats. Plain C++. |
| `dgrad3d_tables.h` | Builds the 27-slot tables of the 3D lower-star kernel from `Explicit3x3x3SmallRegularGrid` iterators and checks them against `TopologicalRegularGrid3D`'s. |
| `dgrad3d_kernel.cu` | 3D lower-star (Robins) kernel: fused, full domain, optional restriction classes, bitmask expansion state, int32 tile-local indexing. |
| `dgrad3d_tiled.cu` | Tiled driver (TransferRing, three streams) and the in-core and constrained entry points, which run the same driver with one tile. |
| `dgrad3d_validate.cxx` | Gate for the above: memcmp against `SlidingWindowRobinsNoalloc::ComputePairing_sliding` and constrained `MyRobinsNoalloc`, tiled == in-core, fixed-footprint procedural runs, timing table. |
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

## 3D discrete Morse theory rules (2026-11-19), bit-exact

`ComputeLowerStarGradient3DGPU` and `ComputeDMFGradient3DGPU` take the same
inputs and produce the same outputs as the CPU raw drivers
`GInt::dmt::ComputeLowerStarGradient3D` / `ComputeDMFGradient3D` (flat
uint32 ranks in, one GradBitfield byte per cell out, optional journals). The
kernels call the `include/dmt_core` functions directly; this directory only
adds the launch structure:

- PLS / PLS-BR: one thread per vertex (3D blocks of 32x4x2), each collapsing
  its 27-slot lower star and writing its cells; then a census kernel counts
  critical and unassigned cells (block reduction + one atomic per block).
- DMF: one top-facet kernel per dimension d = 1..dmax, then one kernel per algorithm per
  dimension d = 0..dmax-1 over the cells of dimension d (enumerated by parity
  class), guard and claim counts by block reduction, and a final sweep that
  makes unpaired cells critical. The optional (dimension, key) order is an
  LSD sequence of stable CUB radix sorts per dimension, one per key component.

`dmt3d_validate` runs 39 shapes (single vertex, degenerate axes, slabs,
cubes, permutation and non-permutation ranks) x 5 configurations, a host-side
cell-by-cell check of the key cache (below) on the 38 shapes up to 2^21
cells, and three large cases: 246 runs, 0 mismatches. `dmt_py` exposes the drivers as
`Gradient3D.pls/pls_br/dmf(..., device='gpu')` when built with
`GPU_DGRAD_ENABLED`.

Timings (RTX 5070 Laptop vs 16 OpenMP threads; GPU = end to end with
pageable transfers, after a warm-up call):

| Case | CPU ms | GPU ms | Speedup | h2d | kernel | d2h |
|---|---|---|---|---|---|---|
| PLS 200^3 | 1024 | 63 | 16.3x | 2.3 | 50.6 | 6.2 |
| PLS-BR 200^3 | 1027 | 69 | 14.9x | 2.6 | 54.1 | 8.1 |
| PLS-BR 200^3 + journal | 1017 | 226 | 4.5x | 10.8 | 140.1 | 55.2 |
| DMF 200^3 | 2623 | 355 | 7.4x | 21.5 | 287.1 | 11.1 |
| DMF 200^3 + order | 31942 | 662 | 48.3x | 21.8 | 459.3 | 73.5 |
| PLS 256^3 | 2168 | 141 | 15.4x | 5.4 | 114.6 | 12.0 |
| PLS-BR 256^3 | 2139 | 142 | 15.0x | 6.6 | 113.5 | 13.9 |
| DMF 256^3 | 5790 | 562 | 10.3x | 45.1 | 420.0 | 25.1 |
| PLS 384x256x192 | 2388 | 152 | 15.7x | 5.5 | 122.6 | 15.0 |
| DMF 384x256x192 | 6801 | 592 | 11.5x | 50.7 | 432.6 | 27.3 |

DMF key cache (`include/dmt_core/dmt_dmf_keycache.h`). The G0 of a cell's
recursive key is its largest facet, which is exactly the g0 direction the
top-facet table stores, and G1 is the facet across the same axis. So the
top-facet bytes already encode every key: once dimension d-1 is filled, the
key of a d-cell is a d-level walk of 2^d - 1 byte reads plus 2^d rank reads,
with no comparisons. The top-facet kernels therefore run one dimension at a
time, and the pass and order kernels read keys back the same way, with no
memory beyond the one byte per cell the table already used. Against
recomputing every key from the ranks this cut the DMF kernel time about 1.6x
(447.6 -> 287.1 ms at 200^3, 698.3 -> 420.0 ms at 256^3) and the order
kernels about 1.7x. The lower-star rows are from the earlier run; that code
is unchanged and repeat runs vary by up to 2x on this laptop GPU.

Profile after the change (Nsight Systems, whole validator): DMF pass kernels
about 60% of DMF kernel time, top-facet kernels about 40%. Both still do
several 64-bit cell-id <-> coordinate divisions per cell (`Cubical3D`);
32-bit or coordinate-based indexing is the next thing to try.

## 3D lower-star (Robins) gradient (2026-11-20): tiled, bit-exact

GPU version of the gradient `ondemandaccurate` computes first:
`SlidingWindowRobinsNoalloc::ComputePairing_sliding` (interior vertices
through the sliding 3x3x3 window, boundary vertices through
`MyRobinsNoalloc<...,5,4>::ComputeLowerStar`), which equals
`MyRobinsNoalloc<...,4,6>::ComputePairing`. The constrained variant
reproduces the 4-argument `MyRobinsNoalloc<...,4,6>` used for the accurate
pass.

API (`dgrad_gpu_api.h`, replaces the old `ComputeDiscreteGradient3D` stub):
`ComputeDiscreteGradient3D`, `ComputeDiscreteGradient3DTiled`,
`ComputeDiscreteGradient3DConstrained(Tiled)`. Every entry point returns
false on a CUDA failure so the caller can fall back to the CPU.

Design:
- **Tiles.** Ownership as in 2D: a cell belongs to the tile owning vertex
  floor(c/2). An owned cell's vertices lie in [o, o+t], so the computed
  vertex box is the owned box plus one vertex on the high side only, and the
  staged values add a one-vertex ring. Each cell is in exactly one lower
  star, so there is one pass and the writes are disjoint. Tile sizes can
  differ per axis (cubes, slabs, one tile for in-core). Host bookkeeping is
  int64. The device view carries int32 per-axis coordinates (each axis below
  2^30) and int32 tile-local pitches.
- **Kernel.** One thread per vertex from a 3D launch, so no id-to-coordinate
  division. Values come from a shared-memory block. A 26-compare mask of
  "window vertex below the center" gives lower-star membership for all 27
  slots, and lowest vertices are computed on the fly, so there are no max/min
  label arrays. Only vertices of one 3x3x3 window are ever compared, and the
  window position w = (ex+1)+3(ey+1)+9(ez+1) orders them like the global
  index, so the tie-break needs no int64 index. All expansion state is
  27-bit masks: lower star, subset, paired, the missing-facet counter as two
  bit planes (a lower-star cell has at most dim <= 3 window facets, exactly
  the ones containing the vertex), and the pair code as three bit planes.
  There are no per-slot arrays and no local memory.
- **CPU orders replicated.** Insertion order is ascending slot (the
  AdjacentCellsIterator order of both meshes). Candidate cofacets are taken
  in CofacetsIterator order, and the earliest candidate wins a tie of lowest
  vertices. The lowest unpaired d-cell becomes critical, and the d-cell scan
  restarts after every pair. Each lower star is split by
  `boundaryValue(cell)`; with classes the key is `class*3+bv`, clamped at
  `LMAXDIM*LMAX_LABEL-1` (23), so the CPU's key collisions (class 1, bv 0 ==
  class 0, bv 3) are reproduced. On degenerate axes
  `Compress6NeighborOffsetToByte` maps two axes to the same code (X == 1
  makes a y step look like an x step), and the driver replicates that per
  volume.
- **Tables.** They are generated from the `Explicit3x3x3SmallRegularGrid`
  iterators and required to equal `TopologicalRegularGrid3D`'s. The small
  grid's facet and cofacet iterators visit one entry past each id list
  (`m_end = 1 + count`), either the zero padding or the next row's count.
  That id is below 7 and never in a lower star, so it does not change the
  CPU result; the table builder skips out-of-window ids.
- **Streaming.** This is the 2D TransferRing: pinned in/out buffers per
  slot, with H2D, kernel and D2H on three streams, overlapped with the
  host's source and sink callbacks. The in-core path is the same driver
  with one tile; its adapters copy on a few host threads.

Gate (`dgrad3d_validate`, exit 0 = all pass; `--quick` skips 256^3, the
footprint runs and the timings):

| Gate | Result |
|---|---|
| tables: small-grid == `TopologicalRegularGrid3D` iterator orders | PASS |
| GPU in-core vs CPU `ComputePairing_sliding`, full cell array, 24 volumes: random / quantized (0..4) / binary / constant / x-ramp / diagonal ramp / z-down / terrace / signed zeros / bumps; odd sizes; 1x1x1, 2x2x2, 3x3x3, 1x1x9, 7x1x1, 1x6x5, 4x1x3, 9x8x1, 2x9x2, 64x3x2, 65x66x67, 256^3 noise, 256^3 qnoise, 301x157x203 binary | **0 mismatches** |
| CPU sliding vs `MyRobinsNoalloc::ComputePairing` (every volume up to 300k vertices) | equal |
| tiled vs in-core: T in {1, 2, 3, (5,3,2), 7, 16, (X,Y,1), (X,4,3)} on small volumes; {64, (100,37,61), (128,128,32), (X,Y,16)} on large ones; slots 1-3, 1-7 tiles per launch, all three block shapes | **0 mismatches** |
| constrained in-core vs CPU 4-arg `MyRobinsNoalloc` (classes 0..3 in blocks plus 6% in 4..7, so keys collide and clamp); constrained tiled vs in-core (2 tile shapes); all-zero classes vs unconstrained; 14 volumes incl. 256^3 | **0 mismatches** |
| procedural 512^3 and **1024^3** (1.07 G vertices: 4.3 GB values + 8.6 GB cells, more than the 8 GB device, never on the host), T=(128,128,64), 3 slots, 2 tiles/launch | every cell assigned, Euler chi = 1 on both; tiled 512^3 hash == in-core hash; device allocated **77.9 MB** for both, device free-min **identical** (0 B difference), free before/after delta **0 B** |

Timing, 256^3 uniform noise (RTX 5070 Laptop vs CPU with 16 OpenMP threads,
best of 3 runs; GPU wall times are end to end from a host array to a host
array, pinned staging copies included; kernel/h2d/d2h are per-batch event
windows summed, so overlapped legs add up to more than the wall time):

| Path | ms | vs CPU pairing | vs CPU labeling + pairing | kernel | h2d | d2h |
|---|---:|---:|---:|---:|---:|---:|
| CPU max/min labeling | 212.5 | | | | | |
| CPU `ComputePairing_sliding` | 688.7 | 1.0x | | | | |
| CPU labeling + pairing | 901.2 | | 1.0x | | | |
| GPU in-core (one tile) | 65.9 | 10.5x | 13.7x | 11.86 | 3.01 | 5.56 |
| GPU tiled T=64^3, 3 slots, 4 tiles/launch | 50.4 | 13.7x | 17.9x | 15.21 | 20.98 | 11.34 |
| GPU tiled T=(128,128,32), 3 slots, 2 tiles/launch | 37.4 | 18.4x | 24.1x | 14.54 | 14.60 | 6.78 |
| GPU tiled T=(256,256,16), 3 slots | **36.6** | **18.8x** | **24.6x** | 14.21 | 21.04 | 6.70 |
| GPU tiled T=(256,256,32), 3 slots | 42.1 | 16.4x | 21.4x | 12.83 | 18.96 | 6.16 |

Notes:
- The GPU needs no labeling pass, so "vs labeling + pairing" is the
  like-for-like number for `ondemandaccurate`'s first gradient.
- The kernel itself is 11.5 ms for 16.7 M vertices (one tile, no overlap).
  The 197 MB that must cross PCIe (64 MB values up, 133 MB cells down) take
  about 9.4 ms at the measured 21 GB/s, so slab tiles end up limited by
  transfers and host copies.
- In-core is slower than tiled because one tile cannot overlap anything, and
  each call allocates and pins about 200 MB of staging buffers.
- Block shapes 32x4x2 (default), 16x8x2 and 8x8x4 are within about 10% of
  each other.
- 1024^3 procedural, tiled with a flat 78 MB device footprint: 4.1 s wall.
  Of that, 2.4 s is spent in the source callback (procedural values on the
  host) and 1.6 s in the sink callback (hashing and counting 8.6 GB).

## Tiled DMT rules: multi-pass DMF, PLS, PLS-BR (2026-11-20), bit-exact

`ComputeDMFGradient3DTiled` (with and without the second pass) and
`ComputeLowerStarGradient3DTiled` (PLS, PLS-BR) in `dmt3d_tiled.cu` are the
blocked versions of the in-core DMT drivers. They follow the measured halos
of the DMT blocking notes:

| Rule | Staged vertices per axis (owned vertices [o, o+t)) | Halo beyond the vertices [o, o+t] the owned cells touch |
|---|---|---|
| PLS, PLS-BR, DMF first pass only | [o-1, o+t+2) | 1 below, 1 above |
| DMF with the second pass | [o-2, o+t+2) | 2 below, 1 above |
| Robins (`dgrad3d_*`, above) | [o-1, o+t+2) | 1 below, 1 above |

All of them use the same ownership rule: a cell belongs to the tile owning
vertex floor(c/2) on every axis.

- **One visit per tile.** Each tile treats its staged box as a cubical complex
  of its own, `dmt::Cubical3DT<int>`. All of its passes run on that local data
  as separate kernels on one stream, so the kernel boundary is the barrier
  between passes. For DMF the passes are: top facets for d = 1..dmax, then
  algorithm 1 for d = 0..dmax-1, then algorithm 2. There is no global sweep
  per pass and no exchange between tiles. Halo cells are recomputed and never
  leave the device. Counters (critical cells, guards, claims, orphans,
  unassigned) count owned cells only.
- **32-bit, no divisions.** Host bookkeeping is int64. Kernels use int32
  tile-local coordinates and launch in 3D over one parity class of cells at a
  time, so a thread knows its cell's coordinates with no division. DMF uses
  new coordinate variants of the key-cache functions
  (`DMFKeyFromTopQ`, `DMFTopFacetsFromTopQ`,
  `DMFLowestSelectingCofacetFromTopQ` in `dmt_core/dmt_dmf_keycache.h`),
  templated on the index type. `dmt::Cubical3D` is now
  `Cubical3DT<int64_t>`, so the CPU drivers are unchanged and
  `dmt3d_validate` still passes. (The key walks were rewritten again later
  in this file, see the DMF kernel speedup section.)
- **Not tiled:** the DMF (dimension, key) order output, which is a global
  sort, and the PLS step/part journals.

Gate (`dmt3d_tiled_validate`, exit 0 = all pass):

| Gate | Result |
|---|---|
| tiled vs in-core GPU vs CPU raw driver (CPU on volumes up to 2^22 cells), gradient + DMF stage memcmp + every stats field; DMF, DMF alg 1 only, PLS, PLS-BR on 14 volumes (1x1x1, 2x2x2, 3x3x3, 1x5x7, 6x1x4, 5x6x1, 2x9x3, 9x2x2, 17x13x11, 24x20x18, 33x7x29, 40x41x39, 96x80x72, 200x130x90, sparse and dense ranks); tiles T = 1, 2, 3, (5,3,2), 7, 16, (X,4,3), (X,Y,4), 32, 64, (100,37,61), (128,128,16), so tiles smaller than the halo and tiles on the boundary are included | **440 tiled runs, 0 mismatches** |
| halo sweep: 24 random 22^3 fields, T=6, staged halo (lo, hi) for lo, hi in 0..2 (each tile a truncated complex, so this measures the halo of the blocking itself) | DMF: every halo below **2/1** fails on 24/24 fields (1/1: 258 cells); 2/1 and 2/2 give 0/24. DMF alg 1, PLS, PLS-BR: everything below **1/1** fails 24/24; 1/1 and above give 0/24. This matches the probe in the notes. |
| DMF on procedural 512^3 and **1024^3** ranks (bijective hash, so ranks are pairwise distinct; neither the volume nor the gradient exists on the host), T=(128,128,64), 3 slots | every cell written, Euler chi = 1 on both; tiled 512^3 hash == in-core hash; device allocated **58.1 MB** for both; free-min identical (0 B); leak 0 B |

## GPU-only benchmark at 1024^3 (2026-11-20)

`gpu3d_bench [N]` runs every tiled 3D GPU path on an N^3 volume. The values
(float) and ranks (uint32, a bijective hash, so pairwise distinct) are held
in host memory, and so is an 8.6 GB host gradient buffer that is reused.
Copies between host and pinned memory run on OpenMP threads, 3 slots. At
1024^3 the volume has 1.07 G vertices and 8.58 G cells, more than the 8 GB
device, so no in-core GPU driver can run it. No CPU gradient is run. Each run
is checked: the driver succeeds, the Euler characteristic of the critical
cells is 1, and device free memory is unchanged afterwards (after a
warm-up; the first launch of each kernel loads its module).

| Method | Tile | Wall ms | G vertices/s | Kernel ms (sum) | H2D GB | D2H GB | Host copies ms | Device MB |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Robins (MSCEER lower star) | (128,128,64) | 1358 | 0.79 | 682 | 4.67 | 8.58 | 766 | 39 |
| Robins (MSCEER lower star) | (256,256,32) | **1311** | **0.82** | 634 | 4.77 | 8.58 | 649 | 79 |
| Robins (MSCEER lower star) | (1024,1024,8) | 1414 | 0.76 | 667 | 5.89 | 8.58 | 611 | 340 |
| DMT PLS | (128,128,64) | 3010 | 0.36 | **1339** | 4.67 | 8.58 | 2367 | 39 |
| DMT PLS | (256,256,32) | 2581 | 0.42 | 1438 | 4.77 | 8.58 | 1431 | 79 |
| DMT PLS | (1024,1024,8) | **2387** | 0.45 | 1465 | 5.89 | 8.58 | 773 | 340 |
| DMT PLS-BR | (128,128,64) | 2250 | 0.48 | 1452 | 4.67 | 8.58 | 1362 | 39 |
| DMT PLS-BR | (256,256,32) | **2174** | 0.49 | 1402 | 4.77 | 8.58 | 864 | 79 |
| DMT PLS-BR | (1024,1024,8) | 2323 | 0.46 | 1516 | 5.89 | 8.58 | 786 | 340 |
| DMT DMF, alg 1 only | (128,128,64) | **1752** | 0.61 | **1037** | 4.67 | 8.58 | 942 | 57 |
| DMT DMF, alg 1 only | (256,256,32) | 2297 | 0.47 | 1599 | 4.77 | 8.58 | 781 | 115 |
| DMT DMF, alg 1 + 2 | (128,128,64) | **2358** | 0.46 | **1641** | 4.80 | 8.58 | 1033 | 58 |
| DMT DMF, alg 1 + 2 | (256,256,32) | 3405 | 0.32 | 2689 | 4.93 | 8.58 | 743 | 118 |

Reading the table:
- **Robins is transfer-bound.** Its kernel is 0.63-0.68 s for 1.07 G
  vertices, but 13.3 GB has to cross PCIe and be copied on the host, so the
  wall time is about 2x the kernel time. The DMT rules are compute-bound:
  wall time is within 10% of the kernel sum, because transfers and host
  copies overlap with the kernels.
- **DMT PLS / PLS-BR** rows are after the lower-star core rewrite below
  (kernel 8.4 s -> 1.34 s, 6.3x). They now take about 2.1-2.3x the Robins
  kernel time and, like Robins, are partly limited by transfers and host
  copies (wall 2.2-2.6 s vs 1.3-1.5 s of kernel).
- **DMF** rows are after the DMF kernel rewrite below (both passes: kernel
  10.6 s -> 1.64 s, wall 11.4 s -> 2.36 s). The second pass adds about 60%
  of kernel time over the first pass alone. Smaller tiles are faster for DMF
  (1.6 s vs 2.7 s of kernel at (256,256,32)): the per-tile passes re-read
  top[] and the ranks around every cell, and a smaller local complex stays
  in L2.
- **Tile shape barely matters** (within about 10%). Small slabs (128,128,64)
  give the smallest device footprint, 39-58 MB. Thin full-plane slabs
  (1024,1024,8) stage the most halo (5.9 GB up instead of 4.7 GB).

## Lower-star core speedup: PLS / PLS-BR (2026-11-20), bit-exact

`include/dmt_core/dmt_lowerstar_pq.h` was 11x slower per vertex than the
Robins kernel (256^3 tiled: 116 ms vs 10.6 ms). `ptxas` showed a 448-byte
stack frame per thread: the packed 64-bit Robins keys, the local ranks, the
27-slot result arrays, all indexed at thread-dependent slots, plus the star
tables read from kernel-parameter space at divergent addresses. The rewrite
keeps every decision of the reference and moves the state into registers:

| Change | 256^3 PLS kernel (tiled, T=128x128x64) |
|---|---:|
| before | 116 ms |
| tables copied to shared memory; `__popc` / `__ffs`; **rank-mask keys** instead of packed keys | 33 ms |
| collapse writes through a recorder; the GPU path keeps pair codes as three bit planes (no `StarResult` in local memory, 0-byte stack) | 25 ms |
| `__launch_bounds__(256, 4)` (84 -> 64 registers, 4 blocks per SM) | **20 ms** |

At 1024^3 the kernel time went from 8.4 s to 1.34 s.

- **Rank-mask keys.** Each lower vertex gets its local rank r, and a cell's
  key is the bitmask of its vertices' ranks: 32 bits and no sort. For two
  distinct vertex sets, the Robins comparison (descending sorted values,
  lexicographic, a proper prefix smaller) is decided by the largest element
  of their symmetric difference. That element is the highest set bit of the
  XOR of the masks, so integer comparison of the masks is exactly the Robins
  order. The masks are built by unrolled loops over the canonical star (a
  cell's vertex set is its far corner plus its facets' vertex sets), so every
  index is static and the keys stay in registers. The argmin is an unrolled
  27-slot select.
- **One collapse, two recorders.** `CollapsePartT<Order, Rec>`:
  `StarResultRecorder` fills `StarResult` (partner, step, part), which the CPU
  drivers and the GPU journal path use; `StarCodeRecorder` keeps pair codes
  and a written-slot mask. `LowerStarCodes` is the journal-free entry point
  the kernels use. The in-core kernel is templated on whether journals are
  requested.
- **Tried and rejected:** a lazy comparator with no ranks, comparing cells
  by the highest-valued vertex of their symmetric difference, with per-thread
  keys in shared memory. It took 39 ms vs 25 ms.
- **The original stays as the reference.** `PLSStarPacked` /
  `PLSBRStarPacked` keep the packed 64-bit keys. The new gate 0 of
  `dmt3d_tiled_validate` compares them star by star with the new `PLSStar` /
  `PLSBRStar` and with `LowerStarCodes`. It covers 2.0 M random stars, with
  every clipping pattern (faces, edges, corners, degenerate axes) and random
  boundary masks, plus every vertex of four volumes, for both rules: **0
  differences** in partner, step, part, orphans and codes.
- Random boundary masks are outside what a real grid produces. They exposed
  a reference behaviour the code path had to match: a part that holds the
  vertex but no edge returns early and leaves its other cells unwritten. The
  code recorder tracks written slots, so the kernels write exactly the cells
  `PLSStar` writes.
- `dmt3d_validate` (CPU vs in-core GPU, 246 runs) and `dmt3d_tiled_validate`
  (440 tiled runs, halo sweep, 1024^3 footprint) still pass with 0
  mismatches.

## DMF kernel speedup (2026-11-20), bit-exact

The DMF pass kernels took 74% of DMF kernel time and the top-facet kernels
25% (Nsight Systems, 256^3). Both had the problem the PLS core had: the key
walks used arrays whose length and axes were only known at run time, so
`ptxas` put them on the stack (240 and 368 bytes per thread). The rewrite
keeps every decision of the reference:

- **Templated on the parity pattern.** The kernels are templated on the
  parity pattern P of their cells (bit i = odd coordinate on axis i), and
  the passes also on the selector. One launch covers one parity class. The
  cell dimension and its odd axes are compile-time constants, so the key
  walk (`DMFKeyWalk<D>`, since replaced by corner-rank keys, see below), the facet keys of `DMFTopFacetsP<P>` and the
  candidates of `DMFLowestSelectingP<P, W>` unroll into straight-line
  register code. All 21 instantiations have a 0-byte stack frame. The top
  facets are kept as a running largest / second-largest with static copies,
  instead of indices into a key array.
- **(cell id, corner vertex) pairs instead of coordinates.** A walk step
  along an axis moves the cell id by the cell stride, and the corner-vertex
  id (floor(q/2)) by the vertex stride when the step lands on the upper
  vertex. A leaf's rank is then a single load, with no coordinate math.
- **Algorithm 1 compares half keys.** Every algorithm-1 candidate cofacet c
  of a cell a has G0(c) = a, so all candidate keys start with the same
  F(a). Comparing them only needs the second halves, F(G1(c)), where G1(c)
  is the facet of c opposite a. Those are keys of a's dimension instead of
  one dimension higher, so the walks are half as long. Algorithm 2
  compares full keys.
- **Same code in both drivers.** The in-core driver (`ComputeDMFGradient3DGPU`,
  int64) now uses the same templated functions on `Cubical3DT<int64_t>`, in
  place of the 1D launch with id -> coordinate divisions. The (dimension,
  key) order output still reads `top[]`, so it is unchanged.

| DMF kernel time | before | after |
|---|---:|---:|
| 256^3 tiled T=(128,128,64), algorithm 1 only | 111 ms | **14 ms** (7.9x) |
| 256^3 tiled T=(128,128,64), algorithms 1 + 2 | 154 ms | **22 ms** (7.0x) |
| 512^3 tiled T=(128,128,64), 1 + 2 | 1357 ms | **195 ms** |
| 512^3 in-core, 1 + 2 | 2698 ms | **311 ms** (8.7x) |
| 1024^3 tiled T=(128,128,64), 1 + 2 | 10573 ms | **1641 ms** (6.4x) |

The id-based functions (`DMFKeyFromTop`, `DMFTopFacetsFromTop`,
`DMFLowestSelectingCofacetFromTop`) are unchanged and stay the reference. The
CPU drivers and the order sort use them. `dmt3d_validate` (CPU vs in-core
GPU, 246 runs, up to 384x256x192, plus the host key-cache check) and
`dmt3d_tiled_validate` (440 tiled runs vs in-core and CPU, halo sweep, 1024^3
footprint) pass with 0 mismatches.

## Robins transfers (2026-11-20)

At 1024^3 the Robins kernel takes 0.63 s, but end to end took 1.31-1.78 s.
Findings and changes:

- **This GPU has one copy engine** (`asyncEngineCount = 1`). Uploads and
  downloads share it and run in submission order across streams. The
  drivers queued each slab's download right behind its own upload. That
  download waits for its kernel, so every later upload sat behind it, and
  upload, kernel and download of consecutive slabs ran nearly back to back.
  **Fix:** `TransferRing::MarkComputed` / `Download` split the old
  `EndComputeAndDownload`, and the 3D drivers queue each download only after
  the next slab's upload. That is the 3D Robins driver and both DMT drivers;
  the 2D drivers are unchanged. The engine now moves one slab out while the
  next comes in, both under the kernel.
- **In-core entry point.** `ComputeDiscreteGradient3D` /
  `ComputeDiscreteGradient3DConstrained` used to stage the whole volume as a
  single tile, with no overlap; a 1024^3 volume did not fit at all. They now
  stream z-slabs of about 64 MB of cells. `ComputeDiscreteGradient3DEx` adds
  `InCoreOptions3D` with two modes:
  - copy: pageable memory, staged through the pinned ring by a persistent
    host worker pool. Spawning threads per slab copy cost more than it
    gained.
  - direct: full-XY z-slabs are contiguous in the caller's arrays, so each
    slab is one DMA straight from the caller's values and one into the
    caller's gradient. There is no host copy and no host wait per slab: slot
    reuse is ordered on the device by events, and uploads are issued
    `slots - 1` slabs ahead. Direct is chosen automatically when the arrays
    are page-locked (`IsHostMemoryPinned`). `PinHostMemory` /
    `UnpinHostMemory` let callers pin once. `hostMode = 2` pins for a single
    call, but releasing 12.9 GB of pinned memory takes about 0.7 s on this
    machine, so per-call pinning does not pay off.
- **The first GPU use of freshly pinned or allocated memory costs about
  0.4 s** (the driver maps it), so first calls are slower. The table shows
  both.
- **Packed downloads (`TiledOptions3D::packOutput`, `InCoreOptions3D::packOutput`)**:
  a device kernel packs two cells per byte (pair code plus assigned bit),
  halving the D2H bytes, and host threads expand it with a 256-entry table.
  Measured here, it does not pay. The host side is limited by memory
  bandwidth: expanding into the caller's array writes the same 8.6 GB a
  memcpy would, and costs what the smaller transfer saves (in-core copy
  1.13-1.24 s packed vs 1.06 s unpacked). It is off by default. It stays
  available, and covered by `dgrad3d_validate`, for machines where PCIe
  rather than host memory is the limit.

1024^3 Robins, steady-state calls (the first call of each mode is shown in
brackets):

| Path | Wall ms | Before |
|---|---:|---:|
| in-core, pageable arrays (copy mode) | **1064** (first 1098) | did not fit (single tile) |
| in-core, arrays pinned once (direct mode) | **879** (first 1256) | - |
| in-core, pinned for the call (`hostMode = 2`) | 1580 | - |
| tiled API with host-array callbacks, T=(256,256,32) | **1040** (best run 932) | 1311-1553 |
| tiled API, T=(128,128,64) | 1152 (best run 997) | 1358-1902 |

With the arrays pinned, the direct path reaches the copy engine's bound.
That is 13.25 GB per run at the measured ~16 GB/s H2D and ~17 GB/s D2H,
about 0.8 s, with the kernel (0.63 s) hidden underneath.

The same reordering helps the DMT drivers at 1024^3: PLS 2.1 -> 1.54 s, DMF
(both passes) 2.36 -> 1.96 s. With downloads deferred, the per-batch
"kernel" event windows include queueing, so their sums now exceed wall time;
`h2d` and `d2h` sums are still pure transfer time.

## DMT transfers (2026-11-20)

New entry points `ComputeLowerStarGradient3DStreamed` and
`ComputeDMFGradient3DStreamed` (`DmtStreamOptions3D`) take a volume of any
size in host memory, so only slabs live on the device:

- **direct** (automatic when ranks and outputs are page-locked;
  `hostMode = 2` pins for the call): full-XY z-slabs are DMAed straight
  between the caller's arrays and the device. Inside a slab the work still
  runs as (tileX, tileY) tiles of the usual tile grid. The kernels read ranks
  through the slab's strides (`TileGeom::vs`) and write owned cells at their
  slab offsets, so DMF keeps the small local complexes it needs for L2
  locality. Slot reuse is event-ordered on the device, and uploads are issued
  `slots - 1` slabs ahead, at least one.
- **copy** (pageable memory): the tiled ring driver over the host arrays.

The library's host copies (the in-core Robins adapters and the streamed copy
mode) now run on OpenMP. The project already requires it outside Apple, and
`gpu_dgrad` links `OpenMP::OpenMP_CXX`. Copies of 256 KB or more are split
across threads. The earlier 8 MB threshold left every 128x128x64 rank box
(4.7 MB) on one thread, and a condition-variable pool was too slow to wake
for ~1000 short copies per run.

1024^3, steady-state calls (first calls in brackets):

| Method | tiled ring, best tile | streamed, pageable (copy) | streamed, pre-pinned (direct) | kernel |
|---|---:|---:|---:|---:|
| PLS | 1532 | 1681 (1668) | 1610-1750 (2420) | ~1510 |
| PLS-BR | 1538 | 1721 (1663) | 1657 (1705) | ~1556 |
| DMF, alg 1 | 1259 | 1507 (2439) | 1358-1371 (1414) | ~1260 |
| DMF, alg 1 + 2 | 1854-1965 | 1905 (1893) | 1890-2023 (1997) | ~1790-1920 |

Result: the DMT rules are kernel-bound. Since the copy-engine reordering
(previous section), their 13.25 GB of transfers already overlap the kernels
fully. Wall time is within 5-10% of kernel time in every mode, so moving
fewer bytes cannot make them faster. What the direct mode does buy:

- no host copies: 0 ms of host CPU vs 0.7-0.9 s of memcpy per run;
- a volume of any size behind a plain array interface.

It costs device memory: 1.6 GB of slab buffers at slabDepth 64, vs 39-58 MB
for the tiled ring. Tile shapes tried inside slabs ((1024,32,32),
(1024,64,64), (256,128,64), (128,128,32)) land within about 10% of each
other, except wide tiles for DMF, which lose L2 locality (2.4-2.7 s).

`dmt3d_tiled_validate` covers both streamed modes for all four rules on
every volume, with tiny and default tiles and slabs, 1-3 slots, and the
single-slab case. The first version of the direct path missed that case:
with one slab the lookahead clamped to 0. That gives 664 runs with 0
mismatches; `dmt3d_validate` and `dgrad3d_validate` also pass.

## DMT kernels, second round (2026-11-20), bit-exact

Nsight Systems at 512^3 (`gpu3d_bench 512 3 4`, every DMF configuration of
the bench; GPU time summed per kernel):

| DMF kernels | before | after |
|---|---:|---:|
| pass kernels (7 patterns x 2 selectors) | 1401 ms | **1100 ms** |
| top-facet kernels (7 patterns) | 562 ms | **376 ms** |
| extract | 205 ms | **110 ms** |
| total | 2168 ms | **1587 ms** (-27%) |

What changed:

- **Keys from corner ranks.** A cell's DMF key depends only on its vertices'
  values: F_w is defined recursively, with G0 the largest-key facet and G1
  the facet opposite it. The key walk read `top[]` to find G0 at every
  level, a chain of dependent byte reads. `DMFCornerKey<P, B>` computes the
  same key from the cell's 2^D corner ranks, which load independently, by
  that recursion unrolled at compile time. Facets are taken in `Facets`
  order and the first strict maximum wins, as in `DMFKeyFromFacets`. For a
  hex's top facets this is 8 independent loads instead of 18 `top[]` and 24
  rank loads. `top[]` is still read for the selection tests. Every
  instantiation has a 0-byte stack frame.
- **Selection first.** Each (d+1)-cell selects exactly one facet, so a cell
  has about one selecting cofacet on average. The pass first collects the
  selecting cofacets (one `top[]` byte each), and computes keys only when
  there are two or more. Algorithm 1 compares only F(G1(c)), the half of
  the key that differs.
- **No `top[]` fill.** Vertex entries are never read, and every other entry
  is written before it is read.
- **Extract** handles 4 cells per thread along x and keeps per-thread
  critical counts.
- **PLS:** the census is fused into the lower-star kernel, so the census
  kernel is gone: crit and written cells are counted with popcounts of the
  written-slot masks, and unassigned = owned - written. The rank matrix
  skips the centre's 26 pairs (its rank is fixed). Together about 2-3%.
- **Tried and dropped:**
  - staging PLS ranks in shared memory: equal speed, because the kernel is
    ALU-bound in the star core;
  - a 3-blocks-per-SM register budget for PLS: 12% slower;
  - warp-ballot counters in the DMF pass: no measurable change.

1024^3 (tiled, T=(128,128,64), steady state):

| Method | before | after |
|---|---:|---:|
| DMF, alg 1 + 2 | 1854-1965 ms | **1478 ms** |
| DMF, alg 1 only | 1259 ms | **1088 ms** |
| DMF streamed, pre-pinned direct, alg 1 + 2 | 1890-2023 ms | **1659 ms** |
| DMF footprint run (`dmt3d_tiled_validate`, procedural ranks) | 2276 ms | 1466 ms |
| PLS / PLS-BR | 1532-1565 ms | 1508-1558 ms |
| in-core DMF 512^3 kernel | 311 ms | 233 ms |

DMF is now about as fast as PLS at 1024^3. PLS remains bound by the
arithmetic of its 27-slot rank matrix and priority-queue collapse.
`dmt3d_validate` (246 runs, against the CPU) and `dmt3d_tiled_validate` (664
runs, the 2.0 M-star core check, halo sweep, footprint) pass with 0
mismatches.
