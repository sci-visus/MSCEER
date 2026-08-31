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
