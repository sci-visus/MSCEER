# CPU-to-GPU Porting Guide for THIS Machine

Environment-specific guide for adding CUDA compute to a C++/CMake project on
this exact development machine. Everything here was verified while GPU-porting
MSCEER's discrete gradient (branch `cuda-gradient`, `gpu_dgrad/`): the
toolchain facts, the CMake recipe, and every bandwidth/timing constant were
measured here on 2026-08-27, not taken from spec sheets. Feed this to a new
project (e.g. diffg) as-is.

---

## 1. The machine

| Component | Value |
|---|---|
| GPU | NVIDIA GeForce RTX 5070 Laptop, **sm_120** (Blackwell), 36 SMs, 8.55 GB VRAM |
| Driver | 596.08 (CUDA driver API 13.2) |
| CUDA toolkits installed | **v13.0 (use this one)** and v12.1 at `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\` |
| Host compiler | Visual Studio 2022 Professional 17.14, MSVC 19.44 |
| CMake | 4.2.0 |
| Generators | Visual Studio 17 2022 (VS CUDA 13.0 build customizations ARE installed). No ninja on PATH. |
| Profilers | **Nsight Compute is NOT installed** — use `ptxas -v` for register/spill stats (§7) |
| PCIe | Gen5, **trains at x8** (x16 physical) — this halves the transfer numbers you'd expect |

### Trap #1 (will burn you first): the PATH `nvcc` is the WRONG one
`nvcc` on PATH is **CUDA 12.1**, which cannot compile for sm_120 (Blackwell
needs CUDA >= 12.8). CUDA 13.0 is installed alongside it. Never rely on PATH:

- Visual Studio generator: select the toolkit with **`-T cuda=13.0`** (this is
  how the VS generator picks toolkits; `CMAKE_CUDA_COMPILER` is ignored there).
- Ninja/Makefiles (if you install ninja): pass
  `-DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.0/bin/nvcc.exe"`.

---

## 2. CMake recipe (proven configuration)

Keep CUDA strictly opt-in so non-CUDA builds of the project are untouched:

```cmake
# root CMakeLists.txt  (project(... LANGUAGES CXX) stays as-is)
option(GPU_ENABLED "Build CUDA compute path" FALSE)
if(GPU_ENABLED)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES native)   # resolves to 120 on this machine
    endif()
    enable_language(CUDA)
    add_subdirectory(gpu)          # all CUDA code lives in its own directory
endif()
```

```cmake
# gpu/CMakeLists.txt
if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.8)
    message(WARNING "CUDA ${CMAKE_CUDA_COMPILER_VERSION} cannot target sm_120; use -T cuda=13.0")
endif()
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

add_library(gpu_compute STATIC api.cu kernels.cu api.h)
target_include_directories(gpu_compute PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(gpu_validate validate.cxx)              # host-only C++, compiled by MSVC
target_link_libraries(gpu_validate PRIVATE <your-cpu-lib> gpu_compute)
```

Configure + build + run (verbatim, adapted paths):

```bash
cmake -S . -B build_cuda -G "Visual Studio 17 2022" -A x64 -T cuda=13.0 -DGPU_ENABLED=ON
```

```bash
cmake --build build_cuda --config Release --target gpu_validate
```

Facts about this configuration, all verified:
- CMake identifies "NVIDIA 13.0.88 with host compiler MSVC 19.44" and it just
  works; a full configure takes ~30 s (CUDA detection), ~5 s without CUDA.
- **`LNK4098: defaultlib 'LIBCMT' conflicts...`** appears on every CUDA exe:
  it is the static cudart (nvcc default) meeting MSVC `/MD` host objects.
  Benign — keep static cudart deliberately, because the CUDA v13 `bin` dir is
  NOT on PATH here, so a shared-cudart exe would fail to find
  `cudart64_13.dll` when run outside a dev shell.
- CUDA 13 removed deprecated `cudaDeviceProp` fields (`clockRate`,
  `memoryClockRate`, `memoryBusWidth`). Query name / `multiProcessorCount` /
  `totalGlobalMem` / `major`/`minor` only, or it won't compile.
- MSVC host flags pass through nvcc cleanly in this setup; no
  `/permissive-`-style conflicts were hit.

### Architectural rule that saved the MSCEER port
**Do not compile your existing C++ header stack under nvcc.** Structure:
- `api.h` — plain C++ header (no CUDA types), the only thing host code sees.
  Contract: raw pointers + dims in, raw output buffer out.
- `kernels.cu` — standalone device code. It re-implements the inner algorithm
  against flat arrays; it includes nothing from the host library.
- `validate.cxx` — host-only translation unit compiled by MSVC. THIS is where
  your existing library headers get included, the CPU reference runs, and
  outputs are compared.

If the kernel needs lookup tables / constants that the CPU library owns,
**generate them at runtime in the validator from the CPU library's own data
structures** and pass them to the GPU launcher (upload to `__constant__`
memory). One source of truth, and the comparison harness proves agreement.

---

## 3. Measured bandwidth & latency constants (use these, not spec sheets)

| Constant | Measured value | Note |
|---|---|---|
| Device DRAM, effective | **338 GB/s** | D2D copy, read+write counted; spec is 448 |
| PCIe H2D, **pinned** (`cudaMallocHost`) | **20.9 GB/s** | Gen5 x8 |
| PCIe D2H, **pinned** | **20.1 GB/s** | |
| PCIe H2D, **pageable** (`std::vector`/`new`) | **~11–12 GB/s** | measured in a real harness |
| PCIe D2H, **pageable** | **~9 GB/s** | pageable D2H is the worst path |
| First CUDA call (context creation) | **~150–300 ms, once** | do a warm-up call before timing |
| Kernel launch overhead | ~10–20 µs | irrelevant above ~1 ms kernels |
| `cudaMalloc`/`cudaFree` | sub-ms each | fine per-call; reuse buffers inside loops |
| Usable VRAM budget | plan for **~7 GB** of the 8.55 | |

Practical consequences:
- Pinned memory is worth exactly **~2x** on transfers here. For anything
  transfer-dominated, `cudaMallocHost` the staging buffers.
- H2D and D2H are separate copy engines: with pinned memory + 2–3 streams you
  can overlap upload, compute, and download; the pipeline then costs
  ~max(t_h2d, t_kernel, t_d2h) instead of the sum.
- PCIe x8 means this machine is more transfer-limited than a desktop. The
  single highest-leverage design decision is **keeping data resident on the
  GPU across processing steps** and downloading compacted results.

---

## 4. Timing expectations: formulas and worked examples

For one offloaded step:

```
t_end_to_end = t_h2d + t_kernel + t_d2h
t_h2d  = input_bytes  / 21 GB/s   (pinned)   or / 11 GB/s (pageable)
t_d2h  = output_bytes / 20 GB/s   (pinned)   or /  9 GB/s (pageable)
t_kernel >= traffic_bytes / 338 GB/s          (bandwidth floor)
```

The kernel floor only binds for streaming/stencil kernels. **Measured
reality: an ALU-heavy or divergent kernel on this GPU lands at 10–100 GB/s
effective, not 338.** (MSCEER's per-vertex combinatorial kernel: ~15 GB/s
effective, entirely ALU/latency-bound — and still 129x the CPU.)

### Worked example (measured end-to-end, MSCEER 2D gradient, 4096x4096)
67 MB float input up, 67 MB byte output down, pageable buffers:

| | time |
|---|---|
| CPU (16 threads, OpenMP) | 1306 ms |
| GPU H2D | 6.1 ms |
| GPU kernel | 10.2 ms |
| GPU D2H | 7.3 ms |
| **GPU end-to-end** | **23.5 ms — 56x** |

### Expectation table for common shapes (pinned transfers)

| Working set (in + out) | Transfer cost | GPU worth it when CPU step costs |
|---|---|---|
| 10 MB | ~1 ms | > ~5 ms |
| 100 MB | ~10 ms | > ~30 ms |
| 1 GB | ~100 ms | > ~300 ms |
| 6 GB (VRAM-resident cap) | ~600 ms | > ~2 s, or amortized across steps |

Rules of thumb calibrated to this machine:
- **Break-even:** don't offload a step whose CPU time is under ~2–3x its
  round-trip transfer cost — unless the data can stay resident for later steps,
  which changes everything (transfer paid once, every subsequent step is
  kernel-only).
- **Speedup ceilings:** memory-bandwidth-bound work: ~5–10x vs a well-threaded
  CPU (338 vs ~40–60 GB/s effective CPU). Compute-bound, massively parallel,
  cache-unfriendly-on-CPU work: **50–150x** (measured 129x). Transfer-dominated
  work: can be ~1x or worse — restructure for residency instead.
- A per-element independent kernel on this GPU sustains on the order of
  **1–2 x 10^9 nontrivial work-items per second** even when heavily divergent
  (measured: 16.7M lower-star expansions, each a data-dependent iterative
  loop, in ~10 ms).
- Larger than 7 GB working set: process in chunks with double-buffered pinned
  staging; throughput degrades to PCIe rate (~20 GB/s of payload).

### Measurement pattern (avoids the classic lies)
- Warm up: run the full call once, discard; report the second run
  (first pays ~150–300 ms context init + JIT).
- Time GPU phases with `cudaEvent_t` pairs around h2d / kernel / d2h;
  time CPU with `std::chrono::steady_clock`.
- Verify before timing — a wrong kernel is often a fast kernel.

---

## 5. Porting method that worked (stage it exactly like this)

1. **Stage 0 — scaffolding + smoke test.** Opt-in CMake as in §2; a smoke
   test that queries the device, round-trips a trivial kernel, and measures
   the §3 bandwidth table on YOUR buffers. Commit before any algorithm work.
   (Copy `MSCEER/gpu_dgrad/smoke_test.cu` — it is self-contained.)
2. **Stage 1 — simplest correct kernel + validation harness.** Naive
   one-thread-per-item kernel, global memory only, no tuning. The harness
   runs the existing CPU code and the kernel on identical inputs and compares
   outputs elementwise, printing the first few mismatches with coordinates.
   Get this to zero/tolerance BEFORE optimizing anything.
3. **Stage 2+ — optimize with the harness as a ratchet.** Shared-memory
   tiling, fusing precomputation passes into the kernel (often the biggest
   structural win — MSCEER's fusion deleted a 474 ms CPU prerequisite pass
   and 2/3 of the upload), batching many small problems into one launch.
   Keep the naive kernel callable behind a flag for A/B.

Validation policy:
- **Integer/combinatorial/comparison-only algorithms: demand bit-exact
  (`memcmp`) equality with the CPU.** Comparisons and integer ops are
  deterministic on GPU; any mismatch is a real bug (usually an iteration-order
  or tie-break difference — replicate the CPU's exact orders).
- Floating-point *arithmetic* results: expect last-ulp differences (FMA
  contraction, reduction order). Compare with a relative tolerance, and
  compile with fast-math OFF (default) while validating.
- Ties/orderings that depend on comparing equal floats: make the CPU's
  tie-break rule (usually index-based) explicit and replicate it exactly.

Kernel-side habits that mattered:
- One thread per independent work item; tolerate divergence in small
  data-dependent loops (it cost far less than expected).
- Compose each output word in a register and store it **once** — never port
  CPU read-modify-write helper patterns.
- Small dynamic-indexed per-thread arrays land in L1-resident local memory —
  that's fine; check `ptxas -v` shows **0 spill bytes** and < ~64 registers.
- `__constant__` memory for small shared tables; shared memory for tiles
  (a 16x16 block with 1-halo tile costs ~1.4 KB — nowhere near limits).

---

## 6. Pitfall checklist (all hit or verified on this machine)

- [ ] PATH `nvcc` is CUDA 12.1 → always `-T cuda=13.0` (VS) or explicit
      `CMAKE_CUDA_COMPILER` (ninja). sm_120 needs >= 12.8.
- [ ] `CMAKE_CUDA_ARCHITECTURES native` must be set BEFORE
      `enable_language(CUDA)`.
- [ ] LNK4098 LIBCMT warning: benign; keep static cudart (no CUDA bin on PATH).
- [ ] Don't read removed `cudaDeviceProp` fields (clockRate etc.) — CUDA 13.
- [ ] Don't include the host template/library headers in `.cu` files.
- [ ] Pageable D2H is the slowest path (~9 GB/s); pin buffers that matter.
- [ ] First-call context cost (~0.2 s) — warm up before benchmarking, and
      don't let it poison a "GPU is slow for small inputs" conclusion.
- [ ] 8 GB VRAM: budget ~7 GB; chunk + stream beyond that.
- [ ] No Nsight Compute installed; register/spill stats via §7.

## 7. Profiling without Nsight (works today)

```bash
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.0\bin\nvcc.exe" -arch=sm_120 -std=c++17 -O3 --ptxas-options=-v -c gpu\kernels.cu -o NUL -ccbin "C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.44.35207\bin\HostX64\x64"
```

Read from the output, per kernel: `Used N registers` (keep < ~64 for good
occupancy at 256–512 threads/block), `spill stores/loads` (want 0), `smem`
per block. For memory-boundedness, compare achieved GB/s (your traffic model
/ measured kernel ms) against 338 GB/s; below ~100 GB/s with real traffic you
are ALU/latency-bound and should look at divergence and redundant recompute,
not at memory. Installing Nsight Compute (free download) is worthwhile the
moment kernel tuning becomes a real work item.

---

*Provenance: all numbers measured 2026-08-27 on this machine via
`MSCEER/gpu_dgrad` (smoke test + Stage 1/2 validation harness), CUDA 13.0.88,
driver 596.08. Reference implementation of every pattern above:
`MSCEER/gpu_dgrad/{CMakeLists.txt, dgrad_gpu_api.h, dgrad2d_kernel.cu,
dgrad2d_validate.cxx, smoke_test.cu}`.*
