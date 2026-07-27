# Prepartitioned MSC Builder Implementation Notes

This document describes the current implemented partitioned approach, including the freeze-barrier refactor and the public consumers (`msc_2d_lib` and `msc_py`).

## Current Architecture

### Partition grid and supported parallelism

- `include/gi_partitioned_topological_regular_grid.h`
  - `PartitionedTopologicalRegularGrid2D` provides partition indexing and per-partition iterators.
  - `PartitionedTopologicalRegularGrid3D` provides 3-axis partition indexing and per-partition iterators.
  - 3D partition specification supports:
    - explicit tuple splits `(sx, sy, sz)`
    - single `requested_partitions` with ratio-aware split inference.
  - 3D split inference enumerates exact factor triples and selects the minimum distortion candidate against axis-size ratios (deterministic tie-breaks).
  - Supported partition counts are fixed to:
    - `{1,2,3,4,6,8,9,12,16}`
  - `msc_2d_lib` clamps requested parallelism down to the nearest supported value.

### Partitioned pipeline

- `include/gi_morse_smale_complex_partitioned.h`
  - `MorseSmaleComplexPartitioned<...>` orchestrates local build, global reconcile, and optional global simplify.
  - `PartitionLocalMsc` performs per-partition construction/simplification with partition-aware validity checks.
  - `DelayedArcRecord` stores cross-partition arc endpoints for reconcile.
  - `LineageTransferRecord` carries local constituent lineage for fast manifold remap in reconciled output.
  - Scheduling is decoupled from partition count:
    - build/simplify loops use explicit `num_threads`
    - partitions are assigned round-robin via OpenMP `schedule(static, 1)`.

## Freeze-Barrier Refactor (Implemented)

`BuildPartitionLocalMSCs(...)` now uses explicit phases:

1. **Phase 1: Build local graphs**
   - `ComputeFromGradInPartition(...)` runs for each partition.
   - Cross-partition detections emit delayed arcs and freeze intents.
2. **Phase 2: Global freeze exchange barrier**
   - Freeze intents are accumulated globally using per-partition/per-writer buckets.
   - Intents are folded, deduplicated, and applied to each partition before local simplify.
3. **Phase 3: Local simplify**
   - `ConfigurePartitionRestrictions(..., true)` then `ComputeHierarchy(local_persistence_abs)`.
   - Local lineage export is captured for reconcile remap.

This guarantees freeze synchronization before local cancellation and supports arcs crossing many partitions.

## `msc_2d_lib` Integration

- `msc_2d_lib/msc_2d_lib.h` exposes:
  - `Msc2D::BuilderMode { Serial, Partitioned }`
  - `Msc2D::ComputeOptions` with:
    - `builderMode`
    - `requestedParallelism`
    - `basePersistenceAbs`
    - `cancelPersistenceAbs`
    - `accurateAsc`, `accurateDsc`
- Default `ComputeOptions` behavior remains serial unless caller selects partitioned mode.
- Partitioned manifold label paths in `msc_2d_lib/msc_2d_lib.cxx` use imported lineage for fast remap.

## `msc_py` Integration (New approach)

- `msc_py` now routes computation through `msc_2d_lib::Msc2D` (not direct `MorseSmaleComplexBasic` construction).
- `ComputeMSC(...)` defaults:
  - `builderMode = Partitioned`
  - `base_persistence_abs = 0.0` (strict zero local simplify unless caller overrides)
  - `num_threads = 8`
- Special case:
  - if `num_threads == 1`, `ComputeMSC(...)` forces `BuilderMode::Serial`.
- Other binding APIs (`GetAsc2Manifolds`, `GetDsc2Manifolds`, `GetGraph`, `GetCriticalPoints`) now read from `Msc2D`.

## `msc_3d_lib` Integration (Phase 1)

- New module:
  - `msc_3d_lib/msc_3d_lib.h`
  - `msc_3d_lib/msc_3d_lib.cxx`
  - `msc_3d_lib/CMakeLists.txt`
- New public API:
  - `Msc3D::BuilderMode { Serial, Partitioned }`
  - `Msc3D::ComputeOptions` with:
    - `requestedPartitions` or explicit `(splitX, splitY, splitZ)`
    - independent `numThreads`
    - `basePersistenceAbs` and `cancelPersistenceAbs`
- New smokes:
  - `msc_3d_lib_smoke`
  - `msc_3d_partitioned_smoke`

### Current 3D limitations (tracked)

- Phase 1 focuses on structural pipeline parity and diagnostics; 3D manifold-label remap parity is still deferred.
- The current 3D path uses a deterministic explicit mesh-function + all-critical gradient bootstrap to keep the partitioned 3D local-build -> freeze-barrier -> reconcile -> simplify flow runnable.
- Single-partition structural metrics in `msc_3d_partitioned_smoke` are reported and compared, but currently emitted as warnings when mismatched instead of hard-failing.

### 3D smoke matrix snapshot

- `msc_3d_partitioned_smoke 1 3`
  - inferred splits `(1,1,1)`, freeze sums `0 -> 0`, single-partition parity currently reports warning-level mismatch.
- `msc_3d_partitioned_smoke 6 4`
  - inferred splits `(3,2,1)`, freeze sums `6669 -> 13182`, round-robin assignment `t0=2 t1=2 t2=1 t3=1`.
- `msc_3d_partitioned_smoke 12 5`
  - inferred splits `(3,2,2)`, freeze sums `10211 -> 19932`, round-robin assignment `t0=3 t1=3 t2=2 t3=2 t4=2`.

## Validation Metrics

When comparing partitioned and serial output, use:

- living node count
- living arc count
- endpoint histogram of living arcs
- ascending/descending manifold unlabeled pixel counts
- freeze exchange diagnostics:
  - pre/post frozen-node counts per partition

## Commands

Build smoke targets:

- `cmake --build build_msc2d --target msc_2d_partitioned_smoke --config Release`
- `cmake --build build_msc2d --target msc_2d_lib_smoke --config Release`
- `cmake --build build_msc2d --target msc_3d_partitioned_smoke --config Release`
- `cmake --build build_msc2d --target msc_3d_lib_smoke --config Release`

Enable and build Python module:

- `cmake -S . -B build_msc2d -DMSC_PY=ON -DMSC_2D_LIB=ON`
- `cmake --build build_msc2d --target msc_py --config Release`

Enable and build 3D library:

- `cmake -S . -B build_msc2d -DMSC_2D_LIB=ON -DMSC_3D_LIB=ON`
- `cmake --build build_msc2d --target msc_3d_lib_smoke msc_3d_partitioned_smoke --config Release`

Python smoke check:

- import `msc_py`
- call `ComputeMSC(msc_id, raw)` and verify partitioned logs
- call `ComputeMSC(msc_id, raw, num_threads=1)` and verify serial logs
