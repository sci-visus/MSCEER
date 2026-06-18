# Prepartitioned MSC Builder Implementation Notes

This document captures what is implemented for the geometric-first partitioned build prototype and how to validate it.

## New Components

- `include/gi_partitioned_topological_regular_grid.h`
  - Adds `PartitionedTopologicalRegularGrid2D`.
  - Preserves global mesh-style accessors (`numCellsAxis`, `numCells`, `cellid2Coords`, `coords2Cellid`, `dimension`, `boundaryValue`).
  - Adds `PartitionCellsIterator` with cached partition bounds/state:
    - `start_id`, `start_x`, `x_range`, `start_y`, `y_range`.
    - Iterates by `(x,y)` progression (no repeated `cellid -> coords` conversion per iteration step).
  - Adds `cell_id_to_partition_num(cellId)` and coordinate-based partition lookup.
  - Enforces supported partition counts:
    - `{1,2,3,4,6,8,9,12,16}`
    - hard error for unsupported counts.
  - Pinned tilings:
    - `1->1x1`, `2->1x2`, `3->1x3`, `4->2x2`, `6->2x3`, `8->2x4`, `9->3x3`, `12->3x4`, `16->4x4`.

- `include/gi_morse_smale_complex_partitioned.h`
  - Adds `MorseSmaleComplexPartitioned<...>` orchestrator.
  - Reuses `MorseSmaleComplexBasic` as the cancellation engine.
  - Adds `PartitionLocalMsc` (derived from `MorseSmaleComplexBasic`) for per-partition local instances.
  - Adds `DelayedArcRecord` for cross-partition arc deferral:
    - lower/upper cell ids,
    - lower/upper/source partition ids,
    - source local upper node id,
    - persistence hint.
  - Local build path:
    - `ComputeFromGradInPartition(...)` builds critical nodes only from partition cells and defers external-landing traces as delayed records.
  - Local cancellation restriction:
    - partition-aware `isValid(...)` override enforces interior-only cancellation.
  - Global gather/reconcile path:
    - `BuildReconciledGlobalBase(...)` gathers living local pieces,
    - lifts delayed endpoints via local lineage (`GatherNodes`) when endpoints were deleted,
    - materializes reconciled arcs in one merged global base.
  - End-to-end helper:
    - `BuildPartitionedThenContinueSerial(...)` (local partition stage + global reconcile + serial continuation).

- `msc_2d_lib/msc_2d_partitioned_smoke.cxx`
  - Added to instantiate/compile the new partitioned pipeline API and provide a parity harness scaffold.

- `msc_2d_lib/CMakeLists.txt`
  - Adds `msc_2d_partitioned_smoke` executable target.

## Validation Metrics

Use the following parity metrics when comparing partitioned-global output against a serial reference at a selected persistence level:

- living node count
- living arc count
- endpoint histogram of living arcs

These are also the fields emitted by the partitioned smoke harness.

## Commands Used

Build:

- `cmake --build build_msc2d --target msc_2d_partitioned_smoke --config Release`
- `cmake --build build_msc2d --target msc_2d_lib_smoke --config Release`

Run baseline serial smoke:

- `build_msc2d/bin/Release/msc_2d_lib_smoke.exe`

## Current Status

- New partitioned infrastructure compiles in Release (`msc_2d_partitioned_smoke` target builds).
- Existing serial baseline smoke (`msc_2d_lib_smoke`) runs successfully and provides reference output.
- Runtime parity execution path in `msc_2d_partitioned_smoke` needs additional stabilization work before treating it as authoritative; keep it as a compile/integration harness for now.
