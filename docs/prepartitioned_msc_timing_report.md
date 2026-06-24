# Prepartitioned vs Serial Timing Report

This report is a full reset. All previous benchmark sections were removed because they were based on incomplete runs.

## Run Configuration

- Date: 2026-06-24
- Executable: `build_msc2d/bin/Release/msc_2d_partitioned_smoke.exe`
- Data size:
  - Grid: `2048 x 2048` samples
  - Total scalar samples: `4,194,304`
- Data description:
  - Synthetic multiscale scalar field generated in `msc_2d_partitioned_smoke.cxx`
  - Coarse component: random noise smoothed for `30` iterations, scaled by `100`
  - Detail component: independent random noise smoothed for `2` iterations, scaled by `0.05 * coarseScale`
  - Final field: `coarse + detail` (ring mask disabled)
- Partition/thread counts tested: `{1, 2, 4, 8, 16}`
- Local threshold: `localPers = 0.34841` (`1%` of scalar range)
- Global threshold: `globalPers = 1.74205` (`5%` of scalar range)
- Command pattern:
  - `build_msc2d/bin/Release/msc_2d_partitioned_smoke.exe <partitions>`

## Timing Results (ms)

`partition_local_total` and `partition_total` are wall-clock times.  
`partition_local_construct` and `partition_local_simplify` are summed per-partition times.

| Partitions | Discrete Grad | Serial Construct | Serial Simplify | Partition Local MSC | Partition Local Construct (Sum) | Partition Local Simplify (Sum) | Partition Reconcile | Partition Global Simplify | Partition Total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 6238 | 599 | 340 | 945 | 616 | 329 | 55 | 12 | 1012 |
| 2  | 3330 | 419 | 333 | 446 | 541 | 333 | 56 | 16 | 518 |
| 4  | 1830 | 363 | 331 | 212 | 515 | 319 | 57 | 16 | 285 |
| 8  | 1361 | 355 | 339 | 151 | 663 | 349 | 62 | 16 | 229 |
| 16 | 1072 | 324 | 336 | 99  | 701 | 356 | 62 | 22 | 183 |

## Summary Timings Table (ms)

- `Serial Simplified MSC = Serial Construct + Serial Simplify`
- `Partitioned Simplified MSC = Partition Local MSC + Partition Reconcile + Partition Global Simplify`

| Partitions | Discrete Grad | Serial Simplified MSC | Partitioned Simplified MSC |
|---:|---:|---:|---:|
| 1  | 6238 | 939 | 1012 |
| 2  | 3330 | 752 | 518 |
| 4  | 1830 | 694 | 285 |
| 8  | 1361 | 694 | 229 |
| 16 | 1072 | 660 | 183 |

## Structure/Parity Snapshot

| Partitions | Serial Nodes | Serial Arcs | Reconciled Nodes | Reconciled Arcs | Same Endpoint Histogram |
|---:|---:|---:|---:|---:|:---|
| 1  | 24499 | 48685 | 24499 | 48685 | true |
| 2  | 24499 | 48685 | 24499 | 48685 | false |
| 4  | 24499 | 48685 | 24497 | 48681 | false |
| 8  | 24499 | 48685 | 24495 | 48677 | false |
| 16 | 24499 | 48685 | 24495 | 48677 | false |

## Freeze-Barrier Diagnostics

These are from the new global freeze-exchange path (`pre_exchange` is local-only freeze count, `post_exchange` is after global intent fold/apply):

| Partitions | Pre-Exchange Frozen Nodes (Sum) | Post-Exchange Frozen Nodes (Sum) | Pre-Exchange Max per Partition | Post-Exchange Max per Partition |
|---:|---:|---:|---:|---:|
| 1  | 0    | 0    | 0   | 0   |
| 2  | 833  | 1570 | 448 | 788 |
| 4  | 1700 | 3129 | 463 | 816 |
| 8  | 3376 | 6159 | 604 | 1001 |
| 16 | 5074 | 9317 | 437 | 834 |

## Summary

- This benchmark is the new baseline after the freeze-barrier refactor.
- Partitioned wall-clock time improves strongly as partition count increases.
- Reconciled structure is identical at `partitions=1` and remains very close (but not identical) for `partitions>1`.
