# Prepartitioned vs Serial Timing Report

This report tracks runtime comparisons between:

- **Serial baseline:** `MorseSmaleComplexBasic` build + simplification
- **Partitioned pipeline:** partition-local build/simplify + reconcile + global simplify

## Benchmark Configuration

- Executable: `build_msc2d/bin/Release/msc_2d_partitioned_smoke.exe`
- Dataset: synthetic 2D field in `msc_2d_partitioned_smoke.cxx` (`rows=1024`, `cols=1024`)
- Partition/thread counts tested: `{1,2,4,8,16}`
- Local per-partition threshold: `localPers = 1%` of scalar range
- Global target threshold for both methods: `globalPers = 5%` of scalar range
- OpenMP environment override removed before run (`OMP_NUM_THREADS` unset)
- `SetParallelism(partitions)` is used for discrete gradient builder

## Timing Columns

- `Discrete Grad`: discrete gradient computation phase
- `Serial Construct`: serial MSC construction (`ComputeFromGrad`)
- `Serial Simplify`: serial hierarchy simplification (`ComputeHierarchy(globalPers)`)
- `Partition Local Construct`: sum of per-partition local graph construction times
- `Partition Local Simplify`: sum of per-partition local simplification times
- `Partition Reconcile`: global merge/reconcile phase
- `Partition Global Simplify`: post-reconcile simplification at `globalPers`
- `Partition Total`: `local_total + reconcile + global_simplify`

## Measured Results (ms)

| Threads | Discrete Grad | Serial Construct | Serial Simplify | Partition Local Construct | Partition Local Simplify | Partition Reconcile | Partition Global Simplify | Partition Total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 5572 | 427 | 103 | 447 | 80 | 194 | 14 | 735 |
| 2  | 2845 | 353 | 98  | 362 | 70 | 163 | 14 | 611 |
| 4  | 1499 | 319 | 98  | 300 | 87 | 164 | 15 | 571 |
| 8  | 909  | 280 | 102 | 256 | 67 | 171 | 16 | 517 |
| 16 | 683  | 267 | 99  | 240 | 42 | 161 | 17 | 475 |

## Notes

- Discrete gradient and serial construction both show expected speedup with thread count after wiring `SetParallelism(...)` to OpenMP.
- Partition reconcile is currently near-constant overhead (~160-194 ms for this input).
- In this test, partitioned total time decreases with threads, but remains above serial construct+simplify.
- Current smoke outputs also report `same_endpoint_histogram=false` at this 1%->5% workflow; timing and output parity should be interpreted separately.

## Append Policy

For future experiments, append a new section with:

1. exact code revision/patch context,
2. benchmark configuration deltas,
3. timing table,
4. brief interpretation and next action.

## Experiment 2: 2048 Multiscale Field

### Data Description

The smoke generator was changed to a larger and smoother multiscale synthetic field:

- Domain: `2048 x 2048`
- Coarse component:
  - random noise
  - smoothed with `30` averaging iterations
  - amplitude scaled by `100`
- Detail component:
  - independent random noise
  - smoothed with `2` averaging iterations
  - amplitude scaled to `< 5%` of coarse scale (implemented as `0.05 * coarseScale`)
- Final field: `coarse + detail`
- Previous ring mask logic remains disabled/commented out.

Thresholds remained:

- Local per-partition simplify: `1%` of scalar range
- Global simplify target (both serial and partitioned): `5%` of scalar range

Measured in this dataset:

- `localPers = 0.34841`
- `globalPers = 1.74205`

### Measured Results (ms)

| Threads | Discrete Grad | Serial Construct | Serial Simplify | Partition Local Construct | Partition Local Simplify | Partition Reconcile | Partition Global Simplify | Partition Total |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 6355 | 522 | 341 | 557 | 249 | 45 | 12 | 863 |
| 2  | 3352 | 419 | 329 | 477 | 229 | 44 | 14 | 766 |
| 4  | 1796 | 345 | 324 | 420 | 190 | 44 | 14 | 672 |
| 8  | 1190 | 316 | 327 | 377 | 142 | 46 | 14 | 587 |
| 16 | 1065 | 336 | 328 | 341 | 102 | 50 | 16 | 526 |

### Divergence Summary

This experiment shows clear output divergence between partitioned and serial for partition counts > 1:

- Partition `1`:
  - `same_endpoint_histogram=true`
  - exact match to serial (expected baseline).
- Partitions `2,4,8,16`:
  - `same_endpoint_histogram=false`
  - reconciled output has **more** living nodes/arcs than serial at the same global threshold.
  - divergence magnitude increases with partition count.

Observed counts from benchmark output:

- Serial (all runs): `nodes=24499`, `arcs=48685`
- Reconciled:
  - `2`: `nodes=24991`, `arcs=48997`
  - `4`: `nodes=25446`, `arcs=49314`
  - `8`: `nodes=26293`, `arcs=49800`
  - `16`: `nodes=27255`, `arcs=50426`

Interpretation:

- The frozen-cross-boundary + direct delayed-arc insertion policy is very fast for reconcile, but increasingly conservative versus serial as partition boundaries increase, leaving additional structure alive after global simplification.

## Experiment 3: Parallel Partition Execution (One Worker per Partition/Chunk)

### Code Change Summary

`BuildPartitionLocalMSCs(...)` was changed to run partition work in parallel:

- partition loop uses OpenMP parallel-for with static scheduling,
- each partition builds/cancels in thread-local local state,
- results are written once into fixed index slots (`results[pid]`),
- a single global merge/reconcile still happens afterward.

### Important Timing Interpretation

In this experiment:

- `partition_local_total` is **wall-clock elapsed** time for the parallel partition stage.
- `partition_local_construct` and `partition_local_simplify` are **summed per-partition CPU times**.

Therefore, `partition_local_total` can be much smaller than the summed construct/simplify columns.

### Measured Results (ms)

| Threads | Discrete Grad | Serial Construct | Serial Simplify | Serial Total (Wall) | Partition Local Total (Wall) | Partition Local Construct (Summed) | Partition Local Simplify (Summed) | Partition Reconcile | Partition Global Simplify | Partition Total (Wall) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1  | 6317 | 553 | 334 | 887 | 840 | 582 | 257 | 43 | 10 | 893 |
| 2  | 3369 | 429 | 333 | 762 | 390 | 515 | 252 | 44 | 13 | 447 |
| 4  | 1946 | 364 | 332 | 696 | 186 | 487 | 248 | 46 | 14 | 246 |
| 8  | 1203 | 319 | 328 | 647 | 97  | 514 | 236 | 47 | 14 | 158 |
| 16 | 1078 | 358 | 337 | 695 | 78  | 720 | 271 | 48 | 15 | 141 |

### Divergence Status

Divergence pattern remains unchanged from Experiment 2:

- partition count `1` matches serial (`same_endpoint_histogram=true`),
- partition counts `>1` diverge (`same_endpoint_histogram=false`),
- reconciled outputs continue to retain more structure than serial at the same global threshold.

### Performance Takeaway

With parallel partition execution, partitioned wall-clock time drops significantly and now outperforms serial construct+simplify for partition counts `>= 2` in this benchmark configuration.
