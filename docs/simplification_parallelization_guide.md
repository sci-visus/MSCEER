# MSC Simplification Parallelization Guide (Low-Thread)

This guide documents how hierarchy simplification (repeated cancellations) currently works in `MSCEER`, why the serial implementation is fast, and how to safely explore approximate parallelism with a small thread count (`< 16`).

The focus is the cancellation phase in:

- `include/gi_morse_smale_complex_basic.h`
- `include/gi_morse_smale_complex_restricted.h`

## 0) Non-Negotiable Product Behavior

The hierarchy is not only a preprocessing structure. It is part of an interactive workflow where users repeatedly change persistence (for example, with a visualization slider) and expect the visible graph/manifolds to update immediately.

Because of that, this project intentionally keeps historical arcs/nodes and uses timestamped liveness instead of physically deleting/rebuilding graph structure each time.

Design consequence for future work:

- Fast persistence scrubbing (`SetSelectPersAbs(...)`) and fast iteration over living entities are first-class requirements.
- Any parallel simplification strategy that speeds up build time but degrades interactive persistence changes is not acceptable.
- Preserve append-only history plus liveness filtering semantics as a core invariant, not an implementation detail.

Exception note:

- In principle, arcs with exactly zero persistence are a special case that could be physically destroyed without harming persistence-level behavior.
- This codebase still keeps them in the historical arc arrays to preserve stable/compact index space and avoid extra complexity in bookkeeping and remapping logic.
- Unless a future refactor explicitly handles index remapping and iterator semantics, treat zero-persistence arcs the same as all others (timestamped liveness, no hard deletion).

## 1) What Simplification Is Doing

`ComputeHierarchy(pers_limit)` repeatedly cancels arcs in a mutable graph:

1. Build an initial priority queue (`edges_to_cancel`) of arc candidates.
2. Pop a candidate.
3. Re-check if it is still alive and still valid in the current graph.
4. If valid, execute `cancel(a)`.
5. `cancel(a)` destroys nodes/arcs and creates new arcs, then inserts new arcs into the queue.
6. Continue until queue exhausted or persistence threshold is reached.

This is serial because each cancellation mutates the structure used to determine validity/order of the next one.

## 2) Critical Correctness Invariant

For a cancellation to be valid, the arc must be the only currently alive arc connecting its endpoint pair.

This is enforced by `isValid()`:

- Boundary parity must match.
- Multiplicity must be exactly one at current cancellation time.

`countMultiplicity(ap, num_cancelled)` scans incident arcs at the upper node and counts alive arcs with same `(lower, upper)` pair.

Implication: any parallel strategy must preserve the "unique alive connector" rule at commit time, not just at initial candidate selection.

## 3) Core Data Model and Lifecycle

### Nodes

`node` stores:

- `firstarc` intrusive adjacency head
- `destroyed` timestamp (`INT_INFTY` means alive)
- degree counters (`numarcs`, `numlower`)
- manifold references (`amanifoldid`, `dmanifoldid`)

### Arcs

`arc` stores:

- endpoint node IDs (`lower`, `upper`)
- intrusive next pointers (`lower_next`, `upper_next`)
- lifespan (`created`, `destroyed`)
- persistence, dim, boundary
- geometry handle (`geom`)

### Liveness model

Liveness is timestamp-based, not pointer-removal-based:

- arc alive if `created <= place && destroyed > place`
- node alive if `destroyed > place`

`place` is usually `num_cancelled`.

This is central to serial speed because adjacency is never physically compacted after each cancel.

This is also central to interactive speed after the hierarchy is built: changing selected persistence should only change the liveness predicate, not trigger structural reconstruction.

## 4) Detailed Cancellation Operation (`cancel(INT_TYPE a)`)

Given candidate arc `a = (lower -> upper)`:

1. Record cancellation metadata (`mCRecords`, update `max_pers_so_far`).
2. Enumerate all alive upward arcs around `lower` (excluding `a`).
3. Enumerate all alive downward arcs around `upper` (excluding `a`).
4. Create cross-product replacement arcs using `createArc(lap, a, uap, num_cancelled + 1)`.
5. Scan lower neighborhood again; mark affected old arcs destroyed and adjust endpoint counters/manifolds.
6. Scan upper neighborhood similarly.
7. Increment global cancellation time (`num_cancelled++`).
8. Mark the canceled arc and its two endpoints destroyed at this new time.

Notes:

- New arcs are inserted into simplification queue via `InsertArcIntoSimplification(...)`.
- `cancel` resets `ap` pointer after arc creation because `arcs.push_back` may reallocate.
- Topology updates are highly local but can alter multiplicity and degree weights for nearby candidates.

## 5) Why Serial Is Fast Today

The implementation gets good serial performance from:

- **Timestamped liveness**: avoids expensive graph delete/rebuild.
- **Intrusive adjacency lists**: neighborhood scans are pointer-chasing over local incident arcs.
- **Lazy queue invalidation**: stale candidates are discarded only when popped.
- **Local mutation**: each cancellation touches only neighborhoods of two endpoints.
- **Append-only vectors**: `arcs`, `arc_merge_geoms`, `mans`, `cancel_num_to_pers`, and records mostly grow.

## 6) Ordering and Heuristics That Are Not Strict Lowest-Persistence

If the goal is exact persistence-pair ordering, the following behaviors are approximate:

1. **Countweight-based deferral in `get_next_to_cancel`**
   - Recomputes `newcountweight`.
   - If `newcountweight > se.countweight`, candidate is requeued with higher weight.

2. **Artificial persistence bump for very high countweight**
   - If `newcountweight > 1500`, does `se.persistence += 1`, then requeues if still under threshold.
   - This explicitly delays heavy candidates.

3. **Comparator includes non-persistence tie-breakers**
   - Queue ordering uses persistence, then countweight, then arc id.

4. **Threshold-limited queue insertion**
   - Base class only inserts arcs if `a.persistence <= gPersThreshold`.

Restricted variant (`gi_morse_smale_complex_restricted.h`) also changes scheduling policy (including a restricted "next cancel" path), so make sure experiments target the intended class/path.

## 7) Allocation and Contention Hotspots to Watch

If parallelized naively, hotspots include:

- `arcs.push_back(...)` and pointer invalidation/reallocation.
- `arc_merge_geoms.push_back(...)`.
- `mans.push_back(...)` from manifold merges.
- shared priority queue mutations.
- shared degree counters (`numarcs`, `numlower`) and `destroyed` timestamps.
- manifold ID rewrites (`amanifoldid`, `dmanifoldid`).

The current code is intentionally not thread-safe in these sections.

## 8) Practical Low-Thread Parallelization Launch Plan (`< 16`)

### Stage A: Keep exact commit serial, parallelize candidate preparation

Safest first step:

- Keep `cancel(a)` and global queue commit serial.
- Parallelize only candidate scoring/validation snapshots over frontier batches.
- Re-validate selected candidate in serial immediately before commit.

This yields small speedups with minimal risk.

### Stage B: Approximate batch cancellation with conflict filtering

For each batch iteration:

1. Pop top `K` candidates by current queue order (serial).
2. Snapshot lightweight metadata needed for quick conflict tests.
3. Build a conflict graph where candidates conflict if they:
   - share any endpoint node,
   - share the same endpoint pair (multiplicity hazard),
   - touch a 1-ring neighborhood that can alter each other's validity.
4. Select an independent subset (greedy by lowest persistence first).
5. Execute selected cancellations in parallel worker-local logs (no direct global writes).
6. Serial commit phase applies logs one by one with full `isValid()` recheck.
7. Rejected logs are dropped/requeued.

This is approximate because selected candidates are from a stale snapshot, but serial recheck at commit preserves graph validity.

### Stage C: Opportunistic direct parallel commit with fine-grained locks (risky)

Only consider after Stage B is stable:

- lock endpoint neighborhoods;
- re-validate under lock;
- commit if still valid.

Expected complexity is high due to manifold updates and queue mutation; deadlock avoidance and lock ordering become critical.

## 9) Invariants to Preserve During Any Experiment

Must always hold:

- No cancellation without `countMultiplicity == 1` in current live graph.
- Alive/dead timestamps remain monotonic and consistent.
- Degree counters (`numarcs`, `numlower`) match effective alive incidence.
- No arc references nodes already destroyed at the same or earlier cancellation time.
- Queue candidates are always treated as hints; full validity must be checked at execution.
- Interactive persistence switching remains lightweight: changing selected persistence must not require topology rebuild/reallocation.

## 10) Suggested Instrumentation for Experiments

Add counters/timers (guarded by macro) for:

- queue pop count, stale-pop count, requeue count
- cancellations attempted vs committed
- multiplicity failures at commit recheck
- max and average neighborhood sizes
- per-batch accepted candidate count
- time split: candidate generation vs commit

The file already uses aggregate timing helpers (chrono); extend those for batch phases.

## 11) Validation Checklist

For any parallel/approx variant:

1. Run same input/seed across serial and variant.
2. Compare:
   - final number of alive nodes/arcs at selected persistence,
   - endpoint histogram of alive arcs,
   - selected manifold outputs used by downstream consumers.
3. Assert no invariant violations in debug builds.
4. Measure speedup on thread counts: `1, 2, 4, 8, 12, 15`.
5. Track regression scenarios with dense local branching and high multiplicity pressure.

## 12) Recommended First Implementation Target

Start with **Stage A** in `MorseSmaleComplexBasic::ComputeHierarchy` path:

- keep core data structure mutations serial,
- parallelize ranking/filtering work around the queue,
- preserve existing cancellation semantics exactly at commit.

If stable, proceed to Stage B with conflict-filtered approximate batches.

---

If you are a future agent picking this up: prioritize correctness guards around multiplicity and liveness first, then optimize. The fastest way to invalidate results is to allow a cancellation whose endpoint pair is not uniquely connected at commit time.
