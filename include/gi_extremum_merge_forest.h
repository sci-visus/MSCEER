#ifndef GI_EXTREMUM_MERGE_FOREST_H
#define GI_EXTREMUM_MERGE_FOREST_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "gi_basic_types.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace GInt {

	// Merge forest over a fixed set of extrema (minima or maxima), the
	// non-destructive replacement for UFMergeGraph::SimplifyToThreshold.
	//
	// Arcs (saddles connecting two extrema) are processed in LOWEST-PERSISTENCE
	// -FIRST order, with the same lazy recompute-and-repush as
	// UFMergeGraph::PerformCancel: when a popped arc's endpoints have been
	// re-parented, its persistence is recomputed against the current
	// representatives and it goes back on the heap. The survivor of a merge is
	// the more extreme representative under the same (value, tieKey, cell)
	// total order as TopologicalMaxVertexMeshFunction::lessThan, matching
	// MergeByVal. Each merge records, for the DYING representative, its
	// surviving parent and the merge persistence |f(saddle) - f(dying)|.
	//
	// This order matters: processing by saddle VALUE instead (the textbook
	// join-tree elder rule) yields the same number of survivors but attaches
	// dying regions to different survivors, and measured only ~80% pixel
	// agreement with the MSC hierarchy versus 99.6% for persistence order.
	// The MSC's own cancellation loop is lowest-persistence-first, so that is
	// the semantics this class reproduces.
	//
	// The difference from the legacy class is WHEN the threshold is applied:
	// legacy stops the drain at the threshold and must re-drain from scratch
	// for another one, while this runs the heap to completion once and records
	// the whole forest. Because the heap pops in nondecreasing persistence and
	// each node dies exactly once, persistence is NONDECREASING along every
	// parent chain, so the living representative of node i at threshold t is
	// found by climbing while mergePers <= t - ONE build answers every
	// threshold.
	//
	// Partitioned build: arcs may be supplied as per-partition "local" lists
	// (both endpoints owned by the arc's partition) plus a "cross" list, with
	// every cross-arc endpoint pre-flagged frozen. Local lists are sorted and
	// merged in parallel; a component containing a frozen node is tainted and
	// never merges locally - a local arc touching a tainted component taints
	// the other side too ("taint-on-contact") and is deferred. This makes each
	// component's applied merges exactly a sorted prefix of its arc stream, so
	// replaying the deferred + cross residue in sorted order afterwards is
	// order-consistent with the fully serial build. Worst case (everything
	// tainted) degenerates to the serial sorted build - correct, just slower.
	class ExtremumMergeForest {
	public:
		struct Arc {
			int32_t e1, e2;         // compact node ids (current reps once popped)
			float saddleValue;      // f(saddle) under the mesh function
			INDEX_TYPE saddleCell;  // unique tie-break + diagnostics
			float pers;             // heap key; recomputed on lazy re-push
		};

		struct BuildStats {
			double sortLocalMs = 0, mergeLocalMs = 0, residueSortMs = 0, residueMergeMs = 0;
			long long numArcs = 0, numCross = 0, numDeferred = 0, numLoops = 0, numMerges = 0;
		};

		// n nodes; cells[i] = extremum cell id (the output identity),
		// values[i] = f(extremum), tieKeys[i] = Cell2HighestVertex(cells[i])
		// (== cells[i] for vertex extrema). minDirection: true = min graph
		// (lower survives), false = max graph (higher survives).
		void SetNodes(int32_t n, const INDEX_TYPE* cells, const float* values,
			const INDEX_TYPE* tieKeys, bool minDirection) {
			mN = n;
			mMinDirection = minDirection;
			mNodeCell.assign(cells, cells + n);
			mNodeValue.assign(values, values + n);
			mNodeTieKey.assign(tieKeys, tieKeys + n);
		}

		// localArcs: numPartitions vectors, each sorted and consumed in place;
		// every arc in localArcs[p] must have BOTH endpoints owned by p.
		// crossArcs: arcs spanning partitions (any order). frozenInit: mN
		// bytes, 1 = node is an endpoint of some cross arc (BOTH endpoints of
		// every cross arc must be 1 - required for correctness, not an
		// optimization). Serial baseline: numPartitions = 1, empty crossArcs,
		// all-zero frozenInit.
		void Build(int numPartitions, std::vector<Arc>* localArcs,
			std::vector<Arc>& crossArcs, const std::vector<uint8_t>& frozenInit,
			BuildStats* stats) {
			mUf.resize(mN);
			mParent.resize(mN);
			mMergePers.assign(mN, std::numeric_limits<float>::infinity());
			for (int32_t i = 0; i < mN; ++i) { mUf[i] = i; mParent[i] = i; }
			mTaint = frozenInit;
			mResolved.assign(mN, -1);
			mStamp.assign(mN, 0u);
			mCurStamp = 0;

			BuildStats local;
			for (int p = 0; p < numPartitions; ++p) local.numArcs += (long long)localArcs[p].size();
			local.numArcs += (long long)crossArcs.size();
			local.numCross = (long long)crossArcs.size();

			std::vector<std::vector<Arc>> deferred((size_t)numPartitions);
			// per-partition counters instead of OpenMP reductions (MSVC's
			// OpenMP 2.0 reduction support is narrow; these are free anyway)
			std::vector<long long> pLoops((size_t)numPartitions, 0);
			std::vector<long long> pMerges((size_t)numPartitions, 0);
			std::vector<long long> pDefers((size_t)numPartitions, 0);
			std::vector<double> pSortMs((size_t)numPartitions, 0.0);
			std::vector<double> pMergeMs((size_t)numPartitions, 0.0);

			// Per-partition sorts ARE the parallel sort. Each partition's arcs
			// touch only partition-owned nodes, so UF state, taint bytes and
			// merge records are partition-disjoint - no synchronization.
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1)
#endif
			for (int p = 0; p < numPartitions; ++p) {
				std::vector<Arc>& arcs = localArcs[p];
				const auto tS = std::chrono::steady_clock::now();
				// seed the heap: persistence w.r.t. the arc's own endpoints,
				// exactly as UFMergeGraph::AddArc computes it
				for (size_t i = 0; i < arcs.size(); ++i) arcs[i].pers = initialPers(arcs[i]);
				std::priority_queue<Arc, std::vector<Arc>, ArcWorse> heap(
					ArcWorse(), std::move(arcs));
				arcs.clear();
				pSortMs[(size_t)p] = msSince(tS);

				const auto tM = std::chrono::steady_clock::now();
				while (!heap.empty()) {
					Arc a = heap.top();
					heap.pop();
					const int32_t r1 = find(a.e1);
					const int32_t r2 = find(a.e2);
					if (r1 == r2) { pLoops[(size_t)p]++; continue; } // loop-skip BEFORE taint check
					if (a.e1 != r1 || a.e2 != r2) {
						// stale: recompute against the current representatives
						// and re-push (legacy PerformCancel's lazy update)
						a.e1 = r1; a.e2 = r2;
						a.pers = currentPers(a, r1, r2);
						heap.push(a);
						continue;
					}
					if (mTaint[(size_t)r1] | mTaint[(size_t)r2]) {
						mTaint[(size_t)r1] = 1;          // taint-on-contact
						mTaint[(size_t)r2] = 1;
						deferred[(size_t)p].push_back(a);
						pDefers[(size_t)p]++;
						continue;
					}
					applyMerge(r1, r2, a);
					pMerges[(size_t)p]++;
				}
				pMergeMs[(size_t)p] = msSince(tM);
			}
			long long defers = 0;
			for (int p = 0; p < numPartitions; ++p) {
				local.numLoops += pLoops[(size_t)p];
				local.numMerges += pMerges[(size_t)p];
				defers += pDefers[(size_t)p];
				// wall-clock proxy: the slowest band gates the parallel phase
				if (pSortMs[(size_t)p] > local.sortLocalMs) local.sortLocalMs = pSortMs[(size_t)p];
				if (pMergeMs[(size_t)p] > local.mergeLocalMs) local.mergeLocalMs = pMergeMs[(size_t)p];
			}
			local.numDeferred = defers;

			// Global residue: deferred + cross, sorted with the SAME comparator
			// (tie-break included - deferred arcs from different partitions can
			// share a saddle value on quantized data), merged with no taint.
			const auto tRes = std::chrono::steady_clock::now();
			std::vector<Arc> residue;
			residue.reserve((size_t)defers + crossArcs.size());
			for (int p = 0; p < numPartitions; ++p)
				residue.insert(residue.end(), deferred[(size_t)p].begin(), deferred[(size_t)p].end());
			for (size_t i = 0; i < crossArcs.size(); ++i) {
				Arc a = crossArcs[i];
				a.pers = initialPers(a);
				residue.push_back(a);
			}
			std::priority_queue<Arc, std::vector<Arc>, ArcWorse> heap(
				ArcWorse(), std::move(residue));
			local.residueSortMs = msSince(tRes);

			const auto tResM = std::chrono::steady_clock::now();
			while (!heap.empty()) {
				Arc a = heap.top();
				heap.pop();
				const int32_t r1 = find(a.e1);
				const int32_t r2 = find(a.e2);
				if (r1 == r2) { local.numLoops++; continue; }
				if (a.e1 != r1 || a.e2 != r2) {
					a.e1 = r1; a.e2 = r2;
					a.pers = currentPers(a, r1, r2);
					heap.push(a);
					continue;
				}
				applyMerge(r1, r2, a);
				local.numMerges++;
			}
			local.residueMergeMs = msSince(tResM);

			if (stats) *stats = local;
		}

		// outRep[i] = compact id of i's living representative at threshold t
		// (climb while mergePers <= t; roots are +inf). Memoized per call so
		// each chain is walked once - O(n) amortized. NOT thread-safe (shared
		// scratch); the min and max forests are separate objects, so the two
		// directions may still run concurrently.
		// Parallel variant if ever needed: pointer doubling on a scratch
		// parent array p'[i] = (pers[i] <= t ? parent[i] : i), race-tolerant
		// like Label2DFlowTerminals - deferred until a profile demands it.
		void BuildRemap(float t, std::vector<int32_t>& outRep) {
			outRep.resize((size_t)mN);
			if (++mCurStamp == 0) { std::fill(mStamp.begin(), mStamp.end(), 0u); mCurStamp = 1; }
			for (int32_t i = 0; i < mN; ++i) {
				int32_t cur = i;
				mPath.clear();
				while (mStamp[(size_t)cur] != mCurStamp && mMergePers[(size_t)cur] <= t) {
					mPath.push_back(cur);
					cur = mParent[(size_t)cur];
				}
				const int32_t root = (mStamp[(size_t)cur] == mCurStamp) ? mResolved[(size_t)cur] : cur;
				mResolved[(size_t)cur] = root;
				mStamp[(size_t)cur] = mCurStamp;
				for (size_t j = 0; j < mPath.size(); ++j) {
					mResolved[(size_t)mPath[j]] = root;
					mStamp[(size_t)mPath[j]] = mCurStamp;
				}
				outRep[(size_t)i] = root;
			}
		}

		int32_t NumNodes() const { return mN; }
		INDEX_TYPE NodeCell(int32_t id) const { return mNodeCell[(size_t)id]; }
		float NodePers(int32_t id) const { return mMergePers[(size_t)id]; }
		int32_t NodeParent(int32_t id) const { return mParent[(size_t)id]; }

		// Diagnostic: persistence must be NONDECREASING along parent chains
		// (>= not >, equal under exact value ties). Returns violation count.
		long long CheckMonotone() const {
			long long bad = 0;
			for (int32_t i = 0; i < mN; ++i) {
				const int32_t par = mParent[(size_t)i];
				if (par == i) continue;
				if (!(mMergePers[(size_t)par] >= mMergePers[(size_t)i])) bad++;
			}
			return bad;
		}

	private:
		// "worse" = lower heap priority, so top() is the smallest persistence,
		// ties broken by smallest saddle cell - identical to the ordering in
		// UFMergeGraph::MergeArc::operator()
		struct ArcWorse {
			bool operator()(const Arc& a, const Arc& b) const {
				if (a.pers < b.pers) return false;
				if (a.pers > b.pers) return true;
				return a.saddleCell > b.saddleCell;
			}
		};

		// persistence w.r.t. the arc's stored endpoints (AddArc's formula)
		float initialPers(const Arc& a) const {
			const float d1 = std::fabs(mNodeValue[(size_t)a.e1] - a.saddleValue);
			const float d2 = std::fabs(mNodeValue[(size_t)a.e2] - a.saddleValue);
			return (d1 < d2) ? d1 : d2;
		}
		// persistence w.r.t. the current representatives (PerformCancel's
		// recompute on a stale arc)
		float currentPers(const Arc& a, int32_t r1, int32_t r2) const {
			const float d1 = std::fabs(a.saddleValue - mNodeValue[(size_t)r1]);
			const float d2 = std::fabs(a.saddleValue - mNodeValue[(size_t)r2]);
			return (d1 < d2) ? d1 : d2;
		}

		static double msSince(const std::chrono::steady_clock::time_point& t0) {
			return std::chrono::duration_cast<std::chrono::microseconds>(
				std::chrono::steady_clock::now() - t0).count() / 1000.0;
		}

		int32_t find(int32_t x) {
			while (mUf[(size_t)x] != x) {
				mUf[(size_t)x] = mUf[(size_t)mUf[(size_t)x]]; // path halving
				x = mUf[(size_t)x];
			}
			return x;
		}

		// exact replica of TopologicalMaxVertexMeshFunction::lessThan on the
		// stored per-node data: value, then Cell2HighestVertex, then cell id
		bool nodeLess(int32_t a, int32_t b) const {
			const float av = mNodeValue[(size_t)a];
			const float bv = mNodeValue[(size_t)b];
			if (av < bv) return true;
			if (bv < av) return false;
			if (mNodeTieKey[(size_t)a] < mNodeTieKey[(size_t)b]) return true;
			if (mNodeTieKey[(size_t)b] < mNodeTieKey[(size_t)a]) return false;
			return mNodeCell[(size_t)a] < mNodeCell[(size_t)b];
		}

		void applyMerge(int32_t r1, int32_t r2, const Arc& a) {
			// survivor = more extreme (min graph: lower; max graph: higher),
			// matching UFMergeGraph::MergeByVal
			int32_t rs, rd;
			if (nodeLess(r1, r2) == mMinDirection) { rs = r1; rd = r2; }
			else { rs = r2; rd = r1; }
			mUf[(size_t)rd] = rs;
			mTaint[(size_t)rs] |= mTaint[(size_t)rd];
			mParent[(size_t)rd] = rs;
			// a.pers is min(|s-r1|, |s-r2|) w.r.t. the current representatives,
			// which is exactly the distance to the dying (less extreme) side -
			// the same quantity legacy thresholds on
			mMergePers[(size_t)rd] = a.pers;
		}

		bool mMinDirection = true;
		int32_t mN = 0;
		std::vector<INDEX_TYPE> mNodeCell;
		std::vector<float> mNodeValue;
		std::vector<INDEX_TYPE> mNodeTieKey;
		std::vector<int32_t> mParent;      // survivor at death; self for roots
		std::vector<float> mMergePers;     // +inf for roots
		std::vector<int32_t> mUf;          // build-time union-find
		std::vector<uint8_t> mTaint;       // authoritative at representatives only
		std::vector<int32_t> mResolved;    // BuildRemap memo
		std::vector<uint32_t> mStamp;
		std::vector<int32_t> mPath;
		uint32_t mCurStamp = 0;
	};

} // namespace GInt

#endif
