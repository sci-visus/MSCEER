#ifndef GI_EXTREMUM_ARC_EXTRACT_2D_H
#define GI_EXTREMUM_ARC_EXTRACT_2D_H

#include <cstdint>
#include <vector>

#include "gi_basic_types.h"
#include "gi_extremum_merge_forest.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace GInt {
namespace ExtNet2D {

	// Extremum-network arc extraction for 2D, the bridge between the flow-terminal
	// base labeling and ExtremumMergeForest.
	//
	// The observation this rests on (gpu_dgrad's Label2DFlowTerminals work): the
	// saddle->extremum V-path traces the MSC builder performs are redundant with
	// the flow-terminal map. Every critical edge is a saddle joining exactly two
	// minima (the terminals of its two facet VERTICES) and two maxima (the
	// terminals of its two cofacet QUADS), so an arc is two dependent loads
	// instead of two traced paths.
	//
	// TERMINAL ENCODING. Terminal maps here are in LATTICE-INDEX space, which is
	// what gpu::Label2DFlowTerminals natively emits and what msc_2d_lib already
	// holds:
	//   ascending  termLat[vx + vy*X]      = terminal VERTEX lattice index, -1 none
	//   descending termLat[qx + qy*(X-1)]  = terminal QUAD   lattice index, -1 none
	// (The extnet2d benchmark's CPU reference walk stores terminal CELL ids
	// instead; convert with CellFromLattice2D / LatticeFromCell2D below.)
	//
	// Cell-lattice conventions, shared with the rest of GInt: a vertex (vx, vy)
	// is cell 2vx + 2vy*(2X-1); a quad (qx, qy) is cell (2qx+1) + (2qy+1)*(2X-1);
	// a cell is an edge iff exactly one of its coordinates is odd.

	// ---------------------------------------------------------------------------
	// Lattice <-> cell conversions
	// ---------------------------------------------------------------------------

	inline INDEX_TYPE CellFromLattice2D(INDEX_TYPE lat, INDEX_TYPE X, bool ascending) {
		const INDEX_TYPE cellRow = 2 * X - 1;
		if (ascending) {
			const INDEX_TYPE vx = lat % X;
			const INDEX_TYPE vy = lat / X;
			return 2 * vx + 2 * vy * cellRow;
		}
		const INDEX_TYPE QX = X - 1;
		const INDEX_TYPE qx = lat % QX;
		const INDEX_TYPE qy = lat / QX;
		return (2 * qx + 1) + (2 * qy + 1) * cellRow;
	}

	inline INDEX_TYPE LatticeFromCell2D(INDEX_TYPE cell, INDEX_TYPE X, bool ascending) {
		const INDEX_TYPE cellRow = 2 * X - 1;
		const INDEX_TYPE cx = cell % cellRow;
		const INDEX_TYPE cy = cell / cellRow;
		if (ascending) return (cx >> 1) + (cy >> 1) * X;
		return ((cx - 1) >> 1) + ((cy - 1) >> 1) * (X - 1);
	}

	// ---------------------------------------------------------------------------
	// Flow terminals (CPU)
	// ---------------------------------------------------------------------------

	// CPU equivalent of gpu::Label2DFlowTerminals, for builds without CUDA.
	// Follows the offset-doubling successor (a cell v paired with p continues at
	// 2p - v; the same identity holds for vertices flowing down and quads flowing
	// up) until it reaches a critical cell, memoizing whole chains.
	//
	// Races on the memo are benign: the flow is a deterministic forest, so two
	// threads racing on one slot can only ever write the same terminal.
	//
	// Output is in LATTICE-INDEX space (-1 = never assigned), matching the GPU.
	template <typename GradType>
	void CpuFlowTerminals2D(GradType* grad, INDEX_TYPE X, INDEX_TYPE Y, bool ascending,
		std::vector<int32_t>& termLat) {
		const INDEX_TYPE cellRow = 2 * X - 1;
		const INDEX_TYPE NX = ascending ? X : (X - 1);
		const INDEX_TYPE NY = ascending ? Y : (Y - 1);
		const INDEX_TYPE n = NX * NY;
		termLat.assign((size_t)(n > 0 ? n : 0), -1);
		if (n <= 0) return;

#if defined(_OPENMP)
#pragma omp parallel
#endif
		{
			std::vector<INDEX_TYPE> chain;
#if defined(_OPENMP)
#pragma omp for schedule(dynamic, 4096)
#endif
			for (INDEX_TYPE i = 0; i < n; ++i) {
				if (termLat[(size_t)i] >= 0) continue;
				const INDEX_TYPE x = i % NX;
				const INDEX_TYPE y = i / NX;
				INDEX_TYPE cur = ascending ? (2 * x + 2 * y * cellRow)
					: ((2 * x + 1) + (2 * y + 1) * cellRow);
				chain.clear();
				int32_t terminal = -1;
				while (true) {
					const INDEX_TYPE lat = LatticeFromCell2D(cur, X, ascending);
					const int32_t memo = termLat[(size_t)lat];
					if (memo >= 0) { terminal = memo; break; }
					if (grad->getCritical(cur)) { terminal = (int32_t)lat; break; }
					chain.push_back(lat);
					const INDEX_TYPE pair = (INDEX_TYPE)grad->getPair(cur);
					cur = 2 * pair - cur;
				}
				for (size_t c = 0; c < chain.size(); ++c) termLat[(size_t)chain[c]] = terminal;
				// the start slot is not in the chain when the walk ended
				// immediately (critical start or memo hit), so stamp it too
				termLat[(size_t)i] = terminal;
			}
		}
	}

	// ---------------------------------------------------------------------------
	// Base regions
	// ---------------------------------------------------------------------------

	// The compact base regions of one direction: one region per distinct flow
	// terminal. Region id IS the ExtremumMergeForest node id, which is what lets
	// BuildRemap's output be used directly as a paintLabels() remap with no
	// indirection anywhere.
	struct BaseRegions2D {
		std::vector<INDEX_TYPE> terminalCells;  // region -> extremum cell id
		std::vector<int32_t> latToNode;         // lattice index -> region id, else -1
		std::vector<long long> regionWeight;    // region -> #lattice points
		long long totalWeight = 0;

		int32_t count() const { return (int32_t)terminalCells.size(); }
	};

	// One pass over the terminal map. Regions are numbered in FIRST-SEEN lattice
	// order, which makes the numbering a deterministic function of the gradient
	// alone -- it must not depend on thread scheduling, because the id space is
	// the one every downstream label carries.
	inline void BuildBaseRegions2D(const int32_t* termLat, INDEX_TYPE X, INDEX_TYPE Y,
		bool ascending, BaseRegions2D& out) {
		const INDEX_TYPE NX = ascending ? X : (X - 1);
		const INDEX_TYPE NY = ascending ? Y : (Y - 1);
		const size_t n = (size_t)(NX * NY);

		out.terminalCells.clear();
		out.regionWeight.clear();
		out.totalWeight = 0;
		out.latToNode.assign(n, -1);
		if (NX <= 0 || NY <= 0) return;

		for (size_t i = 0; i < n; ++i) {
			const int32_t t = termLat[i];
			if (t < 0) continue;
			int32_t r = out.latToNode[(size_t)t];
			if (r < 0) {
				r = (int32_t)out.terminalCells.size();
				out.latToNode[(size_t)t] = r;
				out.terminalCells.push_back(CellFromLattice2D((INDEX_TYPE)t, X, ascending));
				out.regionWeight.push_back(0);
			}
			out.regionWeight[(size_t)r]++;
			out.totalWeight++;
		}
	}

	// Per-pixel compact region id, for a caller that wants the base labeling
	// itself (msc_2d_lib's ensureBaseLabeling). `remap`, when non-null, is applied
	// to the region id first (used to collapse the base partition at
	// basePersistenceAbs); it must be a permutation-or-merge into [0, remapCount).
	inline void PaintBaseRegions2D(const BaseRegions2D& regions, const int32_t* termLat,
		INDEX_TYPE X, INDEX_TYPE Y, bool ascending,
		const int32_t* remap, int32_t* outLabels) {
		const INDEX_TYPE NX = ascending ? X : (X - 1);
		const INDEX_TYPE NY = ascending ? Y : (Y - 1);
		const INDEX_TYPE n = NX * NY;
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
		for (INDEX_TYPE i = 0; i < n; ++i) {
			const int32_t t = termLat[(size_t)i];
			int32_t r = (t < 0) ? -1 : regions.latToNode[(size_t)t];
			if (r >= 0 && remap != NULL) r = remap[(size_t)r];
			outLabels[(size_t)i] = r;
		}
	}

	// ---------------------------------------------------------------------------
	// Node payloads
	// ---------------------------------------------------------------------------

	// value + tie key per node, in the exact total order
	// TopologicalMaxVertexMeshFunction::lessThan uses, so the forest's survivor
	// choice reproduces UFMergeGraph::MergeByVal (and therefore the MSC's).
	template <typename MeshFuncType, typename MaxVLabelingType>
	void BuildNodePayloads2D(const BaseRegions2D& regions, MeshFuncType* meshfunc,
		MaxVLabelingType* maxv, std::vector<float>& values,
		std::vector<INDEX_TYPE>& tieKeys) {
		const INDEX_TYPE n = (INDEX_TYPE)regions.terminalCells.size();
		values.resize((size_t)n);
		tieKeys.resize((size_t)n);
#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
		for (INDEX_TYPE i = 0; i < n; ++i) {
			const INDEX_TYPE cell = regions.terminalCells[(size_t)i];
			values[(size_t)i] = (float)meshfunc->cellValue(cell);
			tieKeys[(size_t)i] = (INDEX_TYPE)maxv->Cell2HighestVertex(cell);
		}
	}

	// ---------------------------------------------------------------------------
	// Arc extraction
	// ---------------------------------------------------------------------------

	// Arc lists shaped for ExtremumMergeForest::Build: per-band local lists (both
	// endpoints owned by the band), one cross list, and the frozen flags. The
	// forest's contract requires BOTH endpoints of every cross arc to be flagged
	// -- that is a correctness requirement, not an optimization -- so the frozen
	// stores happen at emission time next to the classification.
	struct ArcBuffers2D {
		std::vector<std::vector<ExtremumMergeForest::Arc> > local;
		std::vector<ExtremumMergeForest::Arc> cross;
		std::vector<uint8_t> frozen;
	};

	// Bands are row bands of the doubled cell lattice, applied identically to
	// saddles and to extremum cells so that "same band" is well defined across
	// cell dimensions.
	inline int BandOfCellRow2D(INDEX_TYPE cy, INDEX_TYPE cellRows, int numBands) {
		return (int)((cy * (INDEX_TYPE)numBands) / cellRows);
	}

	// More bands cut the parallel local phase but push more arcs into the serial
	// residue heap (measured residue share 15% -> 25% -> 37% at 4 -> 8 -> 16
	// bands). 8 is the measured optimum at both 2048^2 and 3200^2; 4 and 16 are
	// ~15-40% worse. `requested` <= 0 means auto.
	inline int DefaultBandCount2D(INDEX_TYPE Y, int requested) {
		int bands = requested;
		if (bands <= 0) {
#if defined(_OPENMP)
			bands = omp_get_max_threads();
			if (bands > 8) bands = 8;
#else
			bands = 1;
#endif
		}
		if (bands < 1) bands = 1;
		const INDEX_TYPE cellRows = 2 * Y - 1;
		if ((INDEX_TYPE)bands > cellRows) bands = (int)cellRows;
		return bands;
	}

	// One banded parallel pass over the doubled lattice emitting BOTH directions'
	// arcs. Sharing the pass matters: the critical-edge test and the cell
	// arithmetic are the expensive part, and each critical edge contributes to
	// the min graph and (unless it is on the boundary) the max graph.
	//
	// Either output may be skipped by passing a null terminal map for that
	// direction.
	template <typename GradType, typename MeshFuncType>
	void ExtractArcsBanded2D(GradType* grad, MeshFuncType* meshfunc,
		INDEX_TYPE X, INDEX_TYPE Y, int numBands,
		const int32_t* ascTermLat, const int32_t* dscTermLat,
		const BaseRegions2D* ascRegions, const BaseRegions2D* dscRegions,
		ArcBuffers2D* minOut, ArcBuffers2D* maxOut) {
		const INDEX_TYPE cellRow = 2 * X - 1;
		const INDEX_TYPE cellRows = 2 * Y - 1;
		const INDEX_TYPE QX = X - 1;
		const bool doMin = (minOut != NULL && ascRegions != NULL && ascTermLat != NULL);
		const bool doMax = (maxOut != NULL && dscRegions != NULL && dscTermLat != NULL);

		if (doMin) {
			minOut->local.assign((size_t)numBands, std::vector<ExtremumMergeForest::Arc>());
			minOut->frozen.assign((size_t)ascRegions->count(), 0);
			minOut->cross.clear();
		}
		if (doMax) {
			maxOut->local.assign((size_t)numBands, std::vector<ExtremumMergeForest::Arc>());
			maxOut->frozen.assign((size_t)dscRegions->count(), 0);
			maxOut->cross.clear();
		}
		if (!doMin && !doMax) return;

		std::vector<std::vector<ExtremumMergeForest::Arc> > minCrossPer((size_t)numBands);
		std::vector<std::vector<ExtremumMergeForest::Arc> > maxCrossPer((size_t)numBands);

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 1)
#endif
		for (int b = 0; b < numBands; ++b) {
			const INDEX_TYPE y0 = ((INDEX_TYPE)b * cellRows) / numBands;
			const INDEX_TYPE y1 = ((INDEX_TYPE)(b + 1) * cellRows) / numBands;
			std::vector<ExtremumMergeForest::Arc>* minL = doMin ? &minOut->local[(size_t)b] : NULL;
			std::vector<ExtremumMergeForest::Arc>* maxL = doMax ? &maxOut->local[(size_t)b] : NULL;
			std::vector<ExtremumMergeForest::Arc>* minC = doMin ? &minCrossPer[(size_t)b] : NULL;
			std::vector<ExtremumMergeForest::Arc>* maxC = doMax ? &maxCrossPer[(size_t)b] : NULL;

			for (INDEX_TYPE cy = y0; cy < y1; ++cy) {
				for (INDEX_TYPE cx = 0; cx < cellRow; ++cx) {
					if (((cx ^ cy) & 1) == 0) continue;          // not an edge
					const INDEX_TYPE s = cx + cy * cellRow;
					if (!grad->getCritical(s)) continue;
					const bool xEdge = ((cx & 1) == 1);
					const float sv = (float)meshfunc->cellValue(s);

					// --- min side: terminals of the two facet vertices ---
					if (doMin) {
						const INDEX_TYPE v1 = xEdge ? (s - 1) : (s - cellRow);
						const INDEX_TYPE v2 = xEdge ? (s + 1) : (s + cellRow);
						const int32_t t1 = ascTermLat[(size_t)LatticeFromCell2D(v1, X, true)];
						const int32_t t2 = ascTermLat[(size_t)LatticeFromCell2D(v2, X, true)];
						const int32_t e1 = (t1 < 0) ? -1 : ascRegions->latToNode[(size_t)t1];
						const int32_t e2 = (t2 < 0) ? -1 : ascRegions->latToNode[(size_t)t2];
						if (e1 >= 0 && e2 >= 0 && e1 != e2) {
							ExtremumMergeForest::Arc a;
							a.e1 = e1; a.e2 = e2; a.saddleValue = sv; a.saddleCell = s; a.pers = 0.0f;
							// LOCAL iff BOTH extrema live in the saddle's own band
							const int b1 = BandOfCellRow2D(
								ascRegions->terminalCells[(size_t)e1] / cellRow, cellRows, numBands);
							const int b2 = BandOfCellRow2D(
								ascRegions->terminalCells[(size_t)e2] / cellRow, cellRows, numBands);
							if (b1 == b && b2 == b) {
								minL->push_back(a);
							} else {
								minC->push_back(a);
								// idempotent 1-stores; both endpoints, per the
								// forest's frozen-closure requirement
								minOut->frozen[(size_t)e1] = 1;
								minOut->frozen[(size_t)e2] = 1;
							}
						}
					}

					// --- max side: terminals of the two cofacet quads ---
					// Boundary edges have only one cofacet quad and are skipped
					// exactly as trace_up_2saddle does.
					if (!doMax) continue;
					const bool interior = xEdge ? (cy > 0 && cy < cellRows - 1)
						: (cx > 0 && cx < cellRow - 1);
					if (!interior) continue;
					{
						const INDEX_TYPE q1 = xEdge ? (s - cellRow) : (s - 1);
						const INDEX_TYPE q2 = xEdge ? (s + cellRow) : (s + 1);
						const int32_t t1 = dscTermLat[(size_t)LatticeFromCell2D(q1, X, false)];
						const int32_t t2 = dscTermLat[(size_t)LatticeFromCell2D(q2, X, false)];
						const int32_t e1 = (t1 < 0) ? -1 : dscRegions->latToNode[(size_t)t1];
						const int32_t e2 = (t2 < 0) ? -1 : dscRegions->latToNode[(size_t)t2];
						if (e1 >= 0 && e2 >= 0 && e1 != e2) {
							ExtremumMergeForest::Arc a;
							a.e1 = e1; a.e2 = e2; a.saddleValue = sv; a.saddleCell = s; a.pers = 0.0f;
							const int b1 = BandOfCellRow2D(
								dscRegions->terminalCells[(size_t)e1] / cellRow, cellRows, numBands);
							const int b2 = BandOfCellRow2D(
								dscRegions->terminalCells[(size_t)e2] / cellRow, cellRows, numBands);
							if (b1 == b && b2 == b) {
								maxL->push_back(a);
							} else {
								maxC->push_back(a);
								maxOut->frozen[(size_t)e1] = 1;
								maxOut->frozen[(size_t)e2] = 1;
							}
						}
					}
				}
			}
		}

		// Concatenate in band order so the cross list is a deterministic function
		// of the gradient (the forest sorts it, but determinism of the input keeps
		// equal-key ordering reproducible).
		for (int b = 0; b < numBands; ++b) {
			if (doMin) minOut->cross.insert(minOut->cross.end(),
				minCrossPer[(size_t)b].begin(), minCrossPer[(size_t)b].end());
			if (doMax) maxOut->cross.insert(maxOut->cross.end(),
				maxCrossPer[(size_t)b].begin(), maxCrossPer[(size_t)b].end());
		}
	}

	// Diagnostic mirroring the benchmark's extC_frozen_closure gate: every cross
	// arc must have BOTH endpoints frozen, or the partitioned build is not
	// order-consistent with the serial one.
	inline bool CheckFrozenClosure2D(const ArcBuffers2D& bufs) {
		for (size_t i = 0; i < bufs.cross.size(); ++i) {
			const ExtremumMergeForest::Arc& a = bufs.cross[i];
			if (a.e1 < 0 || a.e2 < 0) return false;
			if (!bufs.frozen[(size_t)a.e1] || !bufs.frozen[(size_t)a.e2]) return false;
		}
		return true;
	}

} // namespace ExtNet2D
} // namespace GInt

#endif
