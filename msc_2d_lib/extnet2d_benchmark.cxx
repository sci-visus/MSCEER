// extnet2d_benchmark: time the extremum-network destructive simplification
// (SimplifiedExtremumGraph / UFMergeGraph, built directly from the discrete
// gradient) against the full MSC -> ComputeHierarchy -> SetSelectPersAbs path
// (serial and partitioned) at a sweep of persistence thresholds, and measure
// how closely the two segmentations agree when painted over the same base
// (flow-terminal) regions.
//
// Usage: extnet2d_benchmark [size=2048] [partitions=4] [use_gpu=1] [nthresh=10] [runB_repeat=1]
//   thresholds are k% of the function value range for k = 1..nthresh.

#include "gi_discrete_gradient_computer.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_partitioned_topological_regular_grid.h"
#include "gi_extrema_region_builder.h"

#ifdef EXTNET2D_HAS_GPU
#include "dgrad_gpu_api.h"
#include "dgrad2d_tables.h"
#endif

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

typedef GInt::Accurate2D::MeshType MeshType;
typedef GInt::Accurate2D::MeshFuncType MeshFuncType;
typedef GInt::Accurate2D::GradType GradType;
typedef GInt::MorseSmaleComplexBasic<float, MeshType, MeshFuncType, GradType> SerialMscType;
typedef GInt::MorseSmaleComplexPartitioned<float, MeshType, MeshFuncType, GradType> PartitionedPipelineType;
typedef GInt::SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType> ExtGraphType;
typedef GInt::UFMergeGraph<MeshFuncType, MeshType> UFGraphType;

static double msSince(const std::chrono::steady_clock::time_point& t0) {
	return std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::steady_clock::now() - t0).count() / 1000.0;
}

static void smoothFieldInPlace(std::vector<float>& field, int rows, int cols, int iterations) {
	if (iterations <= 0) return;
	std::vector<float> scratch(field.size(), 0.0f);
	for (int it = 0; it < iterations; ++it) {
		for (int r = 0; r < rows; ++r) {
			const int r0 = (r > 0) ? r - 1 : r;
			const int r1 = (r < rows - 1) ? r + 1 : r;
			for (int c = 0; c < cols; ++c) {
				const int c0 = (c > 0) ? c - 1 : c;
				const int c1 = (c < cols - 1) ? c + 1 : c;
				float sum = 0.0f;
				int n = 0;
				for (int rr = r0; rr <= r1; ++rr) {
					for (int cc = c0; cc <= c1; ++cc) {
						sum += field[(size_t)rr * (size_t)cols + (size_t)cc];
						n++;
					}
				}
				scratch[(size_t)r * (size_t)cols + (size_t)c] = sum / static_cast<float>(n);
			}
		}
		field.swap(scratch);
	}
}

// -----------------------------------------------------------------------------
// Base (flow-terminal) labeling. Every vertex flows to its terminal critical
// vertex (basin-of-minimum identity), every quad to its terminal critical quad
// (region-of-maximum identity). Terminals are stored as CELL ids on the doubled
// lattice so they key directly into the MSC nodemap and the extremum graph.
// CPU walk uses the offset-doubling successor (v paired with edge e continues
// at 2e - v; same identity for quads) with memoization; races on the memo are
// benign because the flow is a deterministic forest.
// -----------------------------------------------------------------------------

static void cpuFlowTerminals(GradType* grad, int64_t X, int64_t Y, bool ascending,
	std::vector<int64_t>& termCell) {
	const int64_t cellRow = 2 * X - 1;
	const int64_t NX = ascending ? X : (X - 1);
	const int64_t NY = ascending ? Y : (Y - 1);
	const int64_t n = NX * NY;
	termCell.assign((size_t)n, -1);
	std::vector<int64_t> chain;
#pragma omp parallel private(chain)
	{
		chain.clear();
#pragma omp for schedule(dynamic, 4096)
		for (int64_t i = 0; i < n; ++i) {
			if (termCell[(size_t)i] >= 0) continue;
			const int64_t x = i % NX;
			const int64_t y = i / NX;
			int64_t cur = ascending ? (2 * x + 2 * y * cellRow)
				: ((2 * x + 1) + (2 * y + 1) * cellRow);
			chain.clear();
			int64_t terminal = -1;
			while (true) {
				const int64_t cx = cur % cellRow;
				const int64_t cy = cur / cellRow;
				const int64_t lat = ascending ? ((cx >> 1) + (cy >> 1) * NX)
					: (((cx - 1) >> 1) + ((cy - 1) >> 1) * NX);
				const int64_t memo = termCell[(size_t)lat];
				if (memo >= 0) { terminal = memo; break; }
				if (grad->getCritical(cur)) { terminal = cur; break; }
				chain.push_back(lat);
				const int64_t pair = grad->getPair(cur);
				cur = 2 * pair - cur;
			}
			for (size_t c = 0; c < chain.size(); ++c) termCell[(size_t)chain[c]] = terminal;
			// the start slot is not in the chain when the walk terminated
			// immediately (critical start or memo hit), so stamp it explicitly
			termCell[(size_t)i] = terminal;
		}
	}
}

struct BaseRegions {
	std::vector<int64_t> terminalCells;              // region -> terminal cell id
	std::vector<long long> regionWeight;             // region -> #pixels (or #quads)
	std::unordered_map<int64_t, int> cellToRegion;   // terminal cell id -> region
	long long totalWeight = 0;
};

static void buildRegions(const std::vector<int64_t>& termCell, BaseRegions& out) {
	out.terminalCells.clear();
	out.regionWeight.clear();
	out.cellToRegion.clear();
	out.totalWeight = 0;
	for (size_t i = 0; i < termCell.size(); ++i) {
		const int64_t t = termCell[i];
		if (t < 0) continue;
		std::unordered_map<int64_t, int>::iterator it = out.cellToRegion.find(t);
		int r;
		if (it == out.cellToRegion.end()) {
			r = (int)out.terminalCells.size();
			out.cellToRegion[t] = r;
			out.terminalCells.push_back(t);
			out.regionWeight.push_back(0);
		}
		else {
			r = it->second;
		}
		out.regionWeight[(size_t)r]++;
		out.totalWeight++;
	}
}

// -----------------------------------------------------------------------------
// Path A per-threshold labeling: parallel walk-down of the merged-manifold
// hierarchy from every living extremum (same traversal as msc_2d_lib's
// buildRemapWalkDown), stamping the living node id onto each base region.
// -----------------------------------------------------------------------------

struct PathAContext {
	std::vector<int> baseNodeOfRegion;  // region -> base MSC node id
	std::vector<int> regionOfBaseNode;  // base MSC node id -> region (-1 = none)
	bool ok = false;
};

static bool buildPathAContext(SerialMscType* msc, const BaseRegions& regions, int wantDim,
	PathAContext& ctx) {
	std::unordered_map<int64_t, int> cellToNode;
	for (INT_TYPE nid = 0; nid < msc->numNodes(); ++nid) {
		if (msc->getNode(nid).dim == wantDim) {
			cellToNode[(int64_t)msc->getNode(nid).cellindex] = (int)nid;
		}
	}
	ctx.baseNodeOfRegion.assign(regions.terminalCells.size(), -1);
	ctx.regionOfBaseNode.assign((size_t)msc->numNodes(), -1);
	for (size_t r = 0; r < regions.terminalCells.size(); ++r) {
		std::unordered_map<int64_t, int>::const_iterator it =
			cellToNode.find(regions.terminalCells[r]);
		if (it == cellToNode.end()) return false;
		ctx.baseNodeOfRegion[r] = it->second;
		ctx.regionOfBaseNode[(size_t)it->second] = (int)r;
	}
	ctx.ok = true;
	return true;
}

// returns number of living roots; fills labelA (region -> living node id)
static size_t pathALabels(SerialMscType* msc, const PathAContext& ctx, bool ascending,
	int wantDim, std::vector<int>& labelA, long long& conflicts) {
	labelA.assign(ctx.baseNodeOfRegion.size(), -1);
	std::vector<int> roots;
	SerialMscType::LivingNodesIterator nit(msc);
	for (nit.begin(); nit.valid(); nit.advance()) {
		const INT_TYPE nid = nit.value();
		if (msc->getNode(nid).dim == wantDim) roots.push_back((int)nid);
	}
	const long long nroots = (long long)roots.size();
	int* dense = labelA.empty() ? NULL : &labelA[0];
	long long conf = 0;
#pragma omp parallel
	{
		std::vector<INT_TYPE> stack;
#pragma omp for schedule(dynamic, 64) reduction(+ : conf)
		for (long long r = 0; r < nroots; ++r) {
			const int nid = roots[(size_t)r];
			stack.clear();
			stack.push_back(msc->activeManifoldForNode(nid, ascending));
			while (!stack.empty()) {
				const INT_TYPE man = stack.back();
				stack.pop_back();
				const GInt::merged_manifold& mm = msc->manifoldAt(man);
				if (mm.merged[0] != -1) {
					stack.push_back(mm.merged[0]);
					stack.push_back(mm.merged[1]);
				}
				const int reg = ctx.regionOfBaseNode[(size_t)mm.basenode];
				if (reg >= 0) {
					if (dense[reg] >= 0 && dense[reg] != nid) conf++;
					dense[reg] = nid;
				}
			}
		}
	}
	conflicts = conf;
	return roots.size();
}

// -----------------------------------------------------------------------------
// Agreement metrics between two region labelings, weighted by region size.
// -----------------------------------------------------------------------------

struct Agreement {
	long long countA = 0, countB = 0;
	double purityAB = 0.0;   // sum over B clusters of max-overlap with an A cluster
	double purityBA = 0.0;   // sum over A clusters of max-overlap with a B cluster
	double mutual = 0.0;     // pixels in (a,b) cells that are each other's argmax
	long long weight = 0;
};

static Agreement computeAgreement(const std::vector<int>& labelA,
	const std::vector<int64_t>& labelB,
	const std::vector<long long>& weight) {
	Agreement ag;
	std::map<std::pair<int, int64_t>, long long> contingency;
	std::map<int, long long> weightA;
	std::map<int64_t, long long> weightB;
	for (size_t r = 0; r < labelA.size(); ++r) {
		if (labelA[r] < 0 || labelB[r] < 0) continue;
		const long long w = weight[r];
		contingency[std::make_pair(labelA[r], labelB[r])] += w;
		weightA[labelA[r]] += w;
		weightB[labelB[r]] += w;
		ag.weight += w;
	}
	ag.countA = (long long)weightA.size();
	ag.countB = (long long)weightB.size();
	if (ag.weight == 0) return ag;

	std::map<int, std::pair<long long, int64_t> > bestB;   // A -> (w, B)
	std::map<int64_t, std::pair<long long, int> > bestA;   // B -> (w, A)
	long long sumMaxPerB = 0;
	for (std::map<std::pair<int, int64_t>, long long>::const_iterator it = contingency.begin();
		it != contingency.end(); ++it) {
		const int a = it->first.first;
		const int64_t b = it->first.second;
		const long long w = it->second;
		if (w > bestB[a].first) bestB[a] = std::make_pair(w, b);
		if (w > bestA[b].first) bestA[b] = std::make_pair(w, a);
	}
	long long sumMaxPerA = 0;
	for (std::map<int, std::pair<long long, int64_t> >::const_iterator it = bestB.begin();
		it != bestB.end(); ++it) sumMaxPerA += it->second.first;
	for (std::map<int64_t, std::pair<long long, int> >::const_iterator it = bestA.begin();
		it != bestA.end(); ++it) sumMaxPerB += it->second.first;

	long long mutualW = 0;
	for (std::map<int, std::pair<long long, int64_t> >::const_iterator it = bestB.begin();
		it != bestB.end(); ++it) {
		const int a = it->first;
		const int64_t b = it->second.second;
		if (bestA[b].second == a) {
			mutualW += contingency[std::make_pair(a, b)];
		}
	}
	ag.purityBA = (double)sumMaxPerA / (double)ag.weight;
	ag.purityAB = (double)sumMaxPerB / (double)ag.weight;
	ag.mutual = (double)mutualW / (double)ag.weight;
	return ag;
}

int main(int argc, char** argv) {
	int size = 2048;
	int partitions = 4;
	int use_gpu = 1;
	int nthresh = 10;
	int runB_repeat = 1;
	if (argc > 1) size = std::atoi(argv[1]);
	if (argc > 2) partitions = std::atoi(argv[2]);
	if (argc > 3) use_gpu = std::atoi(argv[3]);
	if (argc > 4) nthresh = std::atoi(argv[4]);
	if (argc > 5) runB_repeat = std::atoi(argv[5]);
	if (nthresh < 1) nthresh = 1;
	if (runB_repeat < 1) runB_repeat = 1;
	if (!GInt::PartitionedTopologicalRegularGrid2D::IsSupportedPartitionCount(partitions)) {
		fprintf(stderr, "Unsupported partition count %d (allowed: 1,2,3,4,6,8,9,12,16)\n", partitions);
		return 2;
	}
	const int rows = size;
	const int cols = size;

	// --- synthetic multiscale field (identical generator to the partitioned smoke)
	std::vector<float> field((size_t)rows * (size_t)cols, 0.0f);
	{
		std::mt19937 rng(123456789u);
		std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
		std::vector<float> coarseNoise(field.size(), 0.0f);
		std::vector<float> detailNoise(field.size(), 0.0f);
		for (size_t i = 0; i < field.size(); ++i) coarseNoise[i] = dist(rng);
		for (size_t i = 0; i < field.size(); ++i) detailNoise[i] = dist(rng);
		smoothFieldInPlace(coarseNoise, rows, cols, 30);
		smoothFieldInPlace(detailNoise, rows, cols, 2);
		const float coarseScale = 100.0f;
		const float detailScale = 0.05f * coarseScale;
		for (size_t i = 0; i < field.size(); ++i)
			field[i] = coarseScale * coarseNoise[i] + detailScale * detailNoise[i];
	}

	// --- Phase 0: discrete gradient (shared by every path)
	bool gpuActive = false;
#ifdef EXTNET2D_HAS_GPU
	if (use_gpu) {
		GInt::gpu::DeviceInfo info;
		if (GInt::gpu::QueryDevice(info)) {
			gpuActive = true;
			printf("extnet2d_device name=%s sm=%d\n", info.name, info.sm_count);
		}
		else {
			printf("extnet2d_device none; falling back to CPU\n");
		}
	}
#else
	if (use_gpu) printf("extnet2d_device built_without_gpu; using CPU\n");
#endif

	GInt::Accurate2D::DiscreteGradientBuilder dgb;
	dgb.SetFloadArrayAndDims(cols, rows, field.data());
	dgb.SetNeededAccuracy(true, true);
	dgb.SetParallelism(partitions);
#ifdef EXTNET2D_HAS_GPU
	if (gpuActive) {
		dgb.SetGradientOverride(
			[](GInt::Accurate2D::GridType* grid, GInt::Accurate2D::GridFuncType* func,
				GInt::Accurate2D::MeshType* mesh, GInt::Accurate2D::GradType* grad) -> bool {
				GInt::gpu::Dgrad2DTables tables;
				if (!GInt::gpu::BuildDgrad2DTablesFromMesh(*mesh, tables)) return false;
				const auto xy = grid->XY();
				return GInt::gpu::ComputeDiscreteGradient2DFused(
					func->GetImage(), xy[0], xy[1], tables,
					reinterpret_cast<uint8_t*>(grad->m_dgrad->LabelArray()),
					NULL);
			});
	}
#endif
	const auto tGrad = std::chrono::steady_clock::now();
	dgb.ComputeDiscreteGradient();
	const double gradMs = msSince(tGrad);

	MeshType* mesh = dgb.GetTopoMesh();
	MeshFuncType* meshfunc = dgb.GetMeshFunc();
	GradType* grad = dgb.GetGrad();
	GInt::Accurate2D::GridFuncType* gridfunc = dgb.GetGridFunc();

	const float range = gridfunc->GetMaxValue() - gridfunc->GetMinValue();
	std::vector<float> thresholds;
	for (int k = 1; k <= nthresh; ++k) thresholds.push_back(0.01f * (float)k * range);
	const float tMax = thresholds.back();

	printf("extnet2d_config size=%d partitions=%d gpu=%d nthresh=%d repeat=%d range=%f tmax=%f\n",
		size, partitions, gpuActive ? 1 : 0, nthresh, runB_repeat, range, tMax);
	printf("extnet2d_ms phase=gradient ms=%.3f\n", gradMs);

	// --- Phase 1: base labeling (shared, threshold-independent)
	const int64_t X = cols, Y = rows;
	std::vector<int64_t> ascTerm, dscTerm;
	bool gpuLabelUsed = false;
#ifdef EXTNET2D_HAS_GPU
	if (gpuActive) {
		GInt::gpu::Dgrad2DTables tables;
		if (GInt::gpu::BuildDgrad2DTablesFromMesh(*mesh, tables)) {
			const uint8_t* gradBytes =
				reinterpret_cast<const uint8_t*>(grad->m_dgrad->LabelArray());
			std::vector<int32_t> vterm((size_t)(X * Y));
			std::vector<int32_t> qterm((size_t)((X - 1) * (Y - 1)));
			const auto tLab = std::chrono::steady_clock::now();
			if (GInt::gpu::Label2DFlowTerminals(gradBytes, X, Y, tables,
				vterm.data(), qterm.data(), NULL)) {
				const double labMs = msSince(tLab);
				const int64_t cellRow = 2 * X - 1;
				const int64_t QX = X - 1;
				ascTerm.resize(vterm.size());
				dscTerm.resize(qterm.size());
				for (size_t i = 0; i < vterm.size(); ++i) {
					const int64_t t = vterm[i];
					ascTerm[i] = (t < 0) ? -1 : (2 * (t % X) + 2 * (t / X) * cellRow);
				}
				for (size_t i = 0; i < qterm.size(); ++i) {
					const int64_t t = qterm[i];
					dscTerm[i] = (t < 0) ? -1 : ((2 * (t % QX) + 1) + (2 * (t / QX) + 1) * cellRow);
				}
				gpuLabelUsed = true;
				printf("extnet2d_ms phase=base_labeling src=gpu ms=%.3f\n", labMs);
			}
		}
	}
#endif
	{
		std::vector<int64_t> ascCpu, dscCpu;
		const auto tLab = std::chrono::steady_clock::now();
		cpuFlowTerminals(grad, X, Y, true, ascCpu);
		cpuFlowTerminals(grad, X, Y, false, dscCpu);
		const double labMs = msSince(tLab);
		printf("extnet2d_ms phase=base_labeling src=cpu ms=%.3f\n", labMs);
		if (gpuLabelUsed) {
			const bool same = (ascCpu == ascTerm) && (dscCpu == dscTerm);
			printf("extnet2d_gate name=terminal_gpu_cpu_match pass=%d\n", same ? 1 : 0);
			if (!same) return 3;
		}
		else {
			ascTerm.swap(ascCpu);
			dscTerm.swap(dscCpu);
		}
	}

	BaseRegions ascRegions, dscRegions;
	buildRegions(ascTerm, ascRegions);
	buildRegions(dscTerm, dscRegions);
	printf("extnet2d_regions asc=%lld dsc=%lld\n",
		(long long)ascRegions.terminalCells.size(), (long long)dscRegions.terminalCells.size());

	// --- Phase 2: Path A, serial MSC
	const auto tConstruct = std::chrono::steady_clock::now();
	SerialMscType serialMSC(grad, mesh, meshfunc);
	serialMSC.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
	serialMSC.ComputeFromGrad();
	const double constructMs = msSince(tConstruct);
	const auto tHier = std::chrono::steady_clock::now();
	serialMSC.ComputeHierarchy(tMax);
	const double hierMs = msSince(tHier);
	printf("extnet2d_ms phase=mscA_serial_construct ms=%.3f\n", constructMs);
	printf("extnet2d_ms phase=mscA_serial_hierarchy ms=%.3f\n", hierMs);

	PathAContext ctxAsc, ctxDsc;
	long long baseMins = 0, baseMaxs = 0;
	for (INT_TYPE nid = 0; nid < serialMSC.numNodes(); ++nid) {
		const int d = serialMSC.getNode(nid).dim;
		if (d == 0) baseMins++;
		else if (d == 2) baseMaxs++;
	}
	const bool okAsc = buildPathAContext(&serialMSC, ascRegions, 0, ctxAsc);
	const bool okDsc = buildPathAContext(&serialMSC, dscRegions, 2, ctxDsc);
	printf("extnet2d_gate name=terminals_are_base_nodes pass=%d\n", (okAsc && okDsc) ? 1 : 0);
	if (!okAsc || !okDsc) return 4;

	// --- Phase 2b: Path A, partitioned pipeline (timing + living-count cross-check)
	PartitionedPipelineType::TimingBreakdown partTimings;
	PartitionedPipelineType partitioned(grad, mesh, meshfunc);
	const float localPers = thresholds[0];
	std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
		partitioned.BuildPartitionLocalMSCs(partitions, localPers, &partTimings);
	GInt::PartitionedTopologicalRegularGrid2D partitionMesh(mesh, partitions);
	std::unique_ptr<PartitionedPipelineType::ReconciledGlobalMsc> reconciled =
		partitioned.BuildReconciledGlobalBase(partitionMesh, localResults, &partTimings);
	const auto tPartHier = std::chrono::steady_clock::now();
	reconciled->ComputeHierarchy(tMax);
	const double partHierMs = msSince(tPartHier);
	printf("extnet2d_ms phase=mscA_part localPers=%f local_total=%lld local_build=%lld local_simplify=%lld reconcile=%lld global_hier=%.3f\n",
		localPers,
		(long long)partTimings.local_stage_total_ms,
		(long long)partTimings.local_build_ms,
		(long long)partTimings.local_simplify_ms,
		(long long)partTimings.reconcile_ms,
		partHierMs);

	// --- Phase 3: Path B, extremum network
	ExtGraphType extg(mesh, meshfunc, grad);
	extg.SetMode(ExtGraphType::EXTGRAPHMODE::BOTH);
	extg.mMinGraph->SetRetainArcs(true);
	extg.mMaxGraph->SetRetainArcs(true);
	const auto tBuild = std::chrono::steady_clock::now();
	extg.BuildGraphFromGradient();
	const double extBuildMs = msSince(tBuild);
	printf("extnet2d_ms phase=extB_build ms=%.3f\n", extBuildMs);

	const bool countsMatch =
		(baseMins == (long long)extg.mMinGraph->NumExtrema()) &&
		(baseMaxs == (long long)extg.mMaxGraph->NumExtrema()) &&
		(baseMins == (long long)ascRegions.terminalCells.size()) &&
		(baseMaxs == (long long)dscRegions.terminalCells.size());
	printf("extnet2d_gate name=base_counts mscMins=%lld ufMins=%d ascRegions=%lld mscMaxs=%lld ufMaxs=%d dscRegions=%lld pass=%d\n",
		baseMins, extg.mMinGraph->NumExtrema(), (long long)ascRegions.terminalCells.size(),
		baseMaxs, extg.mMaxGraph->NumExtrema(), (long long)dscRegions.terminalCells.size(),
		countsMatch ? 1 : 0);
	if (!countsMatch) return 5;

	// one-time honest datapoint: full rebuild (scan + traces + simplify) with no
	// retained-arc refactor, at the middle threshold
	{
		ExtGraphType extg2(mesh, meshfunc, grad);
		extg2.SetMode(ExtGraphType::EXTGRAPHMODE::BOTH);
		const float tMid = thresholds[(thresholds.size() - 1) / 2];
		const auto tFull = std::chrono::steady_clock::now();
		extg2.ComputeMinMapFromGradient(tMid);
		printf("extnet2d_ms phase=extB_full_rebuild t=%f ms=%.3f\n", tMid, msSince(tFull));
	}

	// --- Phase 3b: Path B2 — extremum network built WITHOUT V-path traces.
	// The saddle -> extremum connections are exactly lookups into the flow
	// terminal maps (a 1-saddle's two facet vertices flow to its two minima; its
	// two cofacet quads flow to its two maxima), so the arcs come from the base
	// labeling we already have instead of walking gradient paths per saddle.
	UFGraphType minG2(meshfunc, mesh, UFGraphType::MINIMAL);
	UFGraphType maxG2(meshfunc, mesh, UFGraphType::MAXIMAL);
	{
		minG2.SetRetainArcs(true);
		maxG2.SetRetainArcs(true);
		const int64_t cellRow = 2 * X - 1;
		const int64_t cellRows = 2 * Y - 1;
		const int64_t numCells = cellRow * cellRows;

		const auto tScan = std::chrono::steady_clock::now();
		std::vector<int64_t> critEdges;
#pragma omp parallel
		{
			std::vector<int64_t> local;
#pragma omp for schedule(static)
			for (int64_t c = 0; c < numCells; ++c) {
				const int64_t cx = c % cellRow;
				const int64_t cy = c / cellRow;
				if (((cx ^ cy) & 1) == 1 && grad->getCritical(c)) local.push_back(c);
			}
#pragma omp critical
			critEdges.insert(critEdges.end(), local.begin(), local.end());
		}
		const double scanMs = msSince(tScan);

		const auto tB2 = std::chrono::steady_clock::now();
		for (size_t r = 0; r < ascRegions.terminalCells.size(); ++r)
			minG2.AddNode((INDEX_TYPE)ascRegions.terminalCells[r]);
		for (size_t r = 0; r < dscRegions.terminalCells.size(); ++r)
			maxG2.AddNode((INDEX_TYPE)dscRegions.terminalCells[r]);
		for (size_t i = 0; i < critEdges.size(); ++i) {
			const int64_t s = critEdges[i];
			const int64_t cx = s % cellRow;
			const int64_t cy = s / cellRow;
			const bool xEdge = ((cx & 1) == 1);
			// min side: the two facet vertices' flow terminals
			const int64_t v1 = xEdge ? (s - 1) : (s - cellRow);
			const int64_t v2 = xEdge ? (s + 1) : (s + cellRow);
			const int64_t m1 = ascTerm[(size_t)((v1 % cellRow) / 2 + (v1 / cellRow) / 2 * X)];
			const int64_t m2 = ascTerm[(size_t)((v2 % cellRow) / 2 + (v2 / cellRow) / 2 * X)];
			if (m1 != m2) minG2.AddArc((INDEX_TYPE)m1, (INDEX_TYPE)m2, (INDEX_TYPE)s);
			// max side: the two cofacet quads' flow terminals; skip boundary
			// saddles exactly like trace_up_2saddle does
			const bool interior = xEdge ? (cy > 0 && cy < cellRows - 1)
				: (cx > 0 && cx < cellRow - 1);
			if (!interior) continue;
			const int64_t q1 = xEdge ? (s - cellRow) : (s - 1);
			const int64_t q2 = xEdge ? (s + cellRow) : (s + 1);
			const int64_t M1 = dscTerm[(size_t)((q1 % cellRow - 1) / 2 + (q1 / cellRow - 1) / 2 * (X - 1))];
			const int64_t M2 = dscTerm[(size_t)((q2 % cellRow - 1) / 2 + (q2 / cellRow - 1) / 2 * (X - 1))];
			if (M1 != M2) maxG2.AddArc((INDEX_TYPE)M1, (INDEX_TYPE)M2, (INDEX_TYPE)s);
		}
		const double b2Ms = msSince(tB2);
		printf("extnet2d_ms phase=extB2_scan ms=%.3f\n", scanMs);
		printf("extnet2d_ms phase=extB2_build_from_terminals ms=%.3f\n", b2Ms);
		printf("extnet2d_extB2 saddles=%lld minNodes=%d maxNodes=%d\n",
			(long long)critEdges.size(), minG2.NumExtrema(), maxG2.NumExtrema());
	}

	// --- per-threshold sweep
	std::vector<int> labelA;
	std::vector<int64_t> labelB;
	bool b2AllMatch = true;
	for (size_t k = 0; k < thresholds.size(); ++k) {
		const float t = thresholds[k];

		// Path A serial: select + walk-down remap, per direction
		long long confAsc = 0, confDsc = 0;
		const auto tA1 = std::chrono::steady_clock::now();
		serialMSC.SetSelectPersAbs(t);
		std::vector<int> labelAasc;
		const size_t livingMins = pathALabels(&serialMSC, ctxAsc, true, 0, labelAasc, confAsc);
		const double aAscMs = msSince(tA1);
		const auto tA2 = std::chrono::steady_clock::now();
		std::vector<int> labelAdsc;
		const size_t livingMaxs = pathALabels(&serialMSC, ctxDsc, false, 2, labelAdsc, confDsc);
		const double aDscMs = msSince(tA2);
		if (confAsc != 0 || confDsc != 0) {
			printf("extnet2d_warn k=%lld walkdown_conflicts asc=%lld dsc=%lld\n",
				(long long)(k + 1), confAsc, confDsc);
		}

		// Path A partitioned: select only (living-count cross-check)
		const auto tAp = std::chrono::steady_clock::now();
		reconciled->SetSelectPersAbs(t);
		long long partLivingMins = 0, partLivingMaxs = 0;
		{
			PartitionedPipelineType::ReconciledGlobalMsc* rm = reconciled.get();
			SerialMscType::LivingNodesIterator nit(rm);
			for (nit.begin(); nit.valid(); nit.advance()) {
				const int d = rm->getNode(nit.value()).dim;
				if (d == 0) partLivingMins++;
				else if (d == 2) partLivingMaxs++;
			}
		}
		const double aPartMs = msSince(tAp);

		// Path B: reset + drain per graph, best of runB_repeat
		double bMinMs = 1e30, bMaxMs = 1e30;
		for (int rep = 0; rep < runB_repeat; ++rep) {
			const auto tB1 = std::chrono::steady_clock::now();
			extg.mMinGraph->ResetSimplification();
			extg.mMinGraph->SimplifyToThreshold(t);
			const double m1 = msSince(tB1);
			if (m1 < bMinMs) bMinMs = m1;
			const auto tB2 = std::chrono::steady_clock::now();
			extg.mMaxGraph->ResetSimplification();
			extg.mMaxGraph->SimplifyToThreshold(t);
			const double m2 = msSince(tB2);
			if (m2 < bMaxMs) bMaxMs = m2;
		}

		// Path B2 (terminal-built graphs): identical drains
		double b2MinMs = 1e30, b2MaxMs = 1e30;
		for (int rep = 0; rep < runB_repeat; ++rep) {
			const auto tB1 = std::chrono::steady_clock::now();
			minG2.ResetSimplification();
			minG2.SimplifyToThreshold(t);
			const double m1 = msSince(tB1);
			if (m1 < b2MinMs) b2MinMs = m1;
			const auto tB2 = std::chrono::steady_clock::now();
			maxG2.ResetSimplification();
			maxG2.SimplifyToThreshold(t);
			const double m2 = msSince(tB2);
			if (m2 < b2MaxMs) b2MaxMs = m2;
		}

		printf("extnet2d_kms k=%lld t=%f A_asc_ms=%.3f A_dsc_ms=%.3f A_part_select_ms=%.3f B_min_ms=%.3f B_max_ms=%.3f B2_min_ms=%.3f B2_max_ms=%.3f\n",
			(long long)(k + 1), t, aAscMs, aDscMs, aPartMs, bMinMs, bMaxMs, b2MinMs, b2MaxMs);

		// Path B labels + agreement
		for (int dir = 0; dir < 2; ++dir) {
			const bool ascending = (dir == 0);
			const BaseRegions& regions = ascending ? ascRegions : dscRegions;
			UFGraphType* g = ascending ? extg.mMinGraph : extg.mMaxGraph;
			labelA = ascending ? labelAasc : labelAdsc;
			labelB.assign(regions.terminalCells.size(), -1);
			bool missing = false;
			for (size_t r = 0; r < regions.terminalCells.size(); ++r) {
				std::unordered_map<INDEX_TYPE, INT_TYPE>::const_iterator it =
					g->mCellIndexToListIdMap.find((INDEX_TYPE)regions.terminalCells[r]);
				if (it == g->mCellIndexToListIdMap.end()) { missing = true; continue; }
				labelB[r] = (int64_t)g->TopoIndex(g->Representative(it->second));
			}
			if (missing) {
				printf("extnet2d_warn k=%lld dir=%s missing_uf_terminals\n",
					(long long)(k + 1), ascending ? "asc" : "dsc");
			}
			// gate: the terminal-built graph must simplify to identical labels
			{
				UFGraphType* g2 = ascending ? &minG2 : &maxG2;
				bool b2same = true;
				for (size_t r = 0; r < regions.terminalCells.size(); ++r) {
					std::unordered_map<INDEX_TYPE, INT_TYPE>::const_iterator it =
						g2->mCellIndexToListIdMap.find((INDEX_TYPE)regions.terminalCells[r]);
					const int64_t lb2 = (it == g2->mCellIndexToListIdMap.end()) ? -1
						: (int64_t)g2->TopoIndex(g2->Representative(it->second));
					if (lb2 != labelB[r]) { b2same = false; break; }
				}
				if (!b2same) {
					printf("extnet2d_warn k=%lld dir=%s extB2_label_mismatch\n",
						(long long)(k + 1), ascending ? "asc" : "dsc");
					b2AllMatch = false;
				}
			}
			const Agreement ag = computeAgreement(labelA, labelB, regions.regionWeight);
			printf("extnet2d_agree k=%lld dir=%s countA=%lld countB=%lld livingA=%lld partLivingA=%lld purityAB=%.5f purityBA=%.5f mutual=%.5f weight=%lld\n",
				(long long)(k + 1), ascending ? "asc" : "dsc",
				ag.countA, ag.countB,
				(long long)(ascending ? livingMins : livingMaxs),
				(ascending ? partLivingMins : partLivingMaxs),
				ag.purityAB, ag.purityBA, ag.mutual, ag.weight);
		}
	}

	printf("extnet2d_gate name=extB2_matches_extB pass=%d\n", b2AllMatch ? 1 : 0);
	printf("extnet2d_done\n");
	return b2AllMatch ? 0 : 6;
}
