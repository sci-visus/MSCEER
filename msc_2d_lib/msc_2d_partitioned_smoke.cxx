#include "gi_discrete_gradient_computer.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_partitioned_topological_regular_grid.h"

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <utility>
#include <vector>

typedef GInt::MorseSmaleComplexBasic<float, GInt::Accurate2D::MeshType, GInt::Accurate2D::MeshFuncType, GInt::Accurate2D::GradType> SerialMscType;
typedef GInt::MorseSmaleComplexPartitioned<float, GInt::Accurate2D::MeshType, GInt::Accurate2D::MeshFuncType, GInt::Accurate2D::GradType> PartitionedPipelineType;

struct AliveArcBlockStats {
	size_t alive_total;
	size_t alive_old_total;
	size_t alive_under_pers;
	size_t alive_old_under_pers;
	size_t blocked_boundary;
	size_t blocked_multiplicity;
	size_t blocked_other;
	size_t valid_under_pers;
	size_t old_blocked_boundary;
	size_t old_blocked_multiplicity;
	size_t old_blocked_other;
	size_t old_valid_under_pers;
	AliveArcBlockStats() :
		alive_total(0),
		alive_old_total(0),
		alive_under_pers(0),
		alive_old_under_pers(0),
		blocked_boundary(0),
		blocked_multiplicity(0),
		blocked_other(0),
		valid_under_pers(0),
		old_blocked_boundary(0),
		old_blocked_multiplicity(0),
		old_blocked_other(0),
		old_valid_under_pers(0) {}
};

static AliveArcBlockStats analyzeAliveArcBlockingReasons(SerialMscType* msc, float persLimit) {
	AliveArcBlockStats stats;
	std::map<std::pair<INT_TYPE, INT_TYPE>, int> multiplicity;
	for (INT_TYPE aid = 0; aid < msc->numArcs(); ++aid) {
		if (!msc->isArcAlive(aid)) continue;
		const GInt::arc<float>& a = msc->getArc(aid);
		++multiplicity[std::make_pair(a.lower, a.upper)];
		stats.alive_total++;
		if (a.created == 0) stats.alive_old_total++;
	}

	for (INT_TYPE aid = 0; aid < msc->numArcs(); ++aid) {
		if (!msc->isArcAlive(aid)) continue;
		const GInt::arc<float>& a = msc->getArc(aid);
		if (a.persistence > persLimit) continue;
		const bool is_old = (a.created == 0);
		stats.alive_under_pers++;
		if (is_old) stats.alive_old_under_pers++;

		const bool boundary_mismatch = (msc->getNode(a.lower).boundary != msc->getNode(a.upper).boundary);
		const int mult = multiplicity[std::make_pair(a.lower, a.upper)];
		if (boundary_mismatch) {
			stats.blocked_boundary++;
			if (is_old) stats.old_blocked_boundary++;
		}
		else if (mult != 1) {
			stats.blocked_multiplicity++;
			if (is_old) stats.old_blocked_multiplicity++;
		}
		else {
			// If these survive under persLimit but still look valid, this is a useful red flag.
			stats.valid_under_pers++;
			if (is_old) stats.old_valid_under_pers++;
		}
	}

	// "other" is kept as explicit residual for sanity accounting.
	stats.blocked_other =
		stats.alive_under_pers - stats.blocked_boundary - stats.blocked_multiplicity - stats.valid_under_pers;
	stats.old_blocked_other =
		stats.alive_old_under_pers - stats.old_blocked_boundary - stats.old_blocked_multiplicity - stats.old_valid_under_pers;
	return stats;
}

static std::map<std::pair<INDEX_TYPE, INDEX_TYPE>, int> livingEndpointHistogram(SerialMscType* msc) {
	std::map<std::pair<INDEX_TYPE, INDEX_TYPE>, int> hist;
	for (INT_TYPE aid = 0; aid < msc->numArcs(); ++aid) {
		if (!msc->isArcAlive(aid)) continue;
		const GInt::arc<float>& a = msc->getArc(aid);
		const INDEX_TYPE lowerCell = msc->getNode(a.lower).cellindex;
		const INDEX_TYPE upperCell = msc->getNode(a.upper).cellindex;
		++hist[std::make_pair(lowerCell, upperCell)];
	}
	return hist;
}

static size_t countLivingNodes(SerialMscType* msc) {
	size_t count = 0;
	for (INT_TYPE nid = 0; nid < msc->numNodes(); ++nid) {
		if (msc->isNodeAlive(nid)) count++;
	}
	return count;
}

static size_t countLivingArcs(SerialMscType* msc) {
	size_t count = 0;
	for (INT_TYPE aid = 0; aid < msc->numArcs(); ++aid) {
		if (msc->isArcAlive(aid)) count++;
	}
	return count;
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

int main(int argc, char** argv) {
	printf("[smoke] start main\n");
	const int rows = 2048;
	const int cols = 2048;
	int partitions = 4;
	if (argc > 1) {
		partitions = std::atoi(argv[1]);
	}
	if (!GInt::PartitionedTopologicalRegularGrid2D::IsSupportedPartitionCount(partitions)) {
		std::cerr << "Unsupported partition/thread count: " << partitions
			<< ". Supported values are {1,2,3,4,6,8,9,12,16}." << std::endl;
		return 2;
	}
	printf("[smoke] config rows=%d cols=%d partitions=%d\n", rows, cols, partitions);

	std::vector<float> field((size_t)rows * (size_t)cols, 0.0f);
	printf("[smoke] allocated field size=%llu\n", (unsigned long long)field.size());
	std::mt19937 rng(123456789u);
	std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

	std::vector<float> coarseNoise((size_t)rows * (size_t)cols, 0.0f);
	std::vector<float> detailNoise((size_t)rows * (size_t)cols, 0.0f);

	printf("[smoke] generating base random fields\n");
	for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < cols; ++c) {
			coarseNoise[(size_t)r * (size_t)cols + (size_t)c] = dist(rng);
			detailNoise[(size_t)r * (size_t)cols + (size_t)c] = dist(rng);
		}
	}
	printf("[smoke] smoothing coarse field (30 iterations)\n");
	smoothFieldInPlace(coarseNoise, rows, cols, 30);
	printf("[smoke] smoothing detail field (2 iterations)\n");
	smoothFieldInPlace(detailNoise, rows, cols, 2);

	const float coarseScale = 100.0f;
	const float detailScale = 0.05f * coarseScale;
	printf("[smoke] composing multiscale field coarseScale=%f detailScale=%f\n", coarseScale, detailScale);
	for (size_t i = 0; i < field.size(); ++i) {
		field[i] = coarseScale * coarseNoise[i] + detailScale * detailNoise[i];
	}
	// Ring masking from earlier experiments intentionally disabled for now:
	// const float cx = 0.5f * static_cast<float>(cols);
	// const float cy = 0.5f * static_cast<float>(rows);
	// const float radius = 400.0f;
	// const float r2 = radius * radius;
	// if (dx * dx + dy * dy > r2) { v = 0.0f; }

	printf("[smoke] creating DiscreteGradientBuilder\n");
	GInt::Accurate2D::DiscreteGradientBuilder dgb;
	printf("[smoke] calling SetFloadArrayAndDims\n");
	dgb.SetFloadArrayAndDims(cols, rows, field.data());
	printf("[smoke] calling SetNeededAccuracy\n");
	dgb.SetNeededAccuracy(true, true);
	printf("[smoke] calling SetParallelism\n");
	dgb.SetParallelism(partitions);
	printf("[smoke] calling ComputeDiscreteGradient\n");
	const auto tDiscreteStart = std::chrono::steady_clock::now();
	dgb.ComputeDiscreteGradient();
	const auto tDiscreteEnd = std::chrono::steady_clock::now();

	printf("[smoke] fetching mesh/func/grad pointers\n");
	auto* mesh = dgb.GetTopoMesh();
	auto* meshfunc = dgb.GetMeshFunc();
	auto* grad = dgb.GetGrad();
	auto* gridfunc = dgb.GetGridFunc();
	printf("[smoke] pointers mesh=%p meshfunc=%p grad=%p\n", (void*)mesh, (void*)meshfunc, (void*)grad);

	const float maxval = gridfunc->GetMaxValue();
	const float minval = gridfunc->GetMinValue();
	const float localPers = 0.01f * (maxval - minval);
	const float globalPers = 0.05f * (maxval - minval);
	printf("[smoke] localPers=%f globalPers=%f\n", localPers, globalPers);

	// Measure wall-clock task times requested in benchmark table.
	const long long discreteGradMs =
		std::chrono::duration_cast<std::chrono::milliseconds>(tDiscreteEnd - tDiscreteStart).count();

	printf("[smoke] constructing serial MSC\n");
	const auto tSerialConstructStart = std::chrono::steady_clock::now();
	SerialMscType serialMSC(grad, mesh, meshfunc);
	printf("[smoke] serial SetBuildArcGeometry\n");
	serialMSC.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
	serialMSC.ComputeFromGrad();
	const INT_TYPE serialPreSimplifyNodes = serialMSC.numNodes();
	const INT_TYPE serialPreSimplifyArcs = serialMSC.numArcs();
	const auto tSerialConstructEnd = std::chrono::steady_clock::now();
	const auto tSerialSimplifyStart = std::chrono::steady_clock::now();
	serialMSC.ComputeHierarchy(globalPers);
	serialMSC.SetSelectPersAbs(globalPers);
	const auto tSerialSimplifyEnd = std::chrono::steady_clock::now();
	const long long serialConstructMs =
		std::chrono::duration_cast<std::chrono::milliseconds>(tSerialConstructEnd - tSerialConstructStart).count();
	const long long serialSimplifyMs =
		std::chrono::duration_cast<std::chrono::milliseconds>(tSerialSimplifyEnd - tSerialSimplifyStart).count();
	printf("[smoke] serial MSC computed nodes=%d arcs=%d\n", (int)serialMSC.numNodes(), (int)serialMSC.numArcs());

	printf("[smoke] constructing partitioned pipeline\n");
	PartitionedPipelineType::TimingBreakdown partitionTimings;
	PartitionedPipelineType partitioned(grad, mesh, meshfunc);
	std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
		partitioned.BuildPartitionLocalMSCs(partitions, localPers, &partitionTimings);
	printf("[smoke] localResults=%llu\n", (unsigned long long)localResults.size());
	long long localPreSimplifyNodesSum = 0;
	long long localPreSimplifyArcsSum = 0;
	long long localPreExchangeFrozenNodesSum = 0;
	long long localPostExchangeFrozenNodesSum = 0;
	INT_TYPE localPreExchangeFrozenNodesMax = 0;
	INT_TYPE localPostExchangeFrozenNodesMax = 0;
	for (size_t i = 0; i < localResults.size(); ++i) {
		localPreSimplifyNodesSum += static_cast<long long>(localResults[i].pre_simplify_num_nodes);
		localPreSimplifyArcsSum += static_cast<long long>(localResults[i].pre_simplify_num_arcs);
		localPreExchangeFrozenNodesSum += static_cast<long long>(localResults[i].pre_exchange_frozen_nodes);
		localPostExchangeFrozenNodesSum += static_cast<long long>(localResults[i].post_exchange_frozen_nodes);
		if (localResults[i].pre_exchange_frozen_nodes > localPreExchangeFrozenNodesMax) {
			localPreExchangeFrozenNodesMax = localResults[i].pre_exchange_frozen_nodes;
		}
		if (localResults[i].post_exchange_frozen_nodes > localPostExchangeFrozenNodesMax) {
			localPostExchangeFrozenNodesMax = localResults[i].post_exchange_frozen_nodes;
		}
	}
	printf("[smoke] constructing partition mesh\n");
	GInt::PartitionedTopologicalRegularGrid2D partitionMesh(mesh, partitions);
	std::unique_ptr<PartitionedPipelineType::ReconciledGlobalMsc> reconciled =
		partitioned.BuildReconciledGlobalBase(partitionMesh, localResults, &partitionTimings);
	printf("[smoke] reconciled base nodes=%d arcs=%d\n",
		(int)reconciled->numNodes(), (int)reconciled->numArcs());
	const auto tPartitionedGlobalSimplifyStart = std::chrono::steady_clock::now();
	reconciled->ComputeHierarchy(globalPers);
	reconciled->SetSelectPersAbs(globalPers);
	const auto tPartitionedGlobalSimplifyEnd = std::chrono::steady_clock::now();
	partitionTimings.global_simplify_ms =
		std::chrono::duration_cast<std::chrono::milliseconds>(tPartitionedGlobalSimplifyEnd - tPartitionedGlobalSimplifyStart).count();
	const long long partitionedTotalMs =
		partitionTimings.local_stage_total_ms + partitionTimings.reconcile_ms + partitionTimings.global_simplify_ms;

	printf("[smoke] counting living nodes/arcs serial\n");
	const size_t serialLivingNodes = countLivingNodes(&serialMSC);
	const size_t serialLivingArcs = countLivingArcs(&serialMSC);
	printf("[smoke] counting living nodes/arcs reconciled\n");
	const size_t reconciledLivingNodes = countLivingNodes(reconciled.get());
	const size_t reconciledLivingArcs = countLivingArcs(reconciled.get());

	printf("[smoke] computing endpoint histograms\n");
	const auto serialHist = livingEndpointHistogram(&serialMSC);
	const auto reconciledHist = livingEndpointHistogram(reconciled.get());
	printf("[smoke] histogram sizes serial=%llu reconciled=%llu\n",
		(unsigned long long)serialHist.size(), (unsigned long long)reconciledHist.size());

	const AliveArcBlockStats serialBlockStats = analyzeAliveArcBlockingReasons(&serialMSC, globalPers);
	const AliveArcBlockStats reconciledBlockStats = analyzeAliveArcBlockingReasons(reconciled.get(), globalPers);
	const SerialMscType::cancel_branch_stats serialBranchStats = serialMSC.GetCancelBranchStats();
	const SerialMscType::cancel_branch_stats reconciledBranchStats = reconciled->GetCancelBranchStats();
	const double serialCancels = static_cast<double>(serialMSC.GetCancellationRecords().size());
	const double reconciledCancels = static_cast<double>(reconciled->GetCancellationRecords().size());
	const double serialTotalMergeUpdates =
		static_cast<double>(serialBranchStats.descending_merge_updates + serialBranchStats.ascending_merge_updates);
	const double serialTotalDestroyOnly =
		static_cast<double>(serialBranchStats.lower_destroy_only + serialBranchStats.upper_destroy_only);
	const double reconciledTotalMergeUpdates =
		static_cast<double>(reconciledBranchStats.descending_merge_updates + reconciledBranchStats.ascending_merge_updates);
	const double reconciledTotalDestroyOnly =
		static_cast<double>(reconciledBranchStats.lower_destroy_only + reconciledBranchStats.upper_destroy_only);
	const double serialMergeUpdatesPerCancel = (serialCancels > 0.0) ? (serialTotalMergeUpdates / serialCancels) : 0.0;
	const double serialDestroyOnlyPerCancel = (serialCancels > 0.0) ? (serialTotalDestroyOnly / serialCancels) : 0.0;
	const double serialDestroyToMergeRatio = (serialTotalMergeUpdates > 0.0) ? (serialTotalDestroyOnly / serialTotalMergeUpdates) : 0.0;
	const double reconciledMergeUpdatesPerCancel = (reconciledCancels > 0.0) ? (reconciledTotalMergeUpdates / reconciledCancels) : 0.0;
	const double reconciledDestroyOnlyPerCancel = (reconciledCancels > 0.0) ? (reconciledTotalDestroyOnly / reconciledCancels) : 0.0;
	const double reconciledDestroyToMergeRatio = (reconciledTotalMergeUpdates > 0.0) ? (reconciledTotalDestroyOnly / reconciledTotalMergeUpdates) : 0.0;
	std::cout << "alive_block_reasons"
		<< " mode=serial"
		<< " pers=" << globalPers
		<< " alive_total=" << serialBlockStats.alive_total
		<< " alive_old_total=" << serialBlockStats.alive_old_total
		<< " alive_under_pers=" << serialBlockStats.alive_under_pers
		<< " blocked_boundary=" << serialBlockStats.blocked_boundary
		<< " blocked_multiplicity=" << serialBlockStats.blocked_multiplicity
		<< " blocked_other=" << serialBlockStats.blocked_other
		<< " valid_under_pers=" << serialBlockStats.valid_under_pers
		<< " old_alive_under_pers=" << serialBlockStats.alive_old_under_pers
		<< " old_blocked_boundary=" << serialBlockStats.old_blocked_boundary
		<< " old_blocked_multiplicity=" << serialBlockStats.old_blocked_multiplicity
		<< " old_blocked_other=" << serialBlockStats.old_blocked_other
		<< " old_valid_under_pers=" << serialBlockStats.old_valid_under_pers
		<< std::endl;
	std::cout << "alive_block_reasons"
		<< " mode=reconciled"
		<< " pers=" << globalPers
		<< " alive_total=" << reconciledBlockStats.alive_total
		<< " alive_old_total=" << reconciledBlockStats.alive_old_total
		<< " alive_under_pers=" << reconciledBlockStats.alive_under_pers
		<< " blocked_boundary=" << reconciledBlockStats.blocked_boundary
		<< " blocked_multiplicity=" << reconciledBlockStats.blocked_multiplicity
		<< " blocked_other=" << reconciledBlockStats.blocked_other
		<< " valid_under_pers=" << reconciledBlockStats.valid_under_pers
		<< " old_alive_under_pers=" << reconciledBlockStats.alive_old_under_pers
		<< " old_blocked_boundary=" << reconciledBlockStats.old_blocked_boundary
		<< " old_blocked_multiplicity=" << reconciledBlockStats.old_blocked_multiplicity
		<< " old_blocked_other=" << reconciledBlockStats.old_blocked_other
		<< " old_valid_under_pers=" << reconciledBlockStats.old_valid_under_pers
		<< std::endl;
	std::cout << "cancel_branch_stats"
		<< " mode=serial"
		<< " desc_merge_updates=" << serialBranchStats.descending_merge_updates
		<< " asc_merge_updates=" << serialBranchStats.ascending_merge_updates
		<< " lower_destroy_only=" << serialBranchStats.lower_destroy_only
		<< " upper_destroy_only=" << serialBranchStats.upper_destroy_only
		<< std::endl;
	std::cout << "cancel_branch_stats"
		<< " mode=reconciled"
		<< " desc_merge_updates=" << reconciledBranchStats.descending_merge_updates
		<< " asc_merge_updates=" << reconciledBranchStats.ascending_merge_updates
		<< " lower_destroy_only=" << reconciledBranchStats.lower_destroy_only
		<< " upper_destroy_only=" << reconciledBranchStats.upper_destroy_only
		<< std::endl;
	std::cout << "cancel_branch_norm"
		<< " mode=serial"
		<< " cancels=" << serialCancels
		<< " merge_updates_per_cancel=" << serialMergeUpdatesPerCancel
		<< " destroy_only_per_cancel=" << serialDestroyOnlyPerCancel
		<< " destroy_to_merge_ratio=" << serialDestroyToMergeRatio
		<< std::endl;
	std::cout << "cancel_branch_norm"
		<< " mode=reconciled"
		<< " cancels=" << reconciledCancels
		<< " merge_updates_per_cancel=" << reconciledMergeUpdatesPerCancel
		<< " destroy_only_per_cancel=" << reconciledDestroyOnlyPerCancel
		<< " destroy_to_merge_ratio=" << reconciledDestroyToMergeRatio
		<< std::endl;

	std::cout << "partitioned_smoke"
		<< " partitions=" << partitions
		<< " localPers=" << localPers
		<< " globalPers=" << globalPers
		<< " serial_pre_nodes=" << serialPreSimplifyNodes
		<< " serial_pre_arcs=" << serialPreSimplifyArcs
		<< " local_pre_nodes_sum=" << localPreSimplifyNodesSum
		<< " local_pre_arcs_sum=" << localPreSimplifyArcsSum
		<< " local_pre_exchange_frozen_nodes_sum=" << localPreExchangeFrozenNodesSum
		<< " local_post_exchange_frozen_nodes_sum=" << localPostExchangeFrozenNodesSum
		<< " local_pre_exchange_frozen_nodes_max=" << localPreExchangeFrozenNodesMax
		<< " local_post_exchange_frozen_nodes_max=" << localPostExchangeFrozenNodesMax
		<< " serial_nodes=" << serialLivingNodes
		<< " serial_arcs=" << serialLivingArcs
		<< " reconciled_nodes=" << reconciledLivingNodes
		<< " reconciled_arcs=" << reconciledLivingArcs
		<< " same_endpoint_histogram=" << (serialHist == reconciledHist ? "true" : "false")
		<< std::endl;
	std::cout << "benchmark_ms"
		<< " partitions=" << partitions
		<< " discrete_grad=" << discreteGradMs
		<< " serial_construct=" << serialConstructMs
		<< " serial_simplify=" << serialSimplifyMs
		<< " partition_local_total=" << partitionTimings.local_stage_total_ms
		<< " partition_local_construct=" << partitionTimings.local_build_ms
		<< " partition_local_simplify=" << partitionTimings.local_simplify_ms
		<< " partition_reconcile=" << partitionTimings.reconcile_ms
		<< " partition_global_simplify=" << partitionTimings.global_simplify_ms
		<< " partition_total=" << partitionedTotalMs
		<< std::endl;
	printf("[smoke] finished\n");

	return 0;
}
