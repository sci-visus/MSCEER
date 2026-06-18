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

	std::cout << "partitioned_smoke"
		<< " partitions=" << partitions
		<< " localPers=" << localPers
		<< " globalPers=" << globalPers
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
