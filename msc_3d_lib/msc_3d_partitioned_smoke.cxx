#include "msc_3d_lib.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <utility>
#include <vector>

#include "gi_discrete_gradient_labeling.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_regular_grid.h"

typedef GInt::RegularGrid3D GridType;
typedef GInt::RegularGridTrilinearFunction GridFuncType;
typedef GInt::TopologicalRegularGrid3D MeshType;
typedef GInt::DiscreteGradientLabeling<MeshType> GradType;
typedef GInt::TopologicalExplicitDenseMeshFunction<MeshType, float> MeshFuncType;
typedef GInt::MorseSmaleComplexBasic<float, MeshType, MeshFuncType, GradType> SerialMscType;
typedef GInt::MorseSmaleComplexPartitioned<float, MeshType, MeshFuncType, GradType> PartitionedPipelineType;

static size_t countLivingNodes(SerialMscType* msc) {
    size_t count = 0;
    for (INT_TYPE nid = 0; nid < msc->numNodes(); ++nid) {
        if (msc->isNodeAlive(nid)) ++count;
    }
    return count;
}

static size_t countLivingArcs(SerialMscType* msc) {
    size_t count = 0;
    for (INT_TYPE aid = 0; aid < msc->numArcs(); ++aid) {
        if (msc->isArcAlive(aid)) ++count;
    }
    return count;
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

int main(int argc, char** argv) {
    const int xdim = 40;
    const int ydim = 24;
    const int zdim = 20;
    int requestedPartitions = 12;
    int requestedThreads = 5;
    if (argc > 1) requestedPartitions = std::atoi(argv[1]);
    if (argc > 2) requestedThreads = std::atoi(argv[2]);
    if (requestedPartitions < 1) requestedPartitions = 1;
    if (requestedThreads < 1) requestedThreads = 1;

    std::vector<float> field((size_t)xdim * (size_t)ydim * (size_t)zdim, 0.0f);
    std::mt19937 rng(123456789u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (size_t i = 0; i < field.size(); ++i) {
        field[i] = dist(rng);
    }

    GridType grid(GInt::Vec3i(xdim, ydim, zdim), GInt::Vec3b(false, false, false));
    GridFuncType gridfunc(&grid, field.data());
    MeshType mesh(&grid);
    MeshFuncType meshfunc;
    meshfunc.setMeshAndAllocate(&mesh);
    for (INDEX_TYPE cid = 0; cid < mesh.numCells(); ++cid) {
        GInt::Vec3l c;
        mesh.cellid2Coords(cid, c);
        const INDEX_TYPE vx = std::min(static_cast<INDEX_TYPE>(xdim - 1), c[0] / 2);
        const INDEX_TYPE vy = std::min(static_cast<INDEX_TYPE>(ydim - 1), c[1] / 2);
        const INDEX_TYPE vz = std::min(static_cast<INDEX_TYPE>(zdim - 1), c[2] / 2);
        const size_t linear = static_cast<size_t>(vx)
            + static_cast<size_t>(xdim) * static_cast<size_t>(vy)
            + static_cast<size_t>(xdim) * static_cast<size_t>(ydim) * static_cast<size_t>(vz);
        meshfunc.setCellValue(cid, field[linear]);
    }
    GradType grad(&mesh);
    grad.ClearAllGradient();
    for (INDEX_TYPE cid = 0; cid < mesh.numCells(); ++cid) {
        grad.setAssigned(cid, 1);
        grad.setCritical(cid, true);
    }

    const float maxval = gridfunc.GetMaxValue();
    const float minval = gridfunc.GetMinValue();
    const float localPers = 0.0f;
    const float globalPers = 0.1f * (maxval - minval);

    SerialMscType serial(&grad, &mesh, &meshfunc);
    serial.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
    serial.ComputeFromGrad();
    serial.ComputeHierarchy(globalPers);
    serial.SetSelectPersAbs(globalPers);

    PartitionedPipelineType partitioned(&grad, &mesh, &meshfunc);
    PartitionedPipelineType::TimingBreakdown timings;
    std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
        partitioned.BuildPartitionLocalMSCs(requestedPartitions, localPers, &timings, requestedThreads);
    GInt::PartitionedTopologicalRegularGrid3D partitionGrid(&mesh, requestedPartitions);
    std::unique_ptr<PartitionedPipelineType::ReconciledGlobalMsc> reconciled =
        partitioned.BuildReconciledGlobalBase(partitionGrid, localResults, &timings);
    reconciled->ComputeHierarchy(globalPers);
    reconciled->SetSelectPersAbs(globalPers);

    long long preFrozenSum = 0;
    long long postFrozenSum = 0;
    std::vector<int> buildAssignments((size_t)requestedThreads, 0);
    for (size_t i = 0; i < localResults.size(); ++i) {
        preFrozenSum += localResults[i].pre_exchange_frozen_nodes;
        postFrozenSum += localResults[i].post_exchange_frozen_nodes;
        if (localResults[i].build_worker_slot >= 0 && localResults[i].build_worker_slot < requestedThreads) {
            buildAssignments[(size_t)localResults[i].build_worker_slot]++;
        }
    }

    const size_t serialNodes = countLivingNodes(&serial);
    const size_t serialArcs = countLivingArcs(&serial);
    const size_t reconciledNodes = countLivingNodes(reconciled.get());
    const size_t reconciledArcs = countLivingArcs(reconciled.get());
    const auto serialHist = livingEndpointHistogram(&serial);
    const auto reconciledHist = livingEndpointHistogram(reconciled.get());
    const GInt::Vec3l splits = partitionGrid.partition_splits();

    std::cout << "partitioned_3d_smoke"
              << " requested_partitions=" << requestedPartitions
              << " effective_partitions=" << partitionGrid.num_partitions()
              << " requested_threads=" << requestedThreads
              << " inferred_splits=(" << splits[0] << "," << splits[1] << "," << splits[2] << ")"
              << " local_pre_exchange_frozen_nodes_sum=" << preFrozenSum
              << " local_post_exchange_frozen_nodes_sum=" << postFrozenSum
              << " serial_nodes=" << serialNodes
              << " serial_arcs=" << serialArcs
              << " reconciled_nodes=" << reconciledNodes
              << " reconciled_arcs=" << reconciledArcs
              << " same_endpoint_histogram=" << (serialHist == reconciledHist ? "true" : "false")
              << std::endl;

    std::cout << "partitioned_3d_thread_assignment";
    for (int t = 0; t < requestedThreads; ++t) {
        std::cout << " t" << t << "=" << buildAssignments[(size_t)t];
    }
    std::cout << std::endl;

    if (partitionGrid.num_partitions() == 1) {
        if (serialNodes != reconciledNodes || serialArcs != reconciledArcs || serialHist != reconciledHist) {
            std::cerr << "partitioned_3d warning: single-partition parity mismatch"
                      << " serial_nodes=" << serialNodes
                      << " serial_arcs=" << serialArcs
                      << " reconciled_nodes=" << reconciledNodes
                      << " reconciled_arcs=" << reconciledArcs
                      << std::endl;
        }
    }
    return 0;
}
