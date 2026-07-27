#include <algorithm>
#include <chrono>
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

static void smoothField3DInPlace(std::vector<float>& field, int nx, int ny, int nz, int iterations) {
    if (iterations <= 0) return;
    std::vector<float> scratch(field.size(), 0.0f);
    const size_t nxSz = static_cast<size_t>(nx);
    const size_t nySz = static_cast<size_t>(ny);
    for (int it = 0; it < iterations; ++it) {
        for (int z = 0; z < nz; ++z) {
            const int z0 = (z > 0) ? (z - 1) : z;
            const int z1 = (z < nz - 1) ? (z + 1) : z;
            for (int y = 0; y < ny; ++y) {
                const int y0 = (y > 0) ? (y - 1) : y;
                const int y1 = (y < ny - 1) ? (y + 1) : y;
                for (int x = 0; x < nx; ++x) {
                    const int x0 = (x > 0) ? (x - 1) : x;
                    const int x1 = (x < nx - 1) ? (x + 1) : x;
                    float sum = 0.0f;
                    int n = 0;
                    for (int zz = z0; zz <= z1; ++zz) {
                        for (int yy = y0; yy <= y1; ++yy) {
                            for (int xx = x0; xx <= x1; ++xx) {
                                const size_t idx = static_cast<size_t>(zz) * nxSz * nySz
                                    + static_cast<size_t>(yy) * nxSz
                                    + static_cast<size_t>(xx);
                                sum += field[idx];
                                ++n;
                            }
                        }
                    }
                    const size_t outIdx = static_cast<size_t>(z) * nxSz * nySz
                        + static_cast<size_t>(y) * nxSz
                        + static_cast<size_t>(x);
                    scratch[outIdx] = sum / static_cast<float>(n);
                }
            }
        }
        field.swap(scratch);
    }
}

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
    int cubeSize = 512;
    int numThreads = 4;
    int sx = 2;
    int sy = 2;
    int sz = 1;
    if (argc > 1) numThreads = std::atoi(argv[1]);
    if (argc > 2) sx = std::atoi(argv[2]);
    if (argc > 3) sy = std::atoi(argv[3]);
    if (argc > 4) sz = std::atoi(argv[4]);
    if (argc > 5) cubeSize = std::atoi(argv[5]);
    if (numThreads < 1) numThreads = 1;
    if (sx < 1 || sy < 1 || sz < 1) {
        std::cerr << "split tuple must be positive" << std::endl;
        return 2;
    }
    if (cubeSize < 8) {
        std::cerr << "cubeSize must be >= 8" << std::endl;
        return 3;
    }
    const int nx = cubeSize;
    const int ny = cubeSize;
    const int nz = cubeSize;

    const int partitions = sx * sy * sz;
    std::cout << "[bench3d] config"
              << " dims=(" << nx << "," << ny << "," << nz << ")"
              << " threads=" << numThreads
              << " splits=(" << sx << "," << sy << "," << sz << ")"
              << " partitions=" << partitions << std::endl;

    const auto tFieldStart = std::chrono::steady_clock::now();
    std::vector<float> coarseNoise((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
    std::vector<float> detailNoise((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);
    std::vector<float> field((size_t)nx * (size_t)ny * (size_t)nz, 0.0f);

    std::mt19937 rng(123456789u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (size_t i = 0; i < field.size(); ++i) {
        coarseNoise[i] = dist(rng);
        detailNoise[i] = dist(rng);
    }
    smoothField3DInPlace(coarseNoise, nx, ny, nz, 30);
    smoothField3DInPlace(detailNoise, nx, ny, nz, 2);
    const float coarseScale = 100.0f;
    const float detailScale = 0.05f * coarseScale;
    for (size_t i = 0; i < field.size(); ++i) {
        field[i] = coarseScale * coarseNoise[i] + detailScale * detailNoise[i];
    }
    const auto tFieldEnd = std::chrono::steady_clock::now();
    const long long fieldGenMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(tFieldEnd - tFieldStart).count();

    const auto tInitStart = std::chrono::steady_clock::now();
    GridType grid(GInt::Vec3i(nx, ny, nz), GInt::Vec3b(false, false, false));
    GridFuncType gridfunc(&grid, field.data());
    MeshType mesh(&grid);
    MeshFuncType meshfunc;
    meshfunc.setMeshAndAllocate(&mesh);
    for (INDEX_TYPE cid = 0; cid < mesh.numCells(); ++cid) {
        GInt::Vec3l c;
        mesh.cellid2Coords(cid, c);
        const INDEX_TYPE vx = std::min(static_cast<INDEX_TYPE>(nx - 1), c[0] / 2);
        const INDEX_TYPE vy = std::min(static_cast<INDEX_TYPE>(ny - 1), c[1] / 2);
        const INDEX_TYPE vz = std::min(static_cast<INDEX_TYPE>(nz - 1), c[2] / 2);
        const size_t linear = static_cast<size_t>(vx)
            + static_cast<size_t>(nx) * static_cast<size_t>(vy)
            + static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(vz);
        meshfunc.setCellValue(cid, field[linear]);
    }
    GradType grad(&mesh);
    grad.ClearAllGradient();
    for (INDEX_TYPE cid = 0; cid < mesh.numCells(); ++cid) {
        grad.setAssigned(cid, 1);
        grad.setCritical(cid, true);
    }
    const auto tInitEnd = std::chrono::steady_clock::now();
    const long long initMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(tInitEnd - tInitStart).count();

    const float maxval = gridfunc.GetMaxValue();
    const float minval = gridfunc.GetMinValue();
    const float localPers = 0.0f;
    const float globalPers = 0.05f * (maxval - minval);

    const auto tSerialConstructStart = std::chrono::steady_clock::now();
    SerialMscType serial(&grad, &mesh, &meshfunc);
    serial.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
    serial.ComputeFromGrad();
    const auto tSerialConstructEnd = std::chrono::steady_clock::now();

    const auto tSerialSimplifyStart = std::chrono::steady_clock::now();
    serial.ComputeHierarchy(globalPers);
    serial.SetSelectPersAbs(globalPers);
    const auto tSerialSimplifyEnd = std::chrono::steady_clock::now();

    const long long serialConstructMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(tSerialConstructEnd - tSerialConstructStart).count();
    const long long serialSimplifyMs =
        std::chrono::duration_cast<std::chrono::milliseconds>(tSerialSimplifyEnd - tSerialSimplifyStart).count();

    PartitionedPipelineType partitioned(&grad, &mesh, &meshfunc);
    PartitionedPipelineType::TimingBreakdown partitionTimings;
    const GInt::Vec3l splitTuple(sx, sy, sz);
    std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
        partitioned.BuildPartitionLocalMSCs(splitTuple, localPers, &partitionTimings, numThreads);
    GInt::PartitionedTopologicalRegularGrid3D partitionGrid(&mesh, splitTuple);
    std::unique_ptr<PartitionedPipelineType::ReconciledGlobalMsc> reconciled =
        partitioned.BuildReconciledGlobalBase(partitionGrid, localResults, &partitionTimings);
    const auto tGlobalSimplifyStart = std::chrono::steady_clock::now();
    reconciled->ComputeHierarchy(globalPers);
    reconciled->SetSelectPersAbs(globalPers);
    const auto tGlobalSimplifyEnd = std::chrono::steady_clock::now();
    partitionTimings.global_simplify_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(tGlobalSimplifyEnd - tGlobalSimplifyStart).count();

    long long preFrozenSum = 0;
    long long postFrozenSum = 0;
    for (size_t i = 0; i < localResults.size(); ++i) {
        preFrozenSum += static_cast<long long>(localResults[i].pre_exchange_frozen_nodes);
        postFrozenSum += static_cast<long long>(localResults[i].post_exchange_frozen_nodes);
    }

    const size_t serialNodes = countLivingNodes(&serial);
    const size_t serialArcs = countLivingArcs(&serial);
    const size_t reconciledNodes = countLivingNodes(reconciled.get());
    const size_t reconciledArcs = countLivingArcs(reconciled.get());
    const auto serialHist = livingEndpointHistogram(&serial);
    const auto reconciledHist = livingEndpointHistogram(reconciled.get());
    const long long partitionedTotalMs =
        partitionTimings.local_stage_total_ms + partitionTimings.reconcile_ms + partitionTimings.global_simplify_ms;

    std::cout << "benchmark3d_ms"
              << " threads=" << numThreads
              << " partitions=" << partitionGrid.num_partitions()
              << " splits=(" << sx << "," << sy << "," << sz << ")"
              << " field_gen=" << fieldGenMs
              << " init=" << initMs
              << " serial_construct=" << serialConstructMs
              << " serial_simplify=" << serialSimplifyMs
              << " partition_local_total=" << partitionTimings.local_stage_total_ms
              << " partition_local_construct=" << partitionTimings.local_build_ms
              << " partition_local_simplify=" << partitionTimings.local_simplify_ms
              << " partition_reconcile=" << partitionTimings.reconcile_ms
              << " partition_global_simplify=" << partitionTimings.global_simplify_ms
              << " partition_total=" << partitionedTotalMs
              << std::endl;

    std::cout << "benchmark3d_parity"
              << " threads=" << numThreads
              << " partitions=" << partitionGrid.num_partitions()
              << " splits=(" << sx << "," << sy << "," << sz << ")"
              << " local_pre_exchange_frozen_nodes_sum=" << preFrozenSum
              << " local_post_exchange_frozen_nodes_sum=" << postFrozenSum
              << " serial_nodes=" << serialNodes
              << " serial_arcs=" << serialArcs
              << " reconciled_nodes=" << reconciledNodes
              << " reconciled_arcs=" << reconciledArcs
              << " same_endpoint_histogram=" << (serialHist == reconciledHist ? "true" : "false")
              << std::endl;

    return 0;
}
