#include "msc_3d_lib.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "gi_discrete_gradient_labeling.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_regular_grid.h"

namespace GInt {
namespace Msc3D {

typedef GInt::RegularGrid3D GridType;
typedef GInt::RegularGridTrilinearFunction GridFuncType;
typedef GInt::TopologicalRegularGrid3D MeshType;
typedef GInt::DiscreteGradientLabeling<MeshType> GradType;
typedef GInt::TopologicalExplicitDenseMeshFunction<MeshType, float> MeshFuncType;
typedef GInt::MorseSmaleComplexBasic<float, MeshType, MeshFuncType, GradType> MyMscType;
typedef GInt::MorseSmaleComplexPartitioned<float, MeshType, MeshFuncType, GradType> PartitionedPipelineType;
typedef PartitionedPipelineType::ReconciledGlobalMsc ReconciledMscType;

struct Msc3D::Impl {
    std::vector<float> rawData;
    std::unique_ptr<GridType> grid;
    std::unique_ptr<GridFuncType> gridfunc;
    std::unique_ptr<MeshType> mesh;
    std::unique_ptr<MeshFuncType> meshfunc;
    std::unique_ptr<GradType> grad;
    std::unique_ptr<MyMscType> serialMsc;
    std::unique_ptr<ReconciledMscType> partitionedMsc;
    MyMscType* activeMsc;
    float selectedPersistence;
    int dimX;
    int dimY;
    int dimZ;
    int effectivePartitionCountValue;
    Point effectivePartitionSplitsValue;
    int effectiveThreadCountValue;
    bool hasCompute;

    Impl()
        : activeMsc(NULL),
          selectedPersistence(0.0f),
          dimX(-1),
          dimY(-1),
          dimZ(-1),
          effectivePartitionCountValue(1),
          effectiveThreadCountValue(1),
          hasCompute(false) {
        effectivePartitionSplitsValue.x = 1.0f;
        effectivePartitionSplitsValue.y = 1.0f;
        effectivePartitionSplitsValue.z = 1.0f;
    }

    void resetComputedState() {
        rawData.clear();
        grid.reset();
        gridfunc.reset();
        mesh.reset();
        meshfunc.reset();
        grad.reset();
        serialMsc.reset();
        partitionedMsc.reset();
        activeMsc = NULL;
        selectedPersistence = 0.0f;
        dimX = -1;
        dimY = -1;
        dimZ = -1;
        effectivePartitionCountValue = 1;
        effectivePartitionSplitsValue.x = 1.0f;
        effectivePartitionSplitsValue.y = 1.0f;
        effectivePartitionSplitsValue.z = 1.0f;
        effectiveThreadCountValue = 1;
        hasCompute = false;
    }

    MyMscType* activeMscOrThrow() const {
        if (activeMsc == NULL) {
            throw std::runtime_error("MSC result is not available. Call compute() first.");
        }
        return activeMsc;
    }
};

Msc3D::Msc3D() : m_impl(new Impl()) {}

Msc3D::~Msc3D() {}

Msc3D::Msc3D(Msc3D&& other) noexcept : m_impl(std::move(other.m_impl)) {}

Msc3D& Msc3D::operator=(Msc3D&& other) noexcept {
    if (this != &other) {
        m_impl = std::move(other.m_impl);
    }
    return *this;
}

void Msc3D::compute(const float* rowMajorValues, int xdim, int ydim, int zdim, const ComputeOptions& options) {
    if (rowMajorValues == NULL) {
        throw std::invalid_argument("compute() received a null input buffer.");
    }
    if (xdim <= 0 || ydim <= 0 || zdim <= 0) {
        throw std::invalid_argument("compute() requires positive xdim/ydim/zdim.");
    }

    m_impl->resetComputedState();
    m_impl->dimX = xdim;
    m_impl->dimY = ydim;
    m_impl->dimZ = zdim;
    m_impl->rawData.assign(rowMajorValues, rowMajorValues + (size_t)xdim * (size_t)ydim * (size_t)zdim);

    m_impl->grid.reset(new GridType(Vec3i(xdim, ydim, zdim), Vec3b(false, false, false)));
    m_impl->gridfunc.reset(new GridFuncType(m_impl->grid.get(), m_impl->rawData.data()));
    m_impl->mesh.reset(new MeshType(m_impl->grid.get()));
    m_impl->meshfunc.reset(new MeshFuncType());
    m_impl->meshfunc->setMeshAndAllocate(m_impl->mesh.get());
    for (INDEX_TYPE cid = 0; cid < m_impl->mesh->numCells(); ++cid) {
        Vec3l c;
        m_impl->mesh->cellid2Coords(cid, c);
        const INDEX_TYPE vx = std::min(static_cast<INDEX_TYPE>(xdim - 1), c[0] / 2);
        const INDEX_TYPE vy = std::min(static_cast<INDEX_TYPE>(ydim - 1), c[1] / 2);
        const INDEX_TYPE vz = std::min(static_cast<INDEX_TYPE>(zdim - 1), c[2] / 2);
        const size_t linear = static_cast<size_t>(vx)
            + static_cast<size_t>(xdim) * static_cast<size_t>(vy)
            + static_cast<size_t>(xdim) * static_cast<size_t>(ydim) * static_cast<size_t>(vz);
        m_impl->meshfunc->setCellValue(cid, m_impl->rawData[linear]);
    }
    m_impl->grad.reset(new GradType(m_impl->mesh.get()));
    m_impl->grad->ClearAllGradient();
    for (INDEX_TYPE cid = 0; cid < m_impl->mesh->numCells(); ++cid) {
        m_impl->grad->setAssigned(cid, 1);
        m_impl->grad->setCritical(cid, true);
    }

    const float maxval = m_impl->gridfunc->GetMaxValue();
    const float minval = m_impl->gridfunc->GetMinValue();
    const float valueRange = maxval - minval;
    const float defaultBasePersistence = 0.01f * valueRange;
    const float defaultCancelPersistence = 0.1f * valueRange;
    const float cancelPersistence = (options.cancelPersistenceAbs >= 0.0f) ? options.cancelPersistenceAbs : defaultCancelPersistence;
    float basePersistence = (options.basePersistenceAbs >= 0.0f) ? options.basePersistenceAbs : defaultBasePersistence;
    if (basePersistence > cancelPersistence) {
        basePersistence = cancelPersistence;
    }

    m_impl->selectedPersistence = cancelPersistence;
    const int requestedThreads = (options.numThreads < 1) ? 1 : options.numThreads;
    m_impl->effectiveThreadCountValue = requestedThreads;

    if (options.builderMode == BuilderMode::Partitioned) {
        PartitionedPipelineType partitioned(m_impl->grad.get(), m_impl->mesh.get(), m_impl->meshfunc.get());
        std::vector<PartitionedPipelineType::PartitionRunResult> localResults;
        PartitionedTopologicalRegularGrid3D partitionGrid(m_impl->mesh.get(), 1);
        if (options.useExplicitSplits) {
            const Vec3l splits(
                static_cast<INDEX_TYPE>(options.splitX),
                static_cast<INDEX_TYPE>(options.splitY),
                static_cast<INDEX_TYPE>(options.splitZ));
            localResults = partitioned.BuildPartitionLocalMSCs(splits, basePersistence, NULL, requestedThreads);
            partitionGrid = PartitionedTopologicalRegularGrid3D(m_impl->mesh.get(), splits);
        } else {
            const int requestedPartitions = (options.requestedPartitions < 1) ? 1 : options.requestedPartitions;
            localResults = partitioned.BuildPartitionLocalMSCs(requestedPartitions, basePersistence, NULL, requestedThreads);
            partitionGrid = PartitionedTopologicalRegularGrid3D(m_impl->mesh.get(), requestedPartitions);
        }
        m_impl->effectivePartitionCountValue = partitionGrid.num_partitions();
        const Vec3l effSplits = partitionGrid.partition_splits();
        m_impl->effectivePartitionSplitsValue.x = static_cast<float>(effSplits[0]);
        m_impl->effectivePartitionSplitsValue.y = static_cast<float>(effSplits[1]);
        m_impl->effectivePartitionSplitsValue.z = static_cast<float>(effSplits[2]);
        m_impl->partitionedMsc = partitioned.BuildReconciledGlobalBase(partitionGrid, localResults, NULL);
        m_impl->partitionedMsc->ComputeHierarchy(cancelPersistence);
        m_impl->partitionedMsc->SetSelectPersMAX();
        m_impl->activeMsc = m_impl->partitionedMsc.get();
    } else {
        m_impl->serialMsc.reset(new MyMscType(m_impl->grad.get(), m_impl->mesh.get(), m_impl->meshfunc.get()));
        m_impl->serialMsc->SetBuildArcGeometry(Vec3b(false, false, false));
        m_impl->serialMsc->ComputeFromGrad();
        m_impl->serialMsc->ComputeHierarchy(cancelPersistence);
        m_impl->serialMsc->SetSelectPersMAX();
        m_impl->activeMsc = m_impl->serialMsc.get();
        m_impl->effectivePartitionCountValue = 1;
        m_impl->effectivePartitionSplitsValue.x = 1.0f;
        m_impl->effectivePartitionSplitsValue.y = 1.0f;
        m_impl->effectivePartitionSplitsValue.z = 1.0f;
    }

    m_impl->hasCompute = true;
}

void Msc3D::setPersistence(float value) {
    m_impl->activeMscOrThrow()->SetSelectPersAbs(value);
    m_impl->selectedPersistence = value;
}

std::vector<CriticalPoint> Msc3D::criticalPoints() const {
    const MyMscType* active = m_impl->activeMscOrThrow();
    std::vector<CriticalPoint> out;
    MyMscType::LivingNodesIterator nit(const_cast<MyMscType*>(active));
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        const node<float>& n = active->getNode(nid);
        Vec3l coords;
        m_impl->mesh->cellid2Coords(n.cellindex, coords);
        CriticalPoint cp;
        cp.id = static_cast<int>(nid);
        cp.x = static_cast<float>(coords[0]) * 0.5f;
        cp.y = static_cast<float>(coords[1]) * 0.5f;
        cp.z = static_cast<float>(coords[2]) * 0.5f;
        cp.dim = static_cast<int>(n.dim);
        cp.value = n.value;
        out.push_back(cp);
    }
    return out;
}

bool Msc3D::hasResult() const {
    return m_impl->hasCompute && m_impl->activeMsc != NULL;
}

int Msc3D::xdim() const { return m_impl->dimX; }
int Msc3D::ydim() const { return m_impl->dimY; }
int Msc3D::zdim() const { return m_impl->dimZ; }
int Msc3D::effectivePartitionCount() const { return m_impl->effectivePartitionCountValue; }
Point Msc3D::effectivePartitionSplits() const { return m_impl->effectivePartitionSplitsValue; }
int Msc3D::effectiveThreadCount() const { return m_impl->effectiveThreadCountValue; }

} // namespace Msc3D
} // namespace GInt
