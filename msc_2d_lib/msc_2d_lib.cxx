#include "msc_2d_lib.h"

#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "gi_discrete_gradient_computer.h"
#include "gi_graphs.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_ms_complex_to_graph.h"

namespace GInt {
namespace Msc2D {

typedef GInt::MorseSmaleComplexBasic<float, Accurate2D::MeshType, Accurate2D::MeshFuncType, Accurate2D::GradType> MyMscType;
typedef GInt::MorseSmaleComplexPartitioned<float, Accurate2D::MeshType, Accurate2D::MeshFuncType, Accurate2D::GradType> PartitionedPipelineType;
typedef PartitionedPipelineType::ReconciledGlobalMsc ReconciledMscType;

namespace {
int clampSupportedPartitionCount(int requested) {
    static const int kSupportedCounts[] = { 1, 2, 3, 4, 6, 8, 9, 12, 16 };
    if (requested <= kSupportedCounts[0]) {
        return kSupportedCounts[0];
    }
    int best = kSupportedCounts[0];
    for (size_t i = 0; i < sizeof(kSupportedCounts) / sizeof(kSupportedCounts[0]); ++i) {
        if (kSupportedCounts[i] > requested) break;
        best = kSupportedCounts[i];
    }
    return best;
}
}

struct Msc2D::Impl {
    std::unique_ptr<Accurate2D::DiscreteGradientBuilder> dgb;
    Accurate2D::GridType* grid;
    Accurate2D::GridFuncType* gridfunc;
    Accurate2D::MeshType* mesh;
    Accurate2D::MeshFuncType* meshfunc;
    Accurate2D::GradType* grad;
    std::unique_ptr<MyMscType> serialMsc;
    std::unique_ptr<ReconciledMscType> partitionedMsc;
    MyMscType* activeMsc;
    std::unique_ptr<GInt::Geometric2DGraph> geomLineGraph;
    int mX;
    int mY;
    std::vector<float> rawData;
    std::vector<int> baseLabelingAsc2;
    std::vector<int> baseLabelingDsc2;
    float selectedPersistence;
    float basePersistence;
    int effectiveParallelismValue;
    Msc2D::BuilderMode builderMode;
    bool hasCompute;

    Impl()
        : dgb(new Accurate2D::DiscreteGradientBuilder()),
          grid(NULL),
          gridfunc(NULL),
          mesh(NULL),
          meshfunc(NULL),
          grad(NULL),
          activeMsc(NULL),
          mX(-1),
          mY(-1),
          selectedPersistence(0.0f),
          basePersistence(0.0f),
          effectiveParallelismValue(1),
          builderMode(Msc2D::BuilderMode::Serial),
          hasCompute(false) {}

    void resetComputedState() {
        serialMsc.reset();
        partitionedMsc.reset();
        activeMsc = NULL;
        geomLineGraph.reset();
        rawData.clear();
        baseLabelingAsc2.clear();
        baseLabelingDsc2.clear();
        grid = NULL;
        gridfunc = NULL;
        mesh = NULL;
        meshfunc = NULL;
        grad = NULL;
        mX = -1;
        mY = -1;
        selectedPersistence = 0.0f;
        basePersistence = 0.0f;
        effectiveParallelismValue = 1;
        builderMode = Msc2D::BuilderMode::Serial;
        hasCompute = false;
    }

    MyMscType* activeMscOrThrow() const {
        if (!activeMsc) {
            throw std::runtime_error("MSC result is not available. Call compute() first.");
        }
        return activeMsc;
    }

    void ensureComputed() const {
        if (!hasCompute || !activeMsc) {
            throw std::runtime_error("MSC result is not available. Call compute() first.");
        }
    }
};

Msc2D::Msc2D() : m_impl(new Impl()) {}

Msc2D::~Msc2D() {}

Msc2D::Msc2D(Msc2D&& other) noexcept : m_impl(std::move(other.m_impl)) {}

Msc2D& Msc2D::operator=(Msc2D&& other) noexcept {
    if (this != &other) {
        m_impl = std::move(other.m_impl);
    }
    return *this;
}

void Msc2D::compute(const float* rowMajorValues, int rows, int cols, bool accurateAsc, bool accurateDsc) {
    ComputeOptions options;
    options.accurateAsc = accurateAsc;
    options.accurateDsc = accurateDsc;
    compute(rowMajorValues, rows, cols, options);
}

void Msc2D::compute(const float* rowMajorValues, int rows, int cols, const ComputeOptions& options) {
    if (rowMajorValues == NULL) {
        throw std::invalid_argument("compute() received a null input buffer.");
    }
    if (rows <= 0 || cols <= 0) {
        throw std::invalid_argument("compute() requires positive rows and cols.");
    }

    m_impl->resetComputedState();

    const int pX = rows;
    const int pY = cols;
    const int mX = pY;
    const int mY = pX;

    m_impl->mX = mX;
    m_impl->mY = mY;
    m_impl->rawData.resize(static_cast<size_t>(pX) * static_cast<size_t>(pY));

    for (int i = 0; i < pX; ++i) {
        for (int j = 0; j < pY; ++j) {
            m_impl->rawData[static_cast<size_t>(j + i * mX)] = rowMajorValues[static_cast<size_t>(i * pY + j)];
        }
    }

    const int effectiveParallelism = clampSupportedPartitionCount(options.requestedParallelism);
    m_impl->effectiveParallelismValue = effectiveParallelism;
    m_impl->builderMode = options.builderMode;

    m_impl->dgb->SetFloadArrayAndDims(mX, mY, m_impl->rawData.data());
    m_impl->dgb->SetNeededAccuracy(options.accurateAsc, options.accurateDsc);
    m_impl->dgb->SetParallelism(effectiveParallelism);
    m_impl->dgb->ComputeDiscreteGradient();

    m_impl->grid = m_impl->dgb->GetGrid();
    m_impl->gridfunc = m_impl->dgb->GetGridFunc();
    m_impl->mesh = m_impl->dgb->GetTopoMesh();
    m_impl->meshfunc = m_impl->dgb->GetMeshFunc();
    m_impl->grad = m_impl->dgb->GetGrad();

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

    m_impl->basePersistence = basePersistence;
    m_impl->selectedPersistence = cancelPersistence;

    if (options.builderMode == BuilderMode::Partitioned) {
        PartitionedPipelineType partitioned(m_impl->grad, m_impl->mesh, m_impl->meshfunc);
        std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
            partitioned.BuildPartitionLocalMSCs(effectiveParallelism, basePersistence, NULL);
        GInt::PartitionedTopologicalRegularGrid2D partitionMesh(m_impl->mesh, effectiveParallelism);
        m_impl->partitionedMsc = partitioned.BuildReconciledGlobalBase(partitionMesh, localResults, NULL);
        m_impl->partitionedMsc->ComputeHierarchy(cancelPersistence);
        m_impl->partitionedMsc->SetSelectPersMAX();
        m_impl->activeMsc = m_impl->partitionedMsc.get();
    } else {
        m_impl->serialMsc.reset(new MyMscType(m_impl->grad, m_impl->mesh, m_impl->meshfunc));
        m_impl->serialMsc->SetBuildArcGeometry(Vec3b(false, false, false));
        m_impl->serialMsc->ComputeFromGrad();
        m_impl->serialMsc->ComputeHierarchy(cancelPersistence);
        m_impl->serialMsc->SetSelectPersMAX();
        m_impl->activeMsc = m_impl->serialMsc.get();
    }

    m_impl->hasCompute = true;
}

void Msc2D::setPersistence(float value) {
    m_impl->ensureComputed();
    m_impl->selectedPersistence = value;
    m_impl->activeMscOrThrow()->SetSelectPersAbs(value);
}

LabelImage Msc2D::ascending2Manifolds() {
    m_impl->ensureComputed();

    if (m_impl->baseLabelingAsc2.empty()) {
        m_impl->baseLabelingAsc2.assign(m_impl->grid->NumElements(), -1);
        MyMscType* activeMsc = m_impl->activeMscOrThrow();
        activeMsc->SetSelectPersAbs(-1);

        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 0) continue;

            std::set<INDEX_TYPE> manifold;
            activeMsc->fillGeometry(nid, manifold, true);
            for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                if (m_impl->mesh->dimension(*it) != 0) continue;
                m_impl->baseLabelingAsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = static_cast<int>(nid);
            }
        }
    }

    MyMscType* activeMsc = m_impl->activeMscOrThrow();
    activeMsc->SetSelectPersAbs(m_impl->selectedPersistence);

    std::unordered_map<INT_TYPE, int> remap;
    MyMscType::LivingNodesIterator nit(activeMsc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (activeMsc->getNode(nid).dim != 0) continue;

        std::set<INT_TYPE> constituents;
        activeMsc->GatherNodes(nid, constituents, true);
        for (std::set<INT_TYPE>::const_iterator it = constituents.begin(); it != constituents.end(); ++it) {
            remap[*it] = static_cast<int>(nid);
        }
    }

    LabelImage out;
    out.width = m_impl->mX;
    out.height = m_impl->mY;
    out.labels.assign(static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY), -1);

    for (int i = 0; i < m_impl->mX * m_impl->mY; ++i) {
        const int base = m_impl->baseLabelingAsc2[i];
        const std::unordered_map<INT_TYPE, int>::const_iterator it = remap.find(base);
        if (it != remap.end()) {
            out.labels[static_cast<size_t>(i)] = it->second;
        }
    }

    return out;
}

LabelImage Msc2D::descending2Manifolds() {
    m_impl->ensureComputed();

    if (m_impl->baseLabelingDsc2.empty()) {
        m_impl->baseLabelingDsc2.assign(m_impl->grid->NumElements(), -1);
        MyMscType* activeMsc = m_impl->activeMscOrThrow();
        activeMsc->SetSelectPersAbs(-1);

        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 2) continue;

            std::set<INDEX_TYPE> manifold;
            activeMsc->fillGeometry(nid, manifold, false);
            for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                if (m_impl->mesh->dimension(*it) != 2) continue;
                m_impl->baseLabelingDsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = static_cast<int>(nid);
            }
        }
    }

    MyMscType* activeMsc = m_impl->activeMscOrThrow();
    activeMsc->SetSelectPersAbs(m_impl->selectedPersistence);

    std::unordered_map<INT_TYPE, int> remap;
    MyMscType::LivingNodesIterator nit(activeMsc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (activeMsc->getNode(nid).dim != 2) continue;

        std::set<INT_TYPE> constituents;
        activeMsc->GatherNodes(nid, constituents, false);
        for (std::set<INT_TYPE>::const_iterator it = constituents.begin(); it != constituents.end(); ++it) {
            remap[*it] = static_cast<int>(nid);
        }
    }

    LabelImage out;
    out.width = m_impl->mX;
    out.height = m_impl->mY;
    out.labels.assign(static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY), -1);

    for (int i = 0; i < m_impl->mX * m_impl->mY; ++i) {
        const int base = m_impl->baseLabelingDsc2[i];
        const std::unordered_map<INT_TYPE, int>::const_iterator it = remap.find(base);
        if (it != remap.end()) {
            out.labels[static_cast<size_t>(i)] = it->second;
        }
    }

    return out;
}

std::vector<CriticalPoint> Msc2D::criticalPoints() const {
    m_impl->ensureComputed();
    MyMscType* activeMsc = m_impl->activeMscOrThrow();

    std::set<INT_TYPE> livingNodeIds;
    MyMscType::LivingNodesIterator nit(activeMsc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        livingNodeIds.insert(nit.value());
    }

    std::vector<CriticalPoint> output;
    output.reserve(livingNodeIds.size());

    for (std::set<INT_TYPE>::const_iterator it = livingNodeIds.begin(); it != livingNodeIds.end(); ++it) {
        const INT_TYPE id = *it;
        const node<float>& nodeRef = activeMsc->getNode(id);
        GInt::Vec2l coords;
        m_impl->mesh->cellid2Coords(nodeRef.cellindex, coords);
        GInt::Vec2f fcoords = coords;
        fcoords *= 0.5f;

        CriticalPoint cp;
        cp.id = static_cast<int>(id);
        cp.x = fcoords[0];
        cp.y = fcoords[1];
        cp.dim = nodeRef.dim;
        cp.value = nodeRef.value;
        output.push_back(cp);
    }

    return output;
}

void Msc2D::computePolylineGraph(bool useValleys) {
    m_impl->ensureComputed();
    MyMscType* activeMsc = m_impl->activeMscOrThrow();

    MeshCellsGraph* graph = NULL;
    if (useValleys) {
        graph = GInt::BuildMeshCellsGraphFromMSCValleys<MyMscType, Accurate2D::MeshType>(activeMsc, m_impl->mesh);
    } else {
        graph = GInt::BuildMeshCellsGraphFromMSCRidges<MyMscType, Accurate2D::MeshType>(activeMsc, m_impl->mesh);
    }

    m_impl->geomLineGraph.reset(GInt::BuildGeometricGraphFromMeshGraph<Accurate2D::MeshType>(graph, m_impl->mesh, 10));
    delete graph;
}

Graph Msc2D::graph() const {
    Graph out;
    if (!m_impl->geomLineGraph) {
        return out;
    }

    GInt::Geometric2DGraph::vertex_iterator vit(m_impl->geomLineGraph.get());
    for (vit.begin(); vit.valid(); vit.advance()) {
        const auto vid = vit.value();
        const auto& v = m_impl->geomLineGraph->GetVertex(vid);

        Node nodeEntry;
        nodeEntry.id = static_cast<int>(v.vid);
        nodeEntry.edges.reserve(v.edges.size());
        for (size_t i = 0; i < v.edges.size(); ++i) {
            nodeEntry.edges.push_back(static_cast<int>(v.edges[i]));
        }
        nodeEntry.geometry.push_back(Point{ v.store[0], v.store[1] });
        out.nodes.push_back(nodeEntry);
    }

    GInt::Geometric2DGraph::edge_iterator eit(m_impl->geomLineGraph.get());
    for (eit.begin(); eit.valid(); eit.advance()) {
        const auto eid = eit.value();
        const auto& ge = m_impl->geomLineGraph->GetEdge(eid);

        Edge edgeEntry;
        edgeEntry.id = static_cast<int>(ge.eid);
        edgeEntry.from = static_cast<int>(ge.v1);
        edgeEntry.to = static_cast<int>(ge.v2);
        const std::vector<GInt::Vec2f>& line = ge.store->GetLine();
        edgeEntry.geometry.reserve(line.size());
        for (size_t i = 0; i < line.size(); ++i) {
            edgeEntry.geometry.push_back(Point{ line[i][0], line[i][1] });
        }
        out.edges.push_back(edgeEntry);
    }

    return out;
}

bool Msc2D::hasResult() const {
    return m_impl->hasCompute;
}

int Msc2D::width() const {
    return m_impl->mX;
}

int Msc2D::height() const {
    return m_impl->mY;
}

int Msc2D::effectiveParallelism() const {
    return m_impl->effectiveParallelismValue;
}

} // namespace Msc2D
} // namespace GInt
