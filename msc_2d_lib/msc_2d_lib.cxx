#include "msc_2d_lib.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "gi_discrete_gradient_computer.h"
#include "gi_graphs.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_partitioned.h"
#include "gi_ms_complex_to_graph.h"

#ifdef MSC2D_HAS_GPU_DGRAD
#include "dgrad_gpu_api.h"
#include "dgrad2d_tables.h"
#endif

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

bool shouldEmitLabelDiagnostics() {
    const char* env = std::getenv("MSC2D_LABEL_DIAGNOSTICS");
    return env != NULL && env[0] != '\0' && env[0] != '0';
}

long long msSince(const std::chrono::steady_clock::time_point& t0) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - t0)
        .count();
}

// Assign (or look up) a dense 0..M-1 id for a raw base identity. The raw
// identities stored per pixel live in a huge range in partitioned mode (they are
// topological cell indices, up to grid->NumElements()), but there are only
// thousands of distinct ones. Compacting them once at base-labeling time turns
// the per-pixel remap lookup into an array index instead of a hash probe --
// without paying for a table sized by the raw range.
int compactIdFor(std::unordered_map<int, int>& toCompact, int raw) {
    const std::unordered_map<int, int>::const_iterator it = toCompact.find(raw);
    if (it != toCompact.end()) return it->second;
    const int compact = static_cast<int>(toCompact.size());
    toCompact.insert(std::make_pair(raw, compact));
    return compact;
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
    // Per-pixel base identity, stored as a compact 0..M-1 id (-1 = unlabeled).
    // The raw identity (node id in serial mode, base cell index in partitioned
    // mode) is recoverable from the matching baseCompact* map.
    std::vector<int> baseLabelingAsc2;
    std::vector<int> baseLabelingDsc2;
    // Raw base identity -> compact id, filled once alongside the base labeling.
    std::unordered_map<int, int> baseCompactAsc2;
    std::unordered_map<int, int> baseCompactDsc2;
    // Per-call scratch: compact base id -> living node id (-1 = no living owner).
    // Kept as members so the per-call allocation is amortized away.
    std::vector<int> remapDenseAsc2;
    std::vector<int> remapDenseDsc2;
    float selectedPersistence;
    float basePersistence;
    int effectiveParallelismValue;
    Msc2D::BuilderMode builderMode;
    bool hasCompute;
    bool builtArcGeometry;   // whether per-arc geometry was constructed at compute
    // Whether compute() ran with ComputeOptions::useGpuGradient. The base
    // manifold labeling reads it to use the GPU flow-terminal path.
    bool useGpuGradient;

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
          hasCompute(false),
          builtArcGeometry(false),
          useGpuGradient(false) {}

    void resetComputedState() {
        serialMsc.reset();
        partitionedMsc.reset();
        activeMsc = NULL;
        geomLineGraph.reset();
        rawData.clear();
        baseLabelingAsc2.clear();
        baseLabelingDsc2.clear();
        baseCompactAsc2.clear();
        baseCompactDsc2.clear();
        remapDenseAsc2.clear();
        remapDenseDsc2.clear();
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
        builtArcGeometry = false;
        useGpuGradient = false;
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
    m_impl->useGpuGradient = options.useGpuGradient;
#ifdef MSC2D_HAS_GPU_DGRAD
    if (options.useGpuGradient) {
        m_impl->dgb->SetGradientOverride(
            [](Accurate2D::GridType* grid, Accurate2D::GridFuncType* func,
               Accurate2D::MeshType* mesh, Accurate2D::GradType* grad) -> bool {
                gpu::Dgrad2DTables tables;
                if (!gpu::BuildDgrad2DTablesFromMesh(*mesh, tables)) return false;
                const auto xy = grid->XY();
                return gpu::ComputeDiscreteGradient2DFused(
                    func->GetImage(), xy[0], xy[1], tables,
                    reinterpret_cast<uint8_t*>(grad->m_dgrad->LabelArray()),
                    NULL);
            });
    } else {
        m_impl->dgb->SetGradientOverride(NULL);
    }
#else
    if (options.useGpuGradient) {
        printf("Msc2D: built without GPU gradient support (GPU_DGRAD_ENABLED off); using CPU path\n");
    }
#endif
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

    const Vec3b arcGeometryFlags(
        options.buildArcGeometry[0], options.buildArcGeometry[1], options.buildArcGeometry[2]);
    m_impl->builtArcGeometry =
        options.buildArcGeometry[0] || options.buildArcGeometry[1] || options.buildArcGeometry[2];

    if (options.builderMode == BuilderMode::Partitioned) {
        PartitionedPipelineType partitioned(m_impl->grad, m_impl->mesh, m_impl->meshfunc);
        partitioned.SetBuildArcGeometry(arcGeometryFlags);
        std::vector<PartitionedPipelineType::PartitionRunResult> localResults =
            partitioned.BuildPartitionLocalMSCs(effectiveParallelism, basePersistence, NULL);
        GInt::PartitionedTopologicalRegularGrid2D partitionMesh(m_impl->mesh, effectiveParallelism);
        m_impl->partitionedMsc = partitioned.BuildReconciledGlobalBase(partitionMesh, localResults, NULL);
        m_impl->partitionedMsc->ComputeHierarchy(cancelPersistence);
        m_impl->partitionedMsc->SetSelectPersMAX();
        m_impl->activeMsc = m_impl->partitionedMsc.get();
    } else {
        m_impl->serialMsc.reset(new MyMscType(m_impl->grad, m_impl->mesh, m_impl->meshfunc));
        m_impl->serialMsc->SetBuildArcGeometry(arcGeometryFlags);
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

#ifdef MSC2D_HAS_GPU_DGRAD
namespace {
// GPU base-manifold labeling via flow terminals (Label2DFlowTerminals):
// per-vertex basin-of-minimum (ascending) or per-quad region-of-maximum
// (descending), computed by offset-doubled flow on the gradient bytes and
// validated cell-for-cell against the serial fillGeometry paint. Compact ids
// are assigned in LivingNodesIterator order, keyed by NODE id -- exactly what
// the CPU serial paint does -- so the living remap downstream is unchanged.
// Returns false (nothing written) to fall back to the CPU paint.
bool gpuFillBaseLabeling(Accurate2D::MeshType* mesh, Accurate2D::GradType* grad,
                         int mX, int mY, MyMscType* activeMsc, bool ascending,
                         std::vector<int>& baseLabeling,
                         std::unordered_map<int, int>& baseCompact) {
    gpu::Dgrad2DTables tables;
    if (!gpu::BuildDgrad2DTablesFromMesh(*mesh, tables)) return false;
    const uint8_t* gradBytes =
        reinterpret_cast<const uint8_t*>(grad->m_dgrad->LabelArray());
    const int64_t X = mX;
    const int64_t Y = mY;
    if (ascending) {
        std::vector<int32_t> term(static_cast<size_t>(X) * static_cast<size_t>(Y));
        if (!gpu::Label2DFlowTerminals(gradBytes, X, Y, tables, term.data(), NULL, NULL))
            return false;
        // Terminal vertex number -> compact id of that minimum's node.
        std::vector<int> termCompact(term.size(), -1);
        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 0) continue;
            const INDEX_TYPE cell = activeMsc->getNode(nid).cellindex;
            const INDEX_TYPE vnum = mesh->VertexNumberFromCellID(cell);
            termCompact[static_cast<size_t>(vnum)] =
                compactIdFor(baseCompact, static_cast<int>(nid));
        }
        for (size_t i = 0; i < term.size(); ++i) {
            const int32_t t = term[i];
            if (t >= 0) baseLabeling[i] = termCompact[static_cast<size_t>(t)];
        }
    } else {
        const int64_t QX = X - 1;
        const int64_t QY = Y - 1;
        if (QX <= 0 || QY <= 0) return false;
        std::vector<int32_t> qterm(static_cast<size_t>(QX) * static_cast<size_t>(QY));
        if (!gpu::Label2DFlowTerminals(gradBytes, X, Y, tables, NULL, qterm.data(), NULL))
            return false;
        // Terminal quad lattice index -> compact id of that maximum's node.
        // A quad cell has odd coordinates (2qx+1, 2qy+1) on the (2X-1)-wide
        // cell lattice; its quad lattice index is qx + qy*QX.
        std::vector<int> termCompact(qterm.size(), -1);
        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 2) continue;
            const INDEX_TYPE cell = activeMsc->getNode(nid).cellindex;
            const int64_t cx = static_cast<int64_t>(cell) % (2 * X - 1);
            const int64_t cy = static_cast<int64_t>(cell) / (2 * X - 1);
            termCompact[static_cast<size_t>((cx >> 1) + (cy >> 1) * QX)] =
                compactIdFor(baseCompact, static_cast<int>(nid));
        }
        // Project each quad's label onto the same vertex the CPU paint stamps
        // (VertexNumberFromCellID of the quad cell) -- the mapping is
        // injective, so stamping order cannot matter.
        for (int64_t qy = 0; qy < QY; ++qy) {
            for (int64_t qx = 0; qx < QX; ++qx) {
                const int32_t t = qterm[static_cast<size_t>(qx + qy * QX)];
                if (t < 0) continue;
                const int compact = termCompact[static_cast<size_t>(t)];
                if (compact < 0) continue;
                const INDEX_TYPE cellid = (2 * qx + 1) + (2 * qy + 1) * (2 * X - 1);
                baseLabeling[static_cast<size_t>(mesh->VertexNumberFromCellID(cellid))] =
                    compact;
            }
        }
    }
    return true;
}
}  // namespace
#endif  // MSC2D_HAS_GPU_DGRAD

LabelImage Msc2D::ascending2Manifolds() {
    m_impl->ensureComputed();
    MyMscType* activeMsc = m_impl->activeMscOrThrow();
    const bool useImportedLineage =
        (m_impl->builderMode == BuilderMode::Partitioned && m_impl->partitionedMsc.get() != NULL);
    const bool emitDiagnostics = shouldEmitLabelDiagnostics();

    if (m_impl->baseLabelingAsc2.empty()) {
        m_impl->baseLabelingAsc2.assign(m_impl->grid->NumElements(), -1);
        m_impl->baseCompactAsc2.clear();
        activeMsc->SetSelectPersAbs(-1);

        const std::chrono::steady_clock::time_point t_base = std::chrono::steady_clock::now();
        bool gpuLabeled = false;
#ifdef MSC2D_HAS_GPU_DGRAD
        if (m_impl->useGpuGradient && !useImportedLineage) {
            gpuLabeled = gpuFillBaseLabeling(m_impl->mesh, m_impl->grad,
                                             m_impl->mX, m_impl->mY, activeMsc,
                                             true, m_impl->baseLabelingAsc2,
                                             m_impl->baseCompactAsc2);
            if (!gpuLabeled)
                printf("MSC2D: GPU base labeling (asc) unavailable; CPU paint\n");
        }
#endif
        if (!gpuLabeled) {
        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 0) continue;

            if (useImportedLineage) {
                const INDEX_TYPE rep_cell = activeMsc->getNode(nid).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!m_impl->partitionedMsc->GetAscendingLineageCells(rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    const int compact =
                        compactIdFor(m_impl->baseCompactAsc2, static_cast<int>(base_cells[ci]));
                    std::set<INDEX_TYPE> manifold;
                    activeMsc->rec_man_trace_up(base_cells[ci], manifold);
                    for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                        if (m_impl->mesh->dimension(*it) != 0) continue;
                        m_impl->baseLabelingAsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = compact;
                    }
                }
            } else {
                const int compact = compactIdFor(m_impl->baseCompactAsc2, static_cast<int>(nid));
                std::set<INDEX_TYPE> manifold;
                activeMsc->fillGeometry(nid, manifold, true);
                for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                    if (m_impl->mesh->dimension(*it) != 0) continue;
                    m_impl->baseLabelingAsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = compact;
                }
            }
        }
        }
        printf("TIMING: MSC2D base labeling asc (%s) ms=%lld\n",
               gpuLabeled ? "gpu" : "cpu", msSince(t_base));
    }

    activeMsc->SetSelectPersAbs(m_impl->selectedPersistence);

    std::unordered_map<int, int> remap;
    MyMscType::LivingNodesIterator nit(activeMsc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (activeMsc->getNode(nid).dim != 0) continue;

        if (useImportedLineage) {
            std::set<INT_TYPE> representative_constituents;
            activeMsc->GatherNodes(nid, representative_constituents, true);
            for (std::set<INT_TYPE>::const_iterator rcit = representative_constituents.begin();
                 rcit != representative_constituents.end(); ++rcit) {
                const INDEX_TYPE rep_cell = activeMsc->getNode(*rcit).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!m_impl->partitionedMsc->GetAscendingLineageCells(rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    remap[static_cast<int>(base_cells[ci])] = static_cast<int>(nid);
                }
            }
        } else {
            std::set<INT_TYPE> constituents;
            activeMsc->GatherNodes(nid, constituents, true);
            for (std::set<INT_TYPE>::const_iterator it = constituents.begin(); it != constituents.end(); ++it) {
                remap[static_cast<int>(*it)] = static_cast<int>(nid);
            }
        }
    }

    // Project the (raw base identity -> living node) map onto the compact id
    // space once, so the per-pixel loop below is an array read rather than a
    // hash probe. raw -> compact is injective, so no two remap keys land on the
    // same slot; a remap key that never labeled a pixel has no compact id and is
    // dropped, which cannot change the output.
    std::vector<int>& remapDense = m_impl->remapDenseAsc2;
    remapDense.assign(m_impl->baseCompactAsc2.size(), -1);
    for (std::unordered_map<int, int>::const_iterator it = remap.begin(); it != remap.end(); ++it) {
        const std::unordered_map<int, int>::const_iterator cit = m_impl->baseCompactAsc2.find(it->first);
        if (cit != m_impl->baseCompactAsc2.end()) {
            remapDense[static_cast<size_t>(cit->second)] = it->second;
        }
    }

    // Diagnostics-only: an unordered_set insert per pixel. Never compute this
    // unless the printf below is actually going to run.
    size_t base_unlabeled = 0;
    std::unordered_set<int> base_unique_ids;
    if (emitDiagnostics) {
        for (size_t i = 0; i < m_impl->baseLabelingAsc2.size(); ++i) {
            const int base = m_impl->baseLabelingAsc2[i];
            if (base < 0) {
                base_unlabeled++;
            } else {
                base_unique_ids.insert(base);
            }
        }
    }

    LabelImage out;
    out.width = m_impl->mX;
    out.height = m_impl->mY;
    out.labels.assign(static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY), -1);

    size_t remap_miss = 0;
    for (int i = 0; i < m_impl->mX * m_impl->mY; ++i) {
        const int base = m_impl->baseLabelingAsc2[i];
        if (base < 0) continue;
        const int living = remapDense[static_cast<size_t>(base)];
        if (living >= 0) {
            out.labels[static_cast<size_t>(i)] = living;
        } else {
            remap_miss++;
        }
    }
    if (emitDiagnostics) {
        const size_t total = static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY);
        const size_t final_unlabeled = base_unlabeled + remap_miss;
        printf("MSC2D asc_diag mode=%s total=%llu base_unlabeled=%llu remap_miss=%llu final_unlabeled=%llu base_unique_ids=%llu remap_keys=%llu\n",
            useImportedLineage ? "partitioned" : "serial",
            (unsigned long long)total,
            (unsigned long long)base_unlabeled,
            (unsigned long long)remap_miss,
            (unsigned long long)final_unlabeled,
            (unsigned long long)base_unique_ids.size(),
            (unsigned long long)remap.size());
    }

    return out;
}

LabelImage Msc2D::descending2Manifolds() {
    m_impl->ensureComputed();
    MyMscType* activeMsc = m_impl->activeMscOrThrow();
    const bool useImportedLineage =
        (m_impl->builderMode == BuilderMode::Partitioned && m_impl->partitionedMsc.get() != NULL);
    const bool emitDiagnostics = shouldEmitLabelDiagnostics();

    if (m_impl->baseLabelingDsc2.empty()) {
        m_impl->baseLabelingDsc2.assign(m_impl->grid->NumElements(), -1);
        m_impl->baseCompactDsc2.clear();
        activeMsc->SetSelectPersAbs(-1);

        const std::chrono::steady_clock::time_point t_base = std::chrono::steady_clock::now();
        bool gpuLabeled = false;
#ifdef MSC2D_HAS_GPU_DGRAD
        if (m_impl->useGpuGradient && !useImportedLineage) {
            gpuLabeled = gpuFillBaseLabeling(m_impl->mesh, m_impl->grad,
                                             m_impl->mX, m_impl->mY, activeMsc,
                                             false, m_impl->baseLabelingDsc2,
                                             m_impl->baseCompactDsc2);
            if (!gpuLabeled)
                printf("MSC2D: GPU base labeling (dsc) unavailable; CPU paint\n");
        }
#endif
        if (!gpuLabeled) {
        MyMscType::LivingNodesIterator nit(activeMsc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (activeMsc->getNode(nid).dim != 2) continue;

            if (useImportedLineage) {
                const INDEX_TYPE rep_cell = activeMsc->getNode(nid).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!m_impl->partitionedMsc->GetDescendingLineageCells(rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    const int compact =
                        compactIdFor(m_impl->baseCompactDsc2, static_cast<int>(base_cells[ci]));
                    std::set<INDEX_TYPE> manifold;
                    activeMsc->rec_man_trace_down(base_cells[ci], manifold);
                    for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                        if (m_impl->mesh->dimension(*it) != 2) continue;
                        m_impl->baseLabelingDsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = compact;
                    }
                }
            } else {
                const int compact = compactIdFor(m_impl->baseCompactDsc2, static_cast<int>(nid));
                std::set<INDEX_TYPE> manifold;
                activeMsc->fillGeometry(nid, manifold, false);
                for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin(); it != manifold.end(); ++it) {
                    if (m_impl->mesh->dimension(*it) != 2) continue;
                    m_impl->baseLabelingDsc2[m_impl->mesh->VertexNumberFromCellID(*it)] = compact;
                }
            }
        }
        }
        printf("TIMING: MSC2D base labeling dsc (%s) ms=%lld\n",
               gpuLabeled ? "gpu" : "cpu", msSince(t_base));
    }

    activeMsc->SetSelectPersAbs(m_impl->selectedPersistence);

    std::unordered_map<int, int> remap;
    MyMscType::LivingNodesIterator nit(activeMsc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (activeMsc->getNode(nid).dim != 2) continue;

        if (useImportedLineage) {
            std::set<INT_TYPE> representative_constituents;
            activeMsc->GatherNodes(nid, representative_constituents, false);
            for (std::set<INT_TYPE>::const_iterator rcit = representative_constituents.begin();
                 rcit != representative_constituents.end(); ++rcit) {
                const INDEX_TYPE rep_cell = activeMsc->getNode(*rcit).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!m_impl->partitionedMsc->GetDescendingLineageCells(rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    remap[static_cast<int>(base_cells[ci])] = static_cast<int>(nid);
                }
            }
        } else {
            std::set<INT_TYPE> constituents;
            activeMsc->GatherNodes(nid, constituents, false);
            for (std::set<INT_TYPE>::const_iterator it = constituents.begin(); it != constituents.end(); ++it) {
                remap[static_cast<int>(*it)] = static_cast<int>(nid);
            }
        }
    }

    // See ascending2Manifolds() for why this projection is output-preserving.
    std::vector<int>& remapDense = m_impl->remapDenseDsc2;
    remapDense.assign(m_impl->baseCompactDsc2.size(), -1);
    for (std::unordered_map<int, int>::const_iterator it = remap.begin(); it != remap.end(); ++it) {
        const std::unordered_map<int, int>::const_iterator cit = m_impl->baseCompactDsc2.find(it->first);
        if (cit != m_impl->baseCompactDsc2.end()) {
            remapDense[static_cast<size_t>(cit->second)] = it->second;
        }
    }

    // Diagnostics-only: an unordered_set insert per pixel. Never compute this
    // unless the printf below is actually going to run.
    size_t base_unlabeled = 0;
    std::unordered_set<int> base_unique_ids;
    if (emitDiagnostics) {
        for (size_t i = 0; i < m_impl->baseLabelingDsc2.size(); ++i) {
            const int base = m_impl->baseLabelingDsc2[i];
            if (base < 0) {
                base_unlabeled++;
            } else {
                base_unique_ids.insert(base);
            }
        }
    }

    LabelImage out;
    out.width = m_impl->mX;
    out.height = m_impl->mY;
    out.labels.assign(static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY), -1);

    size_t remap_miss = 0;
    for (int i = 0; i < m_impl->mX * m_impl->mY; ++i) {
        const int base = m_impl->baseLabelingDsc2[i];
        if (base < 0) continue;
        const int living = remapDense[static_cast<size_t>(base)];
        if (living >= 0) {
            out.labels[static_cast<size_t>(i)] = living;
        } else {
            remap_miss++;
        }
    }
    if (emitDiagnostics) {
        const size_t total = static_cast<size_t>(m_impl->mX) * static_cast<size_t>(m_impl->mY);
        const size_t final_unlabeled = base_unlabeled + remap_miss;
        printf("MSC2D dsc_diag mode=%s total=%llu base_unlabeled=%llu remap_miss=%llu final_unlabeled=%llu base_unique_ids=%llu remap_keys=%llu\n",
            useImportedLineage ? "partitioned" : "serial",
            (unsigned long long)total,
            (unsigned long long)base_unlabeled,
            (unsigned long long)remap_miss,
            (unsigned long long)final_unlabeled,
            (unsigned long long)base_unique_ids.size(),
            (unsigned long long)remap.size());
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

std::vector<ArcGeometry> Msc2D::arcGeometry() const {
    m_impl->ensureComputed();
    MyMscType* activeMsc = m_impl->activeMscOrThrow();

    std::vector<ArcGeometry> output;
    std::vector<INDEX_TYPE> cells;
    MyMscType::LivingArcsIterator ait(activeMsc);
    for (ait.begin(); ait.valid(); ait.advance()) {
        const INT_TYPE aid = ait.value();

        ArcGeometry arcEntry;
        arcEntry.id = static_cast<int>(aid);
        arcEntry.lowerNodeId = static_cast<int>(activeMsc->arcLowerNode(aid));
        arcEntry.upperNodeId = static_cast<int>(activeMsc->arcUpperNode(aid));
        arcEntry.dim = activeMsc->getArc(aid).dim;

        // Only realize the polyline geometry when it was actually constructed at
        // compute time (ComputeOptions.buildArcGeometry). Otherwise leave the
        // line empty -- tracing it here on nogeom arcs is both wrong and O(arc
        // length) per arc (the dominant cost when callers only need connectivity).
        cells.clear();
        if (m_impl->builtArcGeometry) {
            activeMsc->fillArcGeometry(aid, cells);
        } else {
                cells.push_back(activeMsc->getNode(activeMsc->arcLowerNode(aid)).cellindex);
                cells.push_back(activeMsc->getNode(activeMsc->arcUpperNode(aid)).cellindex);
        }
        arcEntry.line.reserve(cells.size());
        for (size_t i = 0; i < cells.size(); ++i) {
            GInt::Vec2l coords;
            m_impl->mesh->cellid2Coords(cells[i], coords);
            GInt::Vec2f fcoords = coords;
            fcoords *= 0.5f;
            arcEntry.line.push_back(Point{ fcoords[0], fcoords[1] });
        }
        output.push_back(arcEntry);
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
