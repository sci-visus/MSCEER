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

// Kill switch for the parallel walk-down remap build (Stage 8). Default on;
// MSC2D_PARALLEL_REMAP=0 restores the original GatherNodes build, which is
// what the equivalence gate compares against. Read per call, like the
// diagnostics flag, so a harness can toggle it between objects.
bool shouldUseParallelRemap() {
    const char* env = std::getenv("MSC2D_PARALLEL_REMAP");
    return !(env != NULL && env[0] == '0' && env[1] == 0);
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
    // Persistent GPU label context (GInt::gpu::LabelCtx2D*), held as void* so
    // the public header never sees gpu_dgrad. Null unless the build has GPU
    // support, compute() ran with useGpuGradient, and a device is usable.
    // Population is deliberately NOT tied to whether the base labeling itself
    // came from the GPU: a CPU-painted labeling still gets GPU relabels.
    void* labelCtx;
    bool baseLabelsOnDevice[2];              // [0] ascending, [1] descending
    // baseToLiving cache: the remap is a function of the base labeling and the
    // selected persistence only, so one build serves every repaint at a
    // threshold instead of one per manifold call.
    bool remapValid[2];
    float remapPersistence[2];
    size_t lastRemapKeys[2];                 // diagnostics only
    // Node id -> compact base region ids, CSR, per direction. Built once with
    // the base labeling; lets the remap build write straight into remapDense.
    std::vector<int> nodeCompactOffsets[2];
    std::vector<int> nodeCompactValues[2];

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
          useGpuGradient(false),
          labelCtx(NULL) {
        for (int d = 0; d < 2; ++d) {
            baseLabelsOnDevice[d] = false;
            remapValid[d] = false;
            remapPersistence[d] = 0.0f;
            lastRemapKeys[d] = 0;
        }
    }

    ~Impl() { destroyLabelCtx(); }

    void resetComputedState() {
        destroyLabelCtx();
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
        for (int d = 0; d < 2; ++d) {
            remapValid[d] = false;
            remapPersistence[d] = 0.0f;
            lastRemapKeys[d] = 0;
            nodeCompactOffsets[d].clear();
            nodeCompactValues[d].clear();
        }
    }

    // Direction-keyed views of the per-direction caches. The ascending and
    // descending paths differ only in these and in a cell-dimension filter.
    std::vector<int>& baseLabelingFor(bool ascending) {
        return ascending ? baseLabelingAsc2 : baseLabelingDsc2;
    }
    const std::vector<int>& baseLabelingFor(bool ascending) const {
        return ascending ? baseLabelingAsc2 : baseLabelingDsc2;
    }
    std::unordered_map<int, int>& baseCompactFor(bool ascending) {
        return ascending ? baseCompactAsc2 : baseCompactDsc2;
    }
    std::vector<int>& remapDenseFor(bool ascending) {
        return ascending ? remapDenseAsc2 : remapDenseDsc2;
    }

    bool useImportedLineage() const {
        return builderMode == Msc2D::BuilderMode::Partitioned &&
               partitionedMsc.get() != NULL;
    }

    bool lineageCells(bool ascending, INDEX_TYPE rep_cell,
                      std::vector<INDEX_TYPE>& out) const {
        if (partitionedMsc.get() == NULL) return false;
        return ascending ? partitionedMsc->GetAscendingLineageCells(rep_cell, out)
                         : partitionedMsc->GetDescendingLineageCells(rep_cell, out);
    }

    // Defined below, after the GPU helpers they call.
    void destroyLabelCtx();
    void uploadBaseLabels(bool ascending);
    void ensureBaseLabeling(bool ascending);
    void buildNodeCompactIndex(bool ascending);
    size_t buildRemapGather(bool ascending, std::vector<int>& remapDense);
    size_t buildRemapWalkDown(bool ascending, std::vector<int>& remapDense,
                              long long& conflicts);
    const std::vector<int>& buildBaseToLiving(bool ascending);
    void paintLabelsInto(bool ascending, const int* remap, int m, int* out_labels);
    void emitLabelDiagnostics(bool ascending, const LabelImage& out) const;
    LabelImage manifolds2D(bool ascending);

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
// validated cell-for-cell against the serial fillGeometry paint.
//
// The raw base identity matches the CPU paint's key space in each builder
// mode: NODE ids in serial mode (`cellKeys == false`; terminals are resolved
// to nodes via the living base extrema), and base CELL indices in partitioned
// mode (`cellKeys == true`; a flow terminal IS the base extremum's cell, the
// same identity GetAscending/DescendingLineageCells hands the remap). Compact
// id numbering may differ from the CPU paint's iteration order, which is
// invisible downstream (outputs carry living node ids, not compacts).
// Returns false (nothing written) to fall back to the CPU paint.
bool gpuFillBaseLabeling(Accurate2D::MeshType* mesh, Accurate2D::GradType* grad,
                         int mX, int mY, MyMscType* activeMsc, bool ascending,
                         bool cellKeys,
                         std::vector<int>& baseLabeling,
                         std::unordered_map<int, int>& baseCompact) {
    gpu::Dgrad2DTables tables;
    if (!gpu::BuildDgrad2DTablesFromMesh(*mesh, tables)) return false;
    const uint8_t* gradBytes =
        reinterpret_cast<const uint8_t*>(grad->m_dgrad->LabelArray());
    const int64_t X = mX;
    const int64_t Y = mY;
    const int64_t cellRow = 2 * X - 1;         // cell-lattice row stride
    if (ascending) {
        std::vector<int32_t> term(static_cast<size_t>(X) * static_cast<size_t>(Y));
        if (!gpu::Label2DFlowTerminals(gradBytes, X, Y, tables, term.data(), NULL, NULL))
            return false;
        // Terminal vertex number -> compact id, resolved lazily per distinct
        // terminal (thousands, not per pixel).
        std::vector<int> termCompact(term.size(), -1);
        if (!cellKeys) {
            MyMscType::LivingNodesIterator nit(activeMsc);
            for (nit.begin(); nit.valid(); nit.advance()) {
                const INT_TYPE nid = nit.value();
                if (activeMsc->getNode(nid).dim != 0) continue;
                const INDEX_TYPE cell = activeMsc->getNode(nid).cellindex;
                const INDEX_TYPE vnum = mesh->VertexNumberFromCellID(cell);
                termCompact[static_cast<size_t>(vnum)] =
                    compactIdFor(baseCompact, static_cast<int>(nid));
            }
        }
        for (size_t i = 0; i < term.size(); ++i) {
            const int32_t t = term[i];
            if (t < 0) continue;
            int compact = termCompact[static_cast<size_t>(t)];
            if (compact < 0 && cellKeys) {
                // A vertex (vx, vy) sits at even cell coords (2vx, 2vy).
                const int64_t vx = t % X;
                const int64_t vy = t / X;
                compact = compactIdFor(baseCompact,
                                       static_cast<int>(2 * vx + 2 * vy * cellRow));
                termCompact[static_cast<size_t>(t)] = compact;
            }
            if (compact >= 0) baseLabeling[i] = compact;
        }
    } else {
        const int64_t QX = X - 1;
        const int64_t QY = Y - 1;
        if (QX <= 0 || QY <= 0) return false;
        std::vector<int32_t> qterm(static_cast<size_t>(QX) * static_cast<size_t>(QY));
        if (!gpu::Label2DFlowTerminals(gradBytes, X, Y, tables, NULL, qterm.data(), NULL))
            return false;
        // Terminal quad lattice index -> compact id. A quad cell has odd
        // coordinates (2qx+1, 2qy+1); its quad lattice index is qx + qy*QX.
        std::vector<int> termCompact(qterm.size(), -1);
        if (!cellKeys) {
            MyMscType::LivingNodesIterator nit(activeMsc);
            for (nit.begin(); nit.valid(); nit.advance()) {
                const INT_TYPE nid = nit.value();
                if (activeMsc->getNode(nid).dim != 2) continue;
                const INDEX_TYPE cell = activeMsc->getNode(nid).cellindex;
                const int64_t cx = static_cast<int64_t>(cell) % cellRow;
                const int64_t cy = static_cast<int64_t>(cell) / cellRow;
                termCompact[static_cast<size_t>((cx >> 1) + (cy >> 1) * QX)] =
                    compactIdFor(baseCompact, static_cast<int>(nid));
            }
        }
        // Project each quad's label onto the same vertex the CPU paint stamps
        // (VertexNumberFromCellID of the quad cell) -- the mapping is
        // injective, so stamping order cannot matter.
        for (int64_t qy = 0; qy < QY; ++qy) {
            for (int64_t qx = 0; qx < QX; ++qx) {
                const int32_t t = qterm[static_cast<size_t>(qx + qy * QX)];
                if (t < 0) continue;
                int compact = termCompact[static_cast<size_t>(t)];
                if (compact < 0 && cellKeys) {
                    const int64_t tqx = t % QX;
                    const int64_t tqy = t / QX;
                    compact = compactIdFor(
                        baseCompact,
                        static_cast<int>((2 * tqx + 1) + (2 * tqy + 1) * cellRow));
                    termCompact[static_cast<size_t>(t)] = compact;
                }
                if (compact < 0) continue;
                const INDEX_TYPE cellid = (2 * qx + 1) + (2 * qy + 1) * cellRow;
                baseLabeling[static_cast<size_t>(mesh->VertexNumberFromCellID(cellid))] =
                    compact;
            }
        }
    }
    return true;
}
}  // namespace
#endif  // MSC2D_HAS_GPU_DGRAD

// ---------------------------------------------------------------------------
// Base labeling, base->living remap, and the per-pixel paint.
//
// ascending2Manifolds() and descending2Manifolds() used to be ~150 lines of
// near-verbatim duplicate differing only in a cell-dimension filter (0 vs 2)
// and in which cached vectors they touched. They are now thin wrappers over
// the four helpers below, which are also what the public region-scale API
// (baseLabeling / baseToLiving / paintLabels) is built from.
//
// The three stages have different lifetimes, and keeping them apart is the
// whole point:
//   * the base labeling depends only on the discrete gradient, so it is built
//     once per compute() and survives every setPersistence();
//   * the base->living remap depends additionally on the selected persistence,
//     so it is cached per (direction, persistence);
//   * the paint depends on the caller's remap, and is the only stage that has
//     to touch every pixel.
// ---------------------------------------------------------------------------

void Msc2D::Impl::destroyLabelCtx() {
#ifdef MSC2D_HAS_GPU_DGRAD
    if (labelCtx) gpu::DestroyLabelCtx2D(static_cast<gpu::LabelCtx2D*>(labelCtx));
#endif
    labelCtx = NULL;
    baseLabelsOnDevice[0] = false;
    baseLabelsOnDevice[1] = false;
}

// Mirrors the freshly built host base labeling onto the device so later
// relabels upload only the remap. Any failure just leaves the context absent,
// and paintLabelsInto falls back to the CPU gather.
void Msc2D::Impl::uploadBaseLabels(bool ascending) {
#ifdef MSC2D_HAS_GPU_DGRAD
    if (!useGpuGradient) return;
    const std::vector<int>& baseLabeling = baseLabelingFor(ascending);
    const long long n = static_cast<long long>(mX) * static_cast<long long>(mY);
    if (n <= 0 || static_cast<long long>(baseLabeling.size()) < n) return;
    if (!labelCtx) {
        labelCtx = gpu::CreateLabelCtx2D(mX, mY);
        if (!labelCtx) {
            printf("MSC2D: GPU label context unavailable; CPU relabel\n");
            return;
        }
    }
    const int d = ascending ? 0 : 1;
    baseLabelsOnDevice[d] =
        gpu::CtxSetBaseLabels(static_cast<gpu::LabelCtx2D*>(labelCtx), ascending,
                              reinterpret_cast<const int32_t*>(baseLabeling.data()));
    if (!baseLabelsOnDevice[d])
        printf("MSC2D: GPU base label upload (%s) failed; CPU relabel\n",
               ascending ? "asc" : "dsc");
#else
    (void)ascending;
#endif
}

// Builds the per-pixel compact base region labeling for one direction, once.
// Walks the BASE (unsimplified) complex, so it brackets the work with
// SetSelectPersAbs(-1) and restores the user's threshold on the way out.
void Msc2D::Impl::ensureBaseLabeling(bool ascending) {
    ensureComputed();
    MyMscType* msc = activeMscOrThrow();
    std::vector<int>& baseLabeling = baseLabelingFor(ascending);
    if (!baseLabeling.empty()) return;

    const bool importedLineage = useImportedLineage();
    const int wantDim = ascending ? 0 : 2;
    std::unordered_map<int, int>& baseCompact = baseCompactFor(ascending);

    baseLabeling.assign(grid->NumElements(), -1);
    baseCompact.clear();
    msc->SetSelectPersAbs(-1);

    const std::chrono::steady_clock::time_point t_base = std::chrono::steady_clock::now();
    bool gpuLabeled = false;
#ifdef MSC2D_HAS_GPU_DGRAD
    if (useGpuGradient) {
        gpuLabeled = gpuFillBaseLabeling(mesh, grad, mX, mY, msc, ascending,
                                         importedLineage, baseLabeling, baseCompact);
        if (!gpuLabeled)
            printf("MSC2D: GPU base labeling (%s) unavailable; CPU paint\n",
                   ascending ? "asc" : "dsc");
    }
#endif
    if (!gpuLabeled) {
        MyMscType::LivingNodesIterator nit(msc);
        for (nit.begin(); nit.valid(); nit.advance()) {
            const INT_TYPE nid = nit.value();
            if (msc->getNode(nid).dim != wantDim) continue;

            if (importedLineage) {
                const INDEX_TYPE rep_cell = msc->getNode(nid).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!lineageCells(ascending, rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    const int compact =
                        compactIdFor(baseCompact, static_cast<int>(base_cells[ci]));
                    std::set<INDEX_TYPE> manifold;
                    if (ascending) {
                        msc->rec_man_trace_up(base_cells[ci], manifold);
                    } else {
                        msc->rec_man_trace_down(base_cells[ci], manifold);
                    }
                    for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin();
                         it != manifold.end(); ++it) {
                        if (mesh->dimension(*it) != wantDim) continue;
                        baseLabeling[mesh->VertexNumberFromCellID(*it)] = compact;
                    }
                }
            } else {
                const int compact = compactIdFor(baseCompact, static_cast<int>(nid));
                std::set<INDEX_TYPE> manifold;
                msc->fillGeometry(nid, manifold, ascending);
                for (std::set<INDEX_TYPE>::const_iterator it = manifold.begin();
                     it != manifold.end(); ++it) {
                    if (mesh->dimension(*it) != wantDim) continue;
                    baseLabeling[mesh->VertexNumberFromCellID(*it)] = compact;
                }
            }
        }
    }
    printf("TIMING: MSC2D base labeling %s (%s) ms=%lld\n", ascending ? "asc" : "dsc",
           gpuLabeled ? "gpu" : "cpu", msSince(t_base));

    uploadBaseLabels(ascending);
    buildNodeCompactIndex(ascending);
    msc->SetSelectPersAbs(selectedPersistence);
}

// Node id -> compact base region ids, as a CSR built once per direction over
// the BASE complex. This is the piece that lets the remap build write straight
// into remapDense instead of accumulating an unordered_map and projecting it
// afterwards, and in partitioned mode it also hoists the per-constituent
// lineage vector copy out of the per-persistence work.
//
// Rows exist only for nodes of the direction's extremal dimension; everything
// else is empty. A raw id with no compact id is dropped, exactly as the
// projection step it replaces does, so this is output-preserving.
void Msc2D::Impl::buildNodeCompactIndex(bool ascending) {
    MyMscType* msc = activeMscOrThrow();
    const int d = ascending ? 0 : 1;
    const int wantDim = ascending ? 0 : 2;
    const bool importedLineage = useImportedLineage();
    const std::unordered_map<int, int>& baseCompact = baseCompactFor(ascending);

    std::vector<int>& offsets = nodeCompactOffsets[d];
    std::vector<int>& values = nodeCompactValues[d];
    const INT_TYPE n = msc->numNodes();
    offsets.assign(static_cast<size_t>(n) + 1, 0);
    values.clear();

    std::vector<INDEX_TYPE> base_cells;
    for (INT_TYPE nid = 0; nid < n; ++nid) {
        offsets[static_cast<size_t>(nid)] = static_cast<int>(values.size());
        if (msc->getNode(nid).dim != wantDim) continue;

        if (importedLineage) {
            const INDEX_TYPE rep_cell = msc->getNode(nid).cellindex;
            base_cells.clear();
            if (!lineageCells(ascending, rep_cell, base_cells)) {
                base_cells.push_back(rep_cell);
            }
            for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                const std::unordered_map<int, int>::const_iterator it =
                    baseCompact.find(static_cast<int>(base_cells[ci]));
                if (it != baseCompact.end()) values.push_back(it->second);
            }
        } else {
            const std::unordered_map<int, int>::const_iterator it =
                baseCompact.find(static_cast<int>(nid));
            if (it != baseCompact.end()) values.push_back(it->second);
        }
    }
    offsets[static_cast<size_t>(n)] = static_cast<int>(values.size());
}

// The original remap build: one GatherNodes into a std::set per living
// extremum, an unordered_map keyed by raw base identity, then a projection
// pass into the compact id space. Kept as the reference semantics and as the
// fallback for MSC2D_PARALLEL_REMAP=0.
size_t Msc2D::Impl::buildRemapGather(bool ascending, std::vector<int>& remapDense) {
    MyMscType* msc = activeMscOrThrow();
    const bool importedLineage = useImportedLineage();
    const int wantDim = ascending ? 0 : 2;
    const std::unordered_map<int, int>& baseCompact = baseCompactFor(ascending);

    std::unordered_map<int, int> remap;
    MyMscType::LivingNodesIterator nit(msc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (msc->getNode(nid).dim != wantDim) continue;

        std::set<INT_TYPE> constituents;
        msc->GatherNodes(nid, constituents, ascending);
        if (importedLineage) {
            for (std::set<INT_TYPE>::const_iterator rcit = constituents.begin();
                 rcit != constituents.end(); ++rcit) {
                const INDEX_TYPE rep_cell = msc->getNode(*rcit).cellindex;
                std::vector<INDEX_TYPE> base_cells;
                if (!lineageCells(ascending, rep_cell, base_cells)) {
                    base_cells.push_back(rep_cell);
                }
                for (size_t ci = 0; ci < base_cells.size(); ++ci) {
                    remap[static_cast<int>(base_cells[ci])] = static_cast<int>(nid);
                }
            }
        } else {
            for (std::set<INT_TYPE>::const_iterator it = constituents.begin();
                 it != constituents.end(); ++it) {
                remap[static_cast<int>(*it)] = static_cast<int>(nid);
            }
        }
    }

    // Project the (raw base identity -> living node) map onto the compact id
    // space. raw -> compact is injective, so no two remap keys land on the same
    // slot; a remap key that never labeled a pixel has no compact id and is
    // dropped, which cannot change the output.
    remapDense.assign(baseCompact.size(), -1);
    for (std::unordered_map<int, int>::const_iterator it = remap.begin();
         it != remap.end(); ++it) {
        const std::unordered_map<int, int>::const_iterator cit =
            baseCompact.find(it->first);
        if (cit != baseCompact.end()) {
            remapDense[static_cast<size_t>(cit->second)] = it->second;
        }
    }
    return remap.size();
}

// Parallel walk-down: for each living extremum, descend the merged-manifold
// hierarchy from its active manifold and stamp the living node id onto every
// base region underneath it. This is the same traversal recGatherNodes
// performs, so it inherits its semantics exactly - what changes is that there
// is no std::set, no unordered_map and no projection pass, the roots run
// concurrently, and the inner loop allocates nothing.
//
// Correctness rests on each base region having exactly one living owner, which
// holds for the ascending manifolds of minima and the descending manifolds of
// maxima. It does NOT hold in general: a cancellation can hand the same
// manifold to several surviving neighbours as merged[1], which is what makes
// the 1-manifolds (2D and 3D) and the 2-manifolds in 3D multi-membership - and
// out of scope here. Rather than trust the argument, the walk counts the
// double-writes it performs; a non-zero count makes the caller redo the build
// with the reference path, so a violation costs time, never correctness.
size_t Msc2D::Impl::buildRemapWalkDown(bool ascending, std::vector<int>& remapDense,
                                       long long& conflicts) {
    MyMscType* msc = activeMscOrThrow();
    const int d = ascending ? 0 : 1;
    const int wantDim = ascending ? 0 : 2;
    const std::vector<int>& offsets = nodeCompactOffsets[d];
    const std::vector<int>& values = nodeCompactValues[d];

    remapDense.assign(baseCompactFor(ascending).size(), -1);
    conflicts = 0;

    std::vector<int> roots;
    MyMscType::LivingNodesIterator nit(msc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        if (msc->getNode(nid).dim == wantDim) roots.push_back(static_cast<int>(nid));
    }

    const long long nroots = static_cast<long long>(roots.size());
    int* dense = remapDense.empty() ? NULL : &remapDense[0];
    long long conf = 0;
#pragma omp parallel
    {
        std::vector<INT_TYPE> stack;
#pragma omp for schedule(dynamic, 64) reduction(+ : conf)
        for (long long r = 0; r < nroots; ++r) {
            const int nid = roots[static_cast<size_t>(r)];
            stack.clear();
            stack.push_back(msc->activeManifoldForNode(nid, ascending));
            while (!stack.empty()) {
                const INT_TYPE man = stack.back();
                stack.pop_back();
                const merged_manifold& mm = msc->manifoldAt(man);
                if (mm.merged[0] != -1) {
                    stack.push_back(mm.merged[0]);
                    stack.push_back(mm.merged[1]);
                }
                const size_t b = static_cast<size_t>(mm.basenode);
                for (int k = offsets[b]; k < offsets[b + 1]; ++k) {
                    const int c = values[static_cast<size_t>(k)];
                    if (dense[c] >= 0 && dense[c] != nid) conf++;
                    dense[c] = nid;
                }
            }
        }
    }
    conflicts = conf;
    return roots.size();
}

// compact base region id -> living node id at the selected persistence.
const std::vector<int>& Msc2D::Impl::buildBaseToLiving(bool ascending) {
    const int d = ascending ? 0 : 1;
    std::vector<int>& remapDense = remapDenseFor(ascending);
    if (remapValid[d] && remapPersistence[d] == selectedPersistence) return remapDense;

    activeMscOrThrow()->SetSelectPersAbs(selectedPersistence);

    const std::chrono::steady_clock::time_point t_remap =
        std::chrono::steady_clock::now();

    const bool wantParallel = shouldUseParallelRemap() && !nodeCompactOffsets[d].empty();
    bool usedParallel = false;
    size_t keys = 0;
    if (wantParallel) {
        long long conflicts = 0;
        keys = buildRemapWalkDown(ascending, remapDense, conflicts);
        usedParallel = (conflicts == 0);
        if (!usedParallel) {
            printf("MSC2D: %lld ambiguous base regions in the %s walk-down; "
                   "rebuilding with GatherNodes\n",
                   conflicts, ascending ? "asc" : "dsc");
        }
    }
    if (!usedParallel) keys = buildRemapGather(ascending, remapDense);

    printf("TIMING: MSC2D remap build %s (%s) ms=%lld\n", ascending ? "asc" : "dsc",
           usedParallel ? "parallel" : "gather", msSince(t_remap));

    lastRemapKeys[d] = keys;
    remapPersistence[d] = selectedPersistence;
    remapValid[d] = true;
    return remapDense;
}

// out_labels[i] = remap[baseLabeling[i]], -1 propagating. A base id outside
// [0, m) also yields -1, matching the GPU contract exactly, so the two paths
// are interchangeable. remap may be null when m == 0.
void Msc2D::Impl::paintLabelsInto(bool ascending, const int* remap, int m,
                                  int* out_labels) {
    const long long n = static_cast<long long>(mX) * static_cast<long long>(mY);
    if (n <= 0 || !out_labels) return;

    bool painted = false;
#ifdef MSC2D_HAS_GPU_DGRAD
    const int d = ascending ? 0 : 1;
    if (labelCtx && baseLabelsOnDevice[d]) {
        painted = gpu::CtxRelabel(static_cast<gpu::LabelCtx2D*>(labelCtx), ascending,
                                  reinterpret_cast<const int32_t*>(remap), m,
                                  reinterpret_cast<int32_t*>(out_labels), NULL);
        if (!painted)
            printf("MSC2D: GPU relabel (%s) failed; CPU gather\n",
                   ascending ? "asc" : "dsc");
    }
#endif
    if (painted) return;

    const std::vector<int>& baseLabeling = baseLabelingFor(ascending);
    const int* base = baseLabeling.data();
#pragma omp parallel for
    for (long long i = 0; i < n; ++i) {
        const int b = base[i];
        out_labels[i] = (b < 0 || b >= m) ? -1 : remap[b];
    }
}

void Msc2D::Impl::emitLabelDiagnostics(bool ascending, const LabelImage& out) const {
    const std::vector<int>& baseLabeling = baseLabelingFor(ascending);
    const size_t total = static_cast<size_t>(mX) * static_cast<size_t>(mY);
    size_t base_unlabeled = 0;
    size_t remap_miss = 0;
    std::unordered_set<int> base_unique_ids;
    for (size_t i = 0; i < total; ++i) {
        const int base = baseLabeling[i];
        if (base < 0) {
            base_unlabeled++;
            continue;
        }
        base_unique_ids.insert(base);
        if (out.labels[i] < 0) remap_miss++;
    }
    printf("MSC2D %s_diag mode=%s total=%llu base_unlabeled=%llu remap_miss=%llu final_unlabeled=%llu base_unique_ids=%llu remap_keys=%llu\n",
        ascending ? "asc" : "dsc",
        useImportedLineage() ? "partitioned" : "serial",
        (unsigned long long)total,
        (unsigned long long)base_unlabeled,
        (unsigned long long)remap_miss,
        (unsigned long long)(base_unlabeled + remap_miss),
        (unsigned long long)base_unique_ids.size(),
        (unsigned long long)lastRemapKeys[ascending ? 0 : 1]);
}

LabelImage Msc2D::Impl::manifolds2D(bool ascending) {
    ensureBaseLabeling(ascending);
    const std::vector<int>& remapDense = buildBaseToLiving(ascending);

    LabelImage out;
    out.width = mX;
    out.height = mY;
    out.labels.assign(static_cast<size_t>(mX) * static_cast<size_t>(mY), -1);
    paintLabelsInto(ascending, remapDense.empty() ? NULL : &remapDense[0],
                    static_cast<int>(remapDense.size()), out.labels.data());

    if (shouldEmitLabelDiagnostics()) emitLabelDiagnostics(ascending, out);
    return out;
}

LabelImage Msc2D::ascending2Manifolds() {
    return m_impl->manifolds2D(true);
}

LabelImage Msc2D::descending2Manifolds() {
    return m_impl->manifolds2D(false);
}

int Msc2D::baseRegionCount(bool ascending) const {
    return static_cast<int>(ascending ? m_impl->baseCompactAsc2.size()
                                      : m_impl->baseCompactDsc2.size());
}

const std::vector<int>& Msc2D::baseLabeling(bool ascending) {
    m_impl->ensureBaseLabeling(ascending);
    return m_impl->baseLabelingFor(ascending);
}

const std::vector<int>& Msc2D::baseToLiving(bool ascending) {
    m_impl->ensureBaseLabeling(ascending);
    return m_impl->buildBaseToLiving(ascending);
}

void Msc2D::paintLabels(bool ascending, const int* remap, int m, int* out_labels) {
    m_impl->ensureBaseLabeling(ascending);
    m_impl->paintLabelsInto(ascending, remap, m, out_labels);
}

const void* Msc2D::deviceBaseLabels(bool ascending) const {
#ifdef MSC2D_HAS_GPU_DGRAD
    if (!m_impl->labelCtx || !m_impl->baseLabelsOnDevice[ascending ? 0 : 1]) return NULL;
    return gpu::CtxBaseLabelsDevice(
        static_cast<const gpu::LabelCtx2D*>(m_impl->labelCtx), ascending);
#else
    (void)ascending;
    return NULL;
#endif
}

const void* Msc2D::paintLabelsDevice(bool ascending, const int* remap, int m) {
#ifdef MSC2D_HAS_GPU_DGRAD
    m_impl->ensureBaseLabeling(ascending);
    if (!m_impl->labelCtx || !m_impl->baseLabelsOnDevice[ascending ? 0 : 1]) return NULL;
    const void* dev = NULL;
    if (!gpu::CtxRelabelDevice(static_cast<gpu::LabelCtx2D*>(m_impl->labelCtx), ascending,
                               reinterpret_cast<const int32_t*>(remap), m, &dev, NULL))
        return NULL;
    return dev;
#else
    (void)ascending;
    (void)remap;
    (void)m;
    return NULL;
#endif
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
