#ifndef MSC_2D_LIB_H
#define MSC_2D_LIB_H

#include <memory>
#include <vector>

namespace GInt {
namespace Msc2D {

struct Point {
    float x;
    float y;
};

struct Node {
    int id;
    std::vector<Point> geometry;
    std::vector<int> edges;
};

struct Edge {
    int id;
    int from;
    int to;
    std::vector<Point> geometry;
};

struct Graph {
    std::vector<Node> nodes;
    std::vector<Edge> edges;
};

struct CriticalPoint {
    int id;
    float x;
    float y;
    int dim;
    float value;
};

struct ArcGeometry {
    int id;
    int lowerNodeId;
    int upperNodeId;
    int dim;
    std::vector<Point> line;
};

struct LabelImage {
    int width;
    int height;
    std::vector<int> labels;
};

// One adjacency between two LIVING regions at the current persistence: the two
// regions are joined by at least one saddle. a < b, never a == b.
struct RegionArc {
    int a;
    int b;
    float saddleValue;  // most extreme saddle joining the pair
    int count;          // number of base arcs collapsed onto this pair
};

class Msc2D {
public:
    enum class BuilderMode {
        Serial,
        Partitioned
    };

    // How the persistence simplification is represented.
    //
    // MscHierarchy: build the Morse-Smale complex and its cancellation
    //   hierarchy; thresholding walks merged_manifold trees. Exact, and the
    //   only mode that can also answer criticalPoints()/arcGeometry()/graph()
    //   without extra work.
    //
    // MergeForest: skip the MSC entirely. Extrema come from the flow-terminal
    //   base labeling, saddles connecting them are read off the same terminal
    //   maps, and one persistence-ordered union-find pass records a merge
    //   forest that answers EVERY threshold by a parent-chain climb. Much
    //   cheaper to build (no MSC nodes/arcs, no cancellation loop, no
    //   partition reconcile), but the segmentation is not bit-identical to
    //   MscHierarchy: the MSC rewires arcs as it cancels while the forest's
    //   arc set is fixed, which measures ~99.7-99.9% pixel agreement.
    //   criticalPoints()/arcGeometry()/computePolylineGraph()/graph() still
    //   work in this mode -- they build the MSC lazily on first call, so a
    //   caller that only segments never pays for it.
    enum class Simplification {
        MscHierarchy,
        MergeForest
    };

    struct ComputeOptions {
        BuilderMode builderMode;
        int requestedParallelism;
        float basePersistenceAbs;
        float cancelPersistenceAbs;
        bool accurateAsc;
        bool accurateDsc;
        // Per-arc-dimension geometry construction (maps to Vec3b). Defaults off.
        // In partitioned mode this yields base-arc geometry (see msc_2d_lib.cxx).
        bool buildArcGeometry[3] = { false, false, false };
        // Produce the base discrete gradient on the GPU (requires the library
        // to be built with GPU_DGRAD_ENABLED and a usable CUDA device; falls
        // back to the CPU path otherwise). The GPU gradient is bit-identical
        // to the CPU one, so all downstream results are unchanged.
        bool useGpuGradient = false;
        // See the Simplification enum. Defaults to the MSC hierarchy, so an
        // existing caller's results are unchanged.
        Simplification simplification = Simplification::MscHierarchy;

        ComputeOptions()
            : builderMode(BuilderMode::Serial),
              requestedParallelism(8),
              basePersistenceAbs(-1.0f),
              cancelPersistenceAbs(-1.0f),
              accurateAsc(true),
              accurateDsc(true) {}
    };

    Msc2D();
    ~Msc2D();

    Msc2D(const Msc2D&) = delete;
    Msc2D& operator=(const Msc2D&) = delete;

    Msc2D(Msc2D&& other) noexcept;
    Msc2D& operator=(Msc2D&& other) noexcept;

    void compute(const float* rowMajorValues, int rows, int cols, bool accurateAsc = true, bool accurateDsc = true);
    void compute(const float* rowMajorValues, int rows, int cols, const ComputeOptions& options);
    void setPersistence(float value);
    LabelImage ascending2Manifolds();
    LabelImage descending2Manifolds();
    std::vector<CriticalPoint> criticalPoints() const;
    // Living MSC arcs with their geometric realization (polyline of world coords).
    // Lines are empty unless arc geometry was requested via ComputeOptions.buildArcGeometry.
    std::vector<ArcGeometry> arcGeometry() const;
    void computePolylineGraph(bool useValleys);
    Graph graph() const;

    // --- Region-scale access -------------------------------------------
    // The manifold calls above hand back a fresh LabelImage every time. These
    // expose the same data at region granularity, so a caller that repaints
    // often (an interactive persistence slider, a segmentation editor) can
    // work with a base-region-indexed array of its own and paint once.

    // Number of distinct base regions in a direction; 0 until the base
    // labeling has been built (baseLabeling() or a manifold call builds it).
    int baseRegionCount(bool ascending) const;

    // Per-pixel compact base region id (-1 = unlabeled), width*height in the
    // same order as LabelImage::labels. Built on first use; it depends only on
    // the discrete gradient, so it survives setPersistence().
    const std::vector<int>& baseLabeling(bool ascending);

    // Compact base region id -> living MSC node id at the CURRENT persistence
    // (-1 = no living owner), sized baseRegionCount(). This is the projection
    // the manifold calls apply, exposed without painting an image.
    const std::vector<int>& baseToLiving(bool ascending);

    // Living-region adjacency at the CURRENT persistence. Ids a/b are in the
    // SAME id space baseToLiving(ascending) returns in this mode -- base region
    // ids in MergeForest mode, living MSC node ids in MscHierarchy mode -- and
    // every id returned is guaranteed to appear there. a < b, no self-loops,
    // one entry per unordered pair, sorted by (a, b). saddleValue is the most
    // extreme saddle joining the pair (LOWEST for ascending/minima, HIGHEST for
    // descending/maxima); count is how many base arcs collapsed onto the pair.
    //
    // The two modes need not agree entry-for-entry: the forest's arc set is
    // fixed while the MSC rewires arcs as it cancels -- the same source as the
    // documented ~99.7-99.9% pixel agreement. Cached per (direction,
    // persistence) like the baseToLiving remap, so an interactive persistence
    // slider pays one build per threshold, not one per query.
    const std::vector<RegionArc>& livingRegionArcs(bool ascending);

    // out_labels[i] = remap[baseLabeling[i]], with -1 propagating and a base
    // id outside [0, m) also yielding -1. remap values are CALLER-DEFINED ints
    // - node ids, class ids, colors - not necessarily node ids. out_labels
    // must hold width()*height() ints. Uses the GPU label context when one is
    // available, otherwise a parallel CPU gather.
    void paintLabels(bool ascending, const int* remap, int m, int* out_labels);

    // Device-resident variants, for callers that consume the labels on the GPU
    // and want to skip the round trip. Both return an opaque device pointer
    // (int32 per pixel) or nullptr when there is no GPU label context. The
    // pointer from paintLabelsDevice() is invalidated by the next paint or
    // base-labeling change in the same direction.
    const void* deviceBaseLabels(bool ascending) const;
    const void* paintLabelsDevice(bool ascending, const int* remap, int m);

    // Frees the GPU label context (device base labels + paint buffers, roughly
    // 2-3 label images of VRAM) while keeping every host-side result. For
    // callers juggling many live Msc2D objects (one per slice of a stack) that
    // want only the active one resident. The next paint on this object lazily
    // re-uploads (~one label image H2D); device pointers previously returned
    // by deviceBaseLabels()/paintLabelsDevice() are invalidated. No-op without
    // GPU support or when no context exists.
    void releaseGpuResources();

    bool hasResult() const;
    int width() const;
    int height() const;
    int effectiveParallelism() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Msc2D
} // namespace GInt

#endif
