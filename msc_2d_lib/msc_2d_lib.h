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

class Msc2D {
public:
    enum class BuilderMode {
        Serial,
        Partitioned
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
