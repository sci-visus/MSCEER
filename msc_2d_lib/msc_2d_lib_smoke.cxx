#include "msc_2d_lib.h"

#include <iostream>
#include <random>
#include <vector>

int main() {
    const int rows = 512;
    const int cols = 512;

    std::vector<float> field(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0.0f);
    std::mt19937 rng(123456789u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    const float cx = 0.5f * static_cast<float>(cols);
    const float cy = 0.5f * static_cast<float>(rows);
    const float radius = 400.0f;
    const float r2 = radius * radius;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const float dx = static_cast<float>(c) - cx;
            const float dy = static_cast<float>(r) - cy;
            float v = dist(rng);
            if (dx * dx + dy * dy > r2) {
                v = 0.0f;
            }
            field[static_cast<size_t>(r) * static_cast<size_t>(cols) + static_cast<size_t>(c)] = v;
        }
    }

    GInt::Msc2D::Msc2D serial;
    serial.compute(field.data(), rows, cols, true, true);
    serial.setPersistence(0.1f);
    serial.computePolylineGraph(true);
    const std::vector<GInt::Msc2D::CriticalPoint> serialCritical = serial.criticalPoints();
    const GInt::Msc2D::Graph serialGraph = serial.graph();
    const GInt::Msc2D::LabelImage serialAsc = serial.ascending2Manifolds();
    const GInt::Msc2D::LabelImage serialDsc = serial.descending2Manifolds();

    std::cout << "serial_summary"
              << " critical_points=" << serialCritical.size()
              << " graph_nodes=" << serialGraph.nodes.size()
              << " graph_edges=" << serialGraph.edges.size()
              << " width=" << serial.width()
              << " height=" << serial.height()
              << std::endl;

    if (!serial.hasResult() || serial.width() != cols || serial.height() != rows) {
        std::cerr << "serial result metadata mismatch" << std::endl;
        return 1;
    }
    if (serialAsc.labels.size() != field.size() || serialDsc.labels.size() != field.size()) {
        std::cerr << "serial manifold label output size mismatch" << std::endl;
        return 1;
    }

    GInt::Msc2D::Msc2D::ComputeOptions partitionedOptions;
    partitionedOptions.builderMode = GInt::Msc2D::Msc2D::BuilderMode::Partitioned;
    partitionedOptions.requestedParallelism = 32; // Should clamp to 16.
    partitionedOptions.basePersistenceAbs = 0.02f;
    partitionedOptions.cancelPersistenceAbs = 0.05f;
    partitionedOptions.accurateAsc = true;
    partitionedOptions.accurateDsc = true;

    GInt::Msc2D::Msc2D partitioned;
    partitioned.compute(field.data(), rows, cols, partitionedOptions);
    partitioned.setPersistence(partitionedOptions.cancelPersistenceAbs);
    partitioned.computePolylineGraph(false);
    const std::vector<GInt::Msc2D::CriticalPoint> partitionedCritical = partitioned.criticalPoints();
    const GInt::Msc2D::Graph partitionedGraph = partitioned.graph();
    const GInt::Msc2D::LabelImage partitionedAsc = partitioned.ascending2Manifolds();
    const GInt::Msc2D::LabelImage partitionedDsc = partitioned.descending2Manifolds();
    const int clamped = partitioned.effectiveParallelism();

    std::cout << "partitioned_summary"
              << " requested_parallelism=" << partitionedOptions.requestedParallelism
              << " effective_parallelism=" << clamped
              << " base_persistence=" << partitionedOptions.basePersistenceAbs
              << " cancel_persistence=" << partitionedOptions.cancelPersistenceAbs
              << " critical_points=" << partitionedCritical.size()
              << " graph_nodes=" << partitionedGraph.nodes.size()
              << " graph_edges=" << partitionedGraph.edges.size()
              << " width=" << partitioned.width()
              << " height=" << partitioned.height()
              << std::endl;

    if (clamped != 16) {
        std::cerr << "partitioned clamping mismatch: expected 16 got " << clamped << std::endl;
        return 2;
    }
    if (!partitioned.hasResult() || partitioned.width() != cols || partitioned.height() != rows) {
        std::cerr << "partitioned result metadata mismatch" << std::endl;
        return 3;
    }
    if (partitionedAsc.labels.size() != field.size() || partitionedDsc.labels.size() != field.size()) {
        std::cerr << "partitioned manifold label output size mismatch" << std::endl;
        return 4;
    }

    return 0;
}
