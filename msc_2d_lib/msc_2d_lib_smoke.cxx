#include "msc_2d_lib.h"

#include <limits>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>

struct RegionStats {
    size_t region_count;
    size_t labeled_pixels;
    size_t unlabeled_pixels;
    size_t min_area;
    size_t max_area;
    double mean_area;
};

static RegionStats computeRegionStats(const GInt::Msc2D::LabelImage& labels) {
    std::unordered_map<int, size_t> areaByLabel;
    size_t unlabeled = 0;
    for (size_t i = 0; i < labels.labels.size(); ++i) {
        const int id = labels.labels[i];
        if (id < 0) {
            ++unlabeled;
            continue;
        }
        ++areaByLabel[id];
    }

    RegionStats stats;
    stats.region_count = areaByLabel.size();
    stats.unlabeled_pixels = unlabeled;
    stats.labeled_pixels = labels.labels.size() - unlabeled;
    stats.min_area = 0;
    stats.max_area = 0;
    stats.mean_area = 0.0;

    if (!areaByLabel.empty()) {
        size_t totalArea = 0;
        size_t minArea = std::numeric_limits<size_t>::max();
        size_t maxArea = 0;
        for (std::unordered_map<int, size_t>::const_iterator it = areaByLabel.begin();
             it != areaByLabel.end(); ++it) {
            const size_t area = it->second;
            totalArea += area;
            if (area < minArea) minArea = area;
            if (area > maxArea) maxArea = area;
        }
        stats.min_area = minArea;
        stats.max_area = maxArea;
        stats.mean_area = static_cast<double>(totalArea) / static_cast<double>(areaByLabel.size());
    }
    return stats;
}

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
    const RegionStats serialAscStats = computeRegionStats(serialAsc);
    std::cout << "serial_asc_regions"
              << " count=" << serialAscStats.region_count
              << " labeled_pixels=" << serialAscStats.labeled_pixels
              << " unlabeled_pixels=" << serialAscStats.unlabeled_pixels
              << " min_area=" << serialAscStats.min_area
              << " max_area=" << serialAscStats.max_area
              << " mean_area=" << serialAscStats.mean_area
              << std::endl;
    if (serialAscStats.region_count == 0 || serialAscStats.labeled_pixels == 0) {
        std::cerr << "serial ascending manifold region stats invalid" << std::endl;
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
    const RegionStats partitionedAscStats = computeRegionStats(partitionedAsc);
    std::cout << "partitioned_asc_regions"
              << " count=" << partitionedAscStats.region_count
              << " labeled_pixels=" << partitionedAscStats.labeled_pixels
              << " unlabeled_pixels=" << partitionedAscStats.unlabeled_pixels
              << " min_area=" << partitionedAscStats.min_area
              << " max_area=" << partitionedAscStats.max_area
              << " mean_area=" << partitionedAscStats.mean_area
              << std::endl;
    if (partitionedAscStats.region_count == 0 || partitionedAscStats.labeled_pixels == 0) {
        std::cerr << "partitioned ascending manifold region stats invalid" << std::endl;
        return 5;
    }

    return 0;
}
