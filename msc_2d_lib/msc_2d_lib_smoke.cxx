#include "msc_2d_lib.h"

#include <cstdint>
#include <iomanip>
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

// FNV-1a over the label image after renumbering labels by first occurrence in
// scan order (-1 stays -1).
//
// The raw node ids are NOT reproducible across processes: the discrete gradient
// and the partitioned build are threaded, so the same input yields the same
// regions with different ids on every run. Canonicalizing first makes the hash a
// stable fingerprint of the actual partition -- which pixels group together and
// where the background is -- so it can be diffed across builds and commits.
static uint64_t canonicalLabelHash(const std::vector<int>& labels) {
    uint64_t h = 14695981039346656037ULL;
    std::unordered_map<int, int> canonical;
    for (size_t i = 0; i < labels.size(); ++i) {
        int id = labels[i];
        if (id >= 0) {
            const std::unordered_map<int, int>::const_iterator it = canonical.find(id);
            if (it != canonical.end()) {
                id = it->second;
            } else {
                const int next = static_cast<int>(canonical.size());
                canonical.insert(std::make_pair(id, next));
                id = next;
            }
        }
        const uint32_t bits = static_cast<uint32_t>(id);
        for (int b = 0; b < 4; ++b) {
            h ^= static_cast<uint64_t>((bits >> (8 * b)) & 0xffu);
            h *= 1099511628211ULL;
        }
    }
    return h;
}

static void emitLabelHash(const char* mode, const char* dir, float pers,
                          const GInt::Msc2D::LabelImage& img) {
    const RegionStats stats = computeRegionStats(img);
    std::cout << "label_hash"
              << " mode=" << mode
              << " dir=" << dir
              << " pers=" << std::fixed << std::setprecision(6) << pers
              << std::defaultfloat
              << " width=" << img.width
              << " height=" << img.height
              << " size=" << img.labels.size()
              << " regions=" << stats.region_count
              << " unlabeled=" << stats.unlabeled_pixels
              << " canonical_hash=0x" << std::hex << std::setw(16) << std::setfill('0')
              << canonicalLabelHash(img.labels)
              << std::dec << std::setfill(' ')
              << std::endl;
}

// Re-threshold repeatedly against the cached base labeling, both directions. The
// non-monotonic tail re-exercises the remap after persistence drops back down.
static void sweepPersistences(GInt::Msc2D::Msc2D& msc, const char* mode,
                              const std::vector<float>& persistences) {
    for (size_t k = 0; k < persistences.size(); ++k) {
        msc.setPersistence(persistences[k]);
        emitLabelHash(mode, "asc", persistences[k], msc.ascending2Manifolds());
        emitLabelHash(mode, "dsc", persistences[k], msc.descending2Manifolds());
    }
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
    if (serialAscStats.unlabeled_pixels != 0) {
        std::cerr << "serial ascending manifolds contain unlabeled pixels: "
                  << serialAscStats.unlabeled_pixels << std::endl;
        return 1;
    }

    GInt::Msc2D::Msc2D::ComputeOptions partitionedOptions;
    partitionedOptions.builderMode = GInt::Msc2D::Msc2D::BuilderMode::Partitioned;
    partitionedOptions.requestedParallelism = 32; // Should clamp to 16.
    partitionedOptions.basePersistenceAbs = 0.0f;
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
    if (partitionedAscStats.unlabeled_pixels != 0) {
        std::cerr << "partitioned ascending manifolds contain unlabeled pixels: "
                  << partitionedAscStats.unlabeled_pixels << std::endl;
        return 6;
    }

    // Regression surface: both builder modes x both directions x several
    // persistences, hashed. The base labelings are cached by the first call in
    // each sweep, so this is exactly the interactive re-thresholding path.
    //
    // These runs pin parallelism to 1 on purpose. Threaded gradient computation
    // and threaded partition reconciliation are not reproducible run to run --
    // same input, same region counts, but a genuinely different partition each
    // time -- so a hash taken from the 8-thread/16-partition objects above would
    // be noise. At parallelism 1 the whole pipeline is deterministic and the
    // hashes below can be diffed across builds and commits.
    //
    // Partitioned mode still stores topological cell indices as base identities
    // at one partition (BuildPartitionLocalMSCs always runs, and the lineage
    // fallback uses node cellindex), so the two-id-space path stays covered.
    GInt::Msc2D::Msc2D::ComputeOptions serialDetOptions;
    serialDetOptions.requestedParallelism = 1;
    GInt::Msc2D::Msc2D serialDet;
    serialDet.compute(field.data(), rows, cols, serialDetOptions);

    std::vector<float> serialSweep;
    serialSweep.push_back(0.0f);
    serialSweep.push_back(0.05f);
    serialSweep.push_back(0.1f);
    serialSweep.push_back(0.2f);
    serialSweep.push_back(0.02f);
    sweepPersistences(serialDet, "serial", serialSweep);

    GInt::Msc2D::Msc2D::ComputeOptions partitionedDetOptions = partitionedOptions;
    partitionedDetOptions.requestedParallelism = 1;
    GInt::Msc2D::Msc2D partitionedDet;
    partitionedDet.compute(field.data(), rows, cols, partitionedDetOptions);

    std::vector<float> partitionedSweep;
    partitionedSweep.push_back(0.0f);
    partitionedSweep.push_back(0.01f);
    partitionedSweep.push_back(0.025f);
    partitionedSweep.push_back(0.05f);
    partitionedSweep.push_back(0.005f);
    sweepPersistences(partitionedDet, "partitioned", partitionedSweep);

    return 0;
}
