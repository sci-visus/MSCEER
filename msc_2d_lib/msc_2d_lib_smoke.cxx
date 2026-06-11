#include "msc_2d_lib.h"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    const int rows = 32;
    const int cols = 32;

    std::vector<float> field(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0.0f);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const float x = (static_cast<float>(c) - 16.0f) / 8.0f;
            const float y = (static_cast<float>(r) - 16.0f) / 8.0f;
            field[static_cast<size_t>(r) * static_cast<size_t>(cols) + static_cast<size_t>(c)] =
                std::exp(-(x * x + y * y)) + 0.15f * std::sin(3.0f * x) - 0.1f * std::cos(2.0f * y);
        }
    }

    GInt::Msc2D::Msc2D msc;
    msc.compute(field.data(), rows, cols, true, true);
    msc.setPersistence(0.0f);

    const auto crit = msc.criticalPoints();
    const auto asc = msc.ascending2Manifolds();
    const auto dsc = msc.descending2Manifolds();
    msc.computePolylineGraph(false);
    const auto graph = msc.graph();

    std::cout << "critical_points=" << crit.size()
              << " asc_labels=" << asc.labels.size()
              << " dsc_labels=" << dsc.labels.size()
              << " graph_nodes=" << graph.nodes.size()
              << " graph_edges=" << graph.edges.size()
              << std::endl;

    return 0;
}
