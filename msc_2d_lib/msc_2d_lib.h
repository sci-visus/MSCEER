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

struct LabelImage {
    int width;
    int height;
    std::vector<int> labels;
};

class Msc2D {
public:
    Msc2D();
    ~Msc2D();

    Msc2D(const Msc2D&) = delete;
    Msc2D& operator=(const Msc2D&) = delete;

    Msc2D(Msc2D&& other) noexcept;
    Msc2D& operator=(Msc2D&& other) noexcept;

    void compute(const float* rowMajorValues, int rows, int cols, bool accurateAsc = true, bool accurateDsc = true);
    void setPersistence(float value);
    LabelImage ascending2Manifolds();
    LabelImage descending2Manifolds();
    std::vector<CriticalPoint> criticalPoints() const;
    void computePolylineGraph(bool useValleys);
    Graph graph() const;

    bool hasResult() const;
    int width() const;
    int height() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Msc2D
} // namespace GInt

#endif
