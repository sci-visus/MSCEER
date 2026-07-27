#ifndef MSC_3D_LIB_H
#define MSC_3D_LIB_H

#include <memory>
#include <vector>

namespace GInt {
namespace Msc3D {

struct Point {
    float x;
    float y;
    float z;
};

struct CriticalPoint {
    int id;
    float x;
    float y;
    float z;
    int dim;
    float value;
};

class Msc3D {
public:
    enum class BuilderMode {
        Serial,
        Partitioned
    };

    struct ComputeOptions {
        BuilderMode builderMode;
        int requestedPartitions;
        bool useExplicitSplits;
        int splitX;
        int splitY;
        int splitZ;
        int numThreads;
        float basePersistenceAbs;
        float cancelPersistenceAbs;

        ComputeOptions()
            : builderMode(BuilderMode::Serial),
              requestedPartitions(8),
              useExplicitSplits(false),
              splitX(1),
              splitY(1),
              splitZ(1),
              numThreads(8),
              basePersistenceAbs(-1.0f),
              cancelPersistenceAbs(-1.0f) {}
    };

    Msc3D();
    ~Msc3D();

    Msc3D(const Msc3D&) = delete;
    Msc3D& operator=(const Msc3D&) = delete;

    Msc3D(Msc3D&& other) noexcept;
    Msc3D& operator=(Msc3D&& other) noexcept;

    void compute(const float* rowMajorValues, int xdim, int ydim, int zdim, const ComputeOptions& options);
    void setPersistence(float value);
    std::vector<CriticalPoint> criticalPoints() const;

    bool hasResult() const;
    int xdim() const;
    int ydim() const;
    int zdim() const;
    int effectivePartitionCount() const;
    Point effectivePartitionSplits() const;
    int effectiveThreadCount() const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Msc3D
} // namespace GInt

#endif
