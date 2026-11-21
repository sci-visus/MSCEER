#include "msc_3d_lib.h"

#include <iostream>
#include <random>
#include <vector>

int main() {
    const int xdim = 40;
    const int ydim = 24;
    const int zdim = 20;
    std::vector<float> field((size_t)xdim * (size_t)ydim * (size_t)zdim, 0.0f);
    std::mt19937 rng(123456789u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (size_t i = 0; i < field.size(); ++i) {
        field[i] = dist(rng);
    }

    GInt::Msc3D::Msc3D serial;
    GInt::Msc3D::Msc3D::ComputeOptions serialOpts;
    serialOpts.builderMode = GInt::Msc3D::Msc3D::BuilderMode::Serial;
    serial.compute(field.data(), xdim, ydim, zdim, serialOpts);
    serial.setPersistence(0.1f);
    const std::vector<GInt::Msc3D::CriticalPoint> serialCritical = serial.criticalPoints();
    if (!serial.hasResult() || serialCritical.empty()) {
        std::cerr << "msc_3d_lib serial compute failed" << std::endl;
        return 1;
    }

    GInt::Msc3D::Msc3D partitioned;
    GInt::Msc3D::Msc3D::ComputeOptions partOpts;
    partOpts.builderMode = GInt::Msc3D::Msc3D::BuilderMode::Partitioned;
    partOpts.requestedPartitions = 12;
    partOpts.numThreads = 5;
    partOpts.basePersistenceAbs = 0.0f;
    partOpts.cancelPersistenceAbs = 0.15f;
    partitioned.compute(field.data(), xdim, ydim, zdim, partOpts);
    partitioned.setPersistence(partOpts.cancelPersistenceAbs);
    if (!partitioned.hasResult() || partitioned.criticalPoints().empty()) {
        std::cerr << "msc_3d_lib partitioned(number) compute failed" << std::endl;
        return 2;
    }
    const GInt::Msc3D::Point inferred = partitioned.effectivePartitionSplits();
    const int inferredProduct = static_cast<int>(inferred.x * inferred.y * inferred.z);
    if (partitioned.effectivePartitionCount() != 12 || inferredProduct != 12) {
        std::cerr << "msc_3d_lib inferred split/count mismatch" << std::endl;
        return 3;
    }

    GInt::Msc3D::Msc3D explicitPartitioned;
    GInt::Msc3D::Msc3D::ComputeOptions explicitOpts;
    explicitOpts.builderMode = GInt::Msc3D::Msc3D::BuilderMode::Partitioned;
    explicitOpts.useExplicitSplits = true;
    explicitOpts.splitX = 4;
    explicitOpts.splitY = 3;
    explicitOpts.splitZ = 1;
    explicitOpts.numThreads = 2;
    explicitOpts.basePersistenceAbs = 0.0f;
    explicitOpts.cancelPersistenceAbs = 0.12f;
    explicitPartitioned.compute(field.data(), xdim, ydim, zdim, explicitOpts);
    if (!explicitPartitioned.hasResult()) {
        std::cerr << "msc_3d_lib partitioned(tuple) compute failed" << std::endl;
        return 4;
    }
    const GInt::Msc3D::Point explicitSplits = explicitPartitioned.effectivePartitionSplits();
    if (explicitPartitioned.effectivePartitionCount() != 12 ||
        explicitSplits.x != 4.0f || explicitSplits.y != 3.0f || explicitSplits.z != 1.0f) {
        std::cerr << "msc_3d_lib explicit split mismatch" << std::endl;
        return 5;
    }

    std::cout << "msc_3d_lib_smoke"
              << " serial_critical_points=" << serialCritical.size()
              << " partitioned_requested_partitions=" << partOpts.requestedPartitions
              << " partitioned_effective_partitions=" << partitioned.effectivePartitionCount()
              << " inferred_splits=(" << inferred.x << "," << inferred.y << "," << inferred.z << ")"
              << " explicit_effective_partitions=" << explicitPartitioned.effectivePartitionCount()
              << " explicit_splits=(" << explicitSplits.x << "," << explicitSplits.y << "," << explicitSplits.z << ")"
              << std::endl;
    return 0;
}
