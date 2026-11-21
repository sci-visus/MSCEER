#ifndef GPU_LINES2D_H
#define GPU_LINES2D_H

// 2D numeric separatrix lines on the GPU (labeling2d/lines2d.cu). Plain C++.

#include <cstdint>
#include "gpu_tables.h"
#include "tile2d.h"

namespace GInt {
namespace gpu {

// Stage A4: TILED numeric separatrix lines. For every critical edge the
// CPU's depth-10 BFS along its V-paths yields the seed vertices (descending
// lines, class 2) and seed quads (ascending lines, class 1); each seed is
// walked by numeric hops (lines2d.cuh) and digitized into per-cell class
// bits. The sink receives, per tile and per round, the owned cell box grown
// by one cell on every side; boxes OVERLAP and the sink MUST OR the bytes
// into its store (marks are idempotent bits). Walks leaving a tile are
// resumed in the neighbouring tile in later rounds until none is pending.
bool ComputeLineClasses2DTiled(const TileSource2D& values,
                               const ByteSource2D& gradSource,
                               const TileGrid2D& grid,
                               const Dgrad2DTables& tables,
                               const ByteSink2D& sinkOr,
                               const LineOptions2D& options,
                               LineStats2D* stats);

} // namespace gpu
} // namespace GInt

#endif
