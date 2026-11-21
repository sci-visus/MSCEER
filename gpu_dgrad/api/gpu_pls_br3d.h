#ifndef GPU_PLS_BR3D_H
#define GPU_PLS_BR3D_H

// PLS-BR (boundary-restricted PLS) on the GPU (dmt3d/). Plain C++. Output is memcmp-identical to the CPU driver
// ComputePLSBRGradient3D.

#include <cstdint>
#include "gpu_dmt3d.h"

namespace GInt {
namespace gpu {

// PLS-BR (boundary-restricted PLS), one thread per vertex. out_grad is required; out_step and out_part are the
// optional per-cell journal arrays. Returns false on a CUDA failure; unlike the CPU
// driver it returns true when cells were left unassigned (stats->unassigned reports
// them, 0 on every valid input).
inline bool ComputePLSBRGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
        uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
        DmtLowerStarStats* stats, Dmt3DTimings* timings) {
    return LowerStarGradient3DGPU(ranks, X, Y, Z, true, out_grad, out_step, out_part, stats, timings);
}

// Tiled version with 32-bit tile-local indexing (see gpu_dmt3d.h).
inline bool ComputePLSBRGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
        const ByteSink3D& gradSink, const DmtTiledOptions3D& options,
        DmtLowerStarStats* stats, TileStats3D* tstats) {
    return LowerStarGradient3DTiled(source, grid, true, gradSink, options, stats, tstats);
}

// Streamed version for a volume in host memory (see DmtStreamOptions3D).
inline bool ComputePLSBRGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
        uint8_t* out_grad, const DmtStreamOptions3D& options,
        DmtLowerStarStats* stats, TileStats3D* tstats) {
    return LowerStarGradient3DStreamed(ranks, X, Y, Z, true, out_grad, options, stats, tstats);
}

} // namespace gpu
} // namespace GInt

#endif
