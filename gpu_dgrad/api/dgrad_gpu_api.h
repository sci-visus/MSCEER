#ifndef DGRAD_GPU_API_H
#define DGRAD_GPU_API_H

// Host-facing interface to the GPU discrete gradient path.
//
// This header is plain C++ (no CUDA types) so any existing MSCEER tool can
// include it and link against the gpu_dgrad library without a CUDA compiler.
//
// Output contract: the gradient is one byte per topological cell of the
// doubled grid (2X-1) x (2Y-1) x (2Z-1), row-major x-fastest, with the exact
// bit layout of GradBitfield in include/gi_discrete_gradient_labeling.h
// (assigned:1, flag:1, pair:3 with 7 = critical, ldir:3). The GPU path must
// produce output that memcmp-matches the CPU path
// (CLSSlidingGradient::ComputePairing_sliding) - the algorithm is
// comparisons-only under the value-then-index total order, so any mismatch is
// a bug, not floating-point noise.

// The interface is split by method; this header includes every part, so code that
// includes dgrad_gpu_api.h keeps compiling. New code can include only what it uses.

#include <cstdint>

#include "tile2d.h"
#include "tile3d.h"
#include "gpu_device.h"
#include "gpu_tables.h"
#include "gpu_cls2d.h"
#include "gpu_labeling2d.h"
#include "gpu_lines2d.h"
#include "gpu_cls3d.h"
#include "gpu_dmt3d.h"
#include "gpu_pls3d.h"
#include "gpu_pls_br3d.h"
#include "gpu_dmf3d.h"

namespace GInt {
namespace gpu {

// Earlier names of the lower-star (CLS) entry points. Kept so that existing callers keep
// compiling; new code calls ComputeCLSGradient2D* and ComputeCLSGradient3D*.
inline bool ComputeDiscreteGradient2D(const float* values,
        int64_t X,
        int64_t Y,
        const uint8_t* max_labels,
        const uint8_t* min_labels,
        const Dgrad2DTables& tables,
        uint8_t* out_grad,
        bool interior_only,
        Dgrad2DTimings* timings) {
    return ComputeCLSGradient2D(values, X, Y, max_labels, min_labels, tables, out_grad, interior_only, timings);
}

inline bool ComputeDiscreteGradient2DBatched(const float* values,
        int64_t X,
        int64_t Y,
        int64_t nslices,
        const Dgrad2DTables& tables,
        uint8_t* out_grad,
        int block_shape,
        Dgrad2DTimings* timings) {
    return ComputeCLSGradient2DBatched(values, X, Y, nslices, tables, out_grad, block_shape, timings);
}

inline bool ComputeDiscreteGradient2DFused(const float* values,
        int64_t X,
        int64_t Y,
        const Dgrad2DTables& tables,
        uint8_t* out_grad,
        Dgrad2DTimings* timings) {
    return ComputeCLSGradient2DFused(values, X, Y, tables, out_grad, timings);
}

inline bool ComputeDiscreteGradient2DConstrained(const float* values,
        const uint8_t* classes,
        int64_t X,
        int64_t Y,
        const Dgrad2DTables& tables,
        uint8_t* out_grad,
        Dgrad2DTimings* timings) {
    return ComputeCLSGradient2DConstrained(values, classes, X, Y, tables, out_grad, timings);
}

inline bool ComputeDiscreteGradient2DConstrainedTiled(const TileSource2D& source,
        const ByteSource2D& classSource,
        const TileGrid2D& grid,
        const Dgrad2DTables& tables,
        const ByteSink2D& sink,
        const TiledOptions2D& options,
        TileStats2D* stats) {
    return ComputeCLSGradient2DConstrainedTiled(source, classSource, grid, tables, sink, options, stats);
}

inline bool ComputeDiscreteGradient2DTiled(const TileSource2D& source,
        const TileGrid2D& grid,
        const Dgrad2DTables& tables,
        const ByteSink2D& sink,
        const TiledOptions2D& options,
        TileStats2D* stats) {
    return ComputeCLSGradient2DTiled(source, grid, tables, sink, options, stats);
}

inline bool ComputeDiscreteGradient3D(const float* values,
        int64_t X,
        int64_t Y,
        int64_t Z,
        const Dgrad3DTables& tables,
        uint8_t* out_grad,
        Dgrad3DTimings* timings) {
    return ComputeCLSGradient3D(values, X, Y, Z, tables, out_grad, timings);
}

inline bool ComputeDiscreteGradient3DEx(const float* values,
        const uint8_t* classes,
        int keyClamp,
        int64_t X,
        int64_t Y,
        int64_t Z,
        const Dgrad3DTables& tables,
        uint8_t* out_grad,
        const InCoreOptions3D& options,
        Dgrad3DTimings* timings) {
    return ComputeCLSGradient3DEx(values, classes, keyClamp, X, Y, Z, tables, out_grad, options, timings);
}

inline bool ComputeDiscreteGradient3DTiled(const TileSource3D& source,
        const TileGrid3D& grid,
        const Dgrad3DTables& tables,
        const ByteSink3D& sink,
        const TiledOptions3D& options,
        TileStats3D* stats) {
    return ComputeCLSGradient3DTiled(source, grid, tables, sink, options, stats);
}

inline bool ComputeDiscreteGradient3DConstrained(const float* values,
        const uint8_t* classes,
        int64_t X,
        int64_t Y,
        int64_t Z,
        const Dgrad3DTables& tables,
        int keyClamp,
        uint8_t* out_grad,
        Dgrad3DTimings* timings) {
    return ComputeCLSGradient3DConstrained(values, classes, X, Y, Z, tables, keyClamp, out_grad, timings);
}

inline bool ComputeDiscreteGradient3DConstrainedTiled(const TileSource3D& source,
        const ByteSource3D& classSource,
        const TileGrid3D& grid,
        const Dgrad3DTables& tables,
        int keyClamp,
        const ByteSink3D& sink,
        const TiledOptions3D& options,
        TileStats3D* stats) {
    return ComputeCLSGradient3DConstrainedTiled(source, classSource, grid, tables, keyClamp, sink, options, stats);
}

// Earlier form of the PLS / PLS-BR entry points, with the rule as an argument. Kept so
// that existing callers keep compiling; new code calls ComputePLSGradient3D* or
// ComputePLSBRGradient3D*.
enum class DmtLowerStarRule { PLS = 0, PLS_BR = 1 };

inline bool ComputeLowerStarGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
        DmtLowerStarRule rule, uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
        DmtLowerStarStats* stats, Dmt3DTimings* timings) {
    return LowerStarGradient3DGPU(ranks, X, Y, Z, rule == DmtLowerStarRule::PLS_BR,
                                  out_grad, out_step, out_part, stats, timings);
}

inline bool ComputeLowerStarGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
        DmtLowerStarRule rule, const ByteSink3D& gradSink, const DmtTiledOptions3D& options,
        DmtLowerStarStats* stats, TileStats3D* tstats) {
    return LowerStarGradient3DTiled(source, grid, rule == DmtLowerStarRule::PLS_BR, gradSink,
                                    options, stats, tstats);
}

inline bool ComputeLowerStarGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
        DmtLowerStarRule rule, uint8_t* out_grad, const DmtStreamOptions3D& options,
        DmtLowerStarStats* stats, TileStats3D* tstats) {
    return LowerStarGradient3DStreamed(ranks, X, Y, Z, rule == DmtLowerStarRule::PLS_BR, out_grad,
                                       options, stats, tstats);
}

} // namespace gpu
} // namespace GInt

#endif
