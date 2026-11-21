#ifndef GPU_DMF3D_H
#define GPU_DMF3D_H

// Shivashankar-Natarajan DMF on the GPU (dmt3d/). Plain C++. Output is memcmp-identical
// to the CPU driver ComputeDMFGradient3D.

#include <cstdint>
#include "gpu_dmt3d.h"

namespace GInt {
namespace gpu {

struct DmtDMFStats {
    int64_t guard1, guard2, alg1, alg2;
    int64_t critical[4];
};

// Shivashankar-Natarajan DMF: the top-facet precompute, then one kernel per
// algorithm per dimension d = 0..dmax-1 (Algorithm 1, then Algorithm 2 if
// second_pass), guard and claim counts by block reductions. out_stage: 1/2 =
// paired by Algorithm 1/2, 0 = critical. out_order: each cell's position in
// the (dimension, key) order, by per-dimension LSD radix sorts on the device.
bool ComputeDMFGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                             bool second_pass, uint8_t* out_grad, uint8_t* out_stage,
                             int64_t* out_order, DmtDMFStats* stats, Dmt3DTimings* timings);

// Tiled version (see gpu_dmt3d.h). stageSink may be null.
bool ComputeDMFGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
                               bool second_pass, const ByteSink3D& gradSink,
                               const ByteSink3D* stageSink, const DmtTiledOptions3D& options,
                               DmtDMFStats* stats, TileStats3D* tstats);

// Streamed version (see DmtStreamOptions3D). out_stage may be null.
bool ComputeDMFGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z, bool second_pass,
                                  uint8_t* out_grad, uint8_t* out_stage, const DmtStreamOptions3D& options,
                                  DmtDMFStats* stats, TileStats3D* tstats);

} // namespace gpu
} // namespace GInt

#endif
