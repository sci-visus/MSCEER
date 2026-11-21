#ifndef GPU_DMT3D_H
#define GPU_DMT3D_H

// Types shared by the PLS, PLS-BR and DMF GPU entry points (dmt3d/). Plain C++.

#include <cstdint>
#include "tile3d.h"

namespace GInt {
namespace gpu {

// GPU versions of the raw host drivers ComputePLSGradient3D, ComputePLSBRGradient3D
// and ComputeDMFGradient3D (include/gradient/). The device code runs the per-star / per-cell
// functions of include/dmt_core unchanged, so every output array is
// memcmp-identical to the CPU driver's (dmt3d_validate checks this). The
// stats structs mirror the CPU ones field for field; they are
// redeclared here because the gi_* headers must not reach nvcc.
//
// ranks: X*Y*Z uint32, x fastest, pairwise distinct (only their order
// matters). Every cell output is (2X-1)(2Y-1)(2Z-1) bytes (or int64 for
// out_order), x fastest, in host memory. Optional outputs and stats/timings
// may be nullptr. Returns false on any CUDA failure (no device, allocation,
// launch), in which case the caller should run the CPU driver instead.

struct DmtLowerStarStats {
    int64_t critical[4];
    int64_t orphans;
    int64_t unassigned;   // cells no star wrote; 0 on success
};

struct Dmt3DTimings {
    float h2d_ms;
    float kernel_ms;   // every kernel of the call, including stats reductions
    float d2h_ms;
};

// Tiled (blocked) versions with 32-bit tile-local indexing (dmt3d_tiled.cu).
// Ownership as in tile3d.h (cell c owned by the tile owning vertex floor(c/2));
// each tile runs every pass of its rule on its staged box (owned + halo, see
// DmtTiledOptions3D) as a cubical complex of its own, and only owned cells
// leave the device. Outputs and stats (counted on owned cells only) are
// byte-identical to the in-core drivers above for every tile shape.
// The DMF (dimension, key) order output is global and is not available tiled.

// Streamed DMT gradients of a volume in host memory (any size: only slabs live
// on the device). hostMode 0 (auto): when ranks and outputs are page-locked
// (cudaMallocHost or PinHostMemory), full-XY z-slabs of slabDepth vertex
// planes are DMAed straight between those arrays and the device, and run as
// (tileX, tileY) tiles inside each slab; otherwise the tiled ring driver
// copies through pinned buffers. 1 = always copy, 2 = pin for the call.
// Output and stats are identical to the tiled and in-core drivers.
struct DmtStreamOptions3D {
    int hostMode = 0;
    int tileX = 128, tileY = 128;   // tiles inside a slab (DMF is fastest on small tiles)
    int slabDepth = 64;             // vertex planes per slab (= tile depth)
    int slots = 2;                  // slabs in flight (direct); ring slots (copy)
};

// Shared implementation of PLS (boundary_restricted = false) and PLS-BR (true): one
// kernel with a flag. Call the named entry points of gpu_pls3d.h and gpu_pls_br3d.h.
bool LowerStarGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                            bool boundary_restricted, uint8_t* out_grad,
                            uint8_t* out_step, uint8_t* out_part,
                            DmtLowerStarStats* stats, Dmt3DTimings* timings);
bool LowerStarGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
                              bool boundary_restricted, const ByteSink3D& gradSink,
                              const DmtTiledOptions3D& options, DmtLowerStarStats* stats,
                              TileStats3D* tstats);
bool LowerStarGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                                 bool boundary_restricted, uint8_t* out_grad,
                                 const DmtStreamOptions3D& options, DmtLowerStarStats* stats,
                                 TileStats3D* tstats);

} // namespace gpu
} // namespace GInt

#endif
