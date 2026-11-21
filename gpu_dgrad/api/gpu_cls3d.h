#ifndef GPU_CLS3D_H
#define GPU_CLS3D_H

// 3D lower-star gradient on the GPU (cls/cls3d_*.cu). Plain C++.
//
// Output contract: one byte per cell of the doubled grid (2X-1)(2Y-1)(2Z-1), x fastest,
// with the bit layout of GradBitfield in include/gi_discrete_gradient_labeling.h. The
// output memcmp-matches the CPU path.

#include <cstdint>
#include "gpu_device.h"
#include "gpu_tables.h"
#include "tile3d.h"

namespace GInt {
namespace gpu {

// Lower-star gradient of an X*Y*Z float volume (x fastest), memcmp-equal to
// CLSSlidingGradient::ComputePairing_sliding (== CLSGradient::
// ComputePairing): one byte per cell of the (2X-1)(2Y-1)(2Z-1) doubled grid,
// every byte written. Membership and lowest-vertex queries are computed on the
// fly (no max/min label arrays), and each lower star is split by
// boundaryValue(cell) exactly as the CPU does (non-periodic grids). This
// in-core entry point is ComputeDiscreteGradient3DEx with default options
// (direct z-slab DMA when the arrays are page-locked, else z-slabs through
// the pinned ring). Returns false on any CUDA failure so the caller can fall
// back to the CPU.
bool ComputeDiscreteGradient3D(const float* values, int64_t X, int64_t Y, int64_t Z,
                               const Dgrad3DTables& tables, uint8_t* out_grad,
                               Dgrad3DTimings* timings);

// In-core options. Both host modes stream the volume in z-slabs (full XY
// planes, slabDepth vertices each), so the transfers overlap the kernel:
//   copy   - slabs are staged through pinned buffers by host-thread copies
//            (any host memory);
//   direct - the slabs are contiguous in the caller's arrays, so they are
//            DMAed straight between those arrays and the device, with no host
//            copy and no host wait per slab. Needs page-locked arrays:
//            cudaMallocHost memory or PinHostMemory below.
// hostMode 0 (auto) = direct when every array is page-locked, else copy;
// 1 = copy; 2 = direct, page-locking the arrays for the duration of the call
// (registering is cheap, unregistering costs ~0.7 s per 12 GB on this
// machine, so pin once and reuse when the call repeats).
struct InCoreOptions3D {
    int hostMode = 0;
    int slabDepth = 0;   // vertices per slab; 0 = about 64 MB of cells per slab
    int slots = 3;       // slabs in flight
    bool packOutput = false;   // copy mode: packed download (TiledOptions3D::packOutput)
};

bool ComputeDiscreteGradient3DEx(const float* values, const uint8_t* classes, int keyClamp,
                                 int64_t X, int64_t Y, int64_t Z, const Dgrad3DTables& tables,
                                 uint8_t* out_grad, const InCoreOptions3D& options,
                                 Dgrad3DTimings* timings);

// Tiled variant (tile3d.h): values and gradient move tile by tile through a
// ring of pinned + device buffers on three streams (H2D, kernel, D2H), so
// neither array has to fit in device or host memory. Byte-identical to
// ComputeDiscreteGradient3D for every tile size, residency and batching.
bool ComputeDiscreteGradient3DTiled(const TileSource3D& source, const TileGrid3D& grid,
                                    const Dgrad3DTables& tables, const ByteSink3D& sink,
                                    const TiledOptions3D& options, TileStats3D* stats);

// Constrained pairing: classes holds one restriction byte per cell; lower-star
// subsets are keyed by (char)(class*3 + boundaryValue), clamped to
// keyClamp (= LMAXDIM*LMAX_LABEL-1 of the CPU instantiation, 23 for the
// CLSGradient<...,4,6> of ondemandaccurate), collisions of that key
// included, as the 4-argument CLSGradient does. All-zero classes reproduce
// ComputeDiscreteGradient3D byte for byte.
bool ComputeDiscreteGradient3DConstrained(const float* values, const uint8_t* classes,
                                          int64_t X, int64_t Y, int64_t Z,
                                          const Dgrad3DTables& tables, int keyClamp,
                                          uint8_t* out_grad, Dgrad3DTimings* timings);

bool ComputeDiscreteGradient3DConstrainedTiled(const TileSource3D& source,
                                               const ByteSource3D& classSource,
                                               const TileGrid3D& grid,
                                               const Dgrad3DTables& tables, int keyClamp,
                                               const ByteSink3D& sink,
                                               const TiledOptions3D& options,
                                               TileStats3D* stats);

} // namespace gpu
} // namespace GInt

#endif
