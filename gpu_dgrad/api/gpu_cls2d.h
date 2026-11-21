#ifndef GPU_CLS2D_H
#define GPU_CLS2D_H

// 2D lower-star gradient on the GPU (cls/cls2d_*.cu). Plain C++.
//
// Output contract: one byte per cell of the doubled grid (2X-1) x (2Y-1), row-major
// x-fastest, with the bit layout of GradBitfield in include/gi_discrete_gradient_labeling.h.
// The output memcmp-matches the CPU path.

#include <cstdint>
#include "gpu_tables.h"
#include "tile2d.h"

namespace GInt {
namespace gpu {

// Stage 1 A/B variant: Robins lower-star gradient of an X*Y float image
// (row-major, x-fastest) for INTERIOR vertices only (grid coords in
// [1,X-2] x [1,Y-2]), reading the CPU-precomputed per-cell label bytes from
// RegularGridMaxMinVertexLabeling2D ((2X-1)*(2Y-1) bytes each). out_grad
// must be (2X-1)*(2Y-1) bytes, zero-initialized by the caller; only cells
// owned by interior vertices are written. interior_only must be true and
// the label pointers non-null. timings is optional.
bool ComputeCLSGradient2D(const float* values, int64_t X, int64_t Y,
                               const uint8_t* max_labels,
                               const uint8_t* min_labels,
                               const Dgrad2DTables& tables,
                               uint8_t* out_grad,
                               bool interior_only,
                               Dgrad2DTimings* timings);

// Stage 3 production path: FUSED (no label arrays - membership and
// lowest-vertex queries computed on the fly from a shared-memory tile) and
// FULL-DOMAIN (boundary vertices handled in-kernel by expanding each
// boundary-class subset independently, matching the CPU's boundaryValue
// stratification - non-periodic grids). Batched: values holds nslices
// independent X*Y images back to back; out_grad receives nslices gradients
// of (2X-1)*(2Y-1) bytes each, every byte written, memcmp-equal to the CPU
// path per slice. block_shape: 0 = 16x16 (default), 1 = 32x8, 2 = 8x8.
bool ComputeCLSGradient2DBatched(const float* values, int64_t X, int64_t Y,
                                      int64_t nslices,
                                      const Dgrad2DTables& tables,
                                      uint8_t* out_grad, int block_shape,
                                      Dgrad2DTimings* timings);

// Single-image convenience wrapper for the batched full-domain fused path.
bool ComputeCLSGradient2DFused(const float* values, int64_t X, int64_t Y,
                                    const Dgrad2DTables& tables,
                                    uint8_t* out_grad,
                                    Dgrad2DTimings* timings);

// Stage A5: CONSTRAINED fused pairing. classes holds one restriction class
// byte per doubled-lattice cell ((2X-1)*(2Y-1)); each lower star is split
// into subsets keyed by class*maxDim() + boundaryValue and every subset is
// homotopy-expanded independently in ascending key order, exactly as the
// 4-argument CLSGradient(mesh, maxv, restriction, grad) does
// (gradient/gi_cls_gradient.h, collisions of that key included), so the
// output memcmp-matches the CPU's constrained pairing. An all-zero class
// array reproduces ComputeCLSGradient2DFused byte for byte.
bool ComputeCLSGradient2DConstrained(const float* values, const uint8_t* classes,
                                          int64_t X, int64_t Y,
                                          const Dgrad2DTables& tables,
                                          uint8_t* out_grad,
                                          Dgrad2DTimings* timings);

// Tiled variant: the class bytes come from a source callback (cell boxes).
bool ComputeCLSGradient2DConstrainedTiled(const TileSource2D& source,
                                               const ByteSource2D& classSource,
                                               const TileGrid2D& grid,
                                               const Dgrad2DTables& tables,
                                               const ByteSink2D& sink,
                                               const TiledOptions2D& options,
                                               TileStats2D* stats);

// Stage A0: TILED full-domain fused path (tile2d.h). Streams the image tile
// by tile through a ring of pinned + device buffers on three CUDA streams, so
// neither the image nor the gradient has to exist in host or device memory
// as a whole: the source callback supplies a box of values, the sink receives
// a box of gradient bytes. Output is byte-identical to
// ComputeCLSGradient2DFused for every tile size and residency (the same
// kernel runs on "views"; see cls2d_kernel.cu). grid comes from
// MakeTileGrid2D(X, Y, T); T >= 1. stats is optional.
bool ComputeCLSGradient2DTiled(const TileSource2D& source,
                                    const TileGrid2D& grid,
                                    const Dgrad2DTables& tables,
                                    const ByteSink2D& sink,
                                    const TiledOptions2D& options,
                                    TileStats2D* stats);

} // namespace gpu
} // namespace GInt

#endif
