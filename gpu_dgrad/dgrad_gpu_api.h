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
// (SlidingWindowRobinsNoalloc::ComputePairing_sliding) - the algorithm is
// comparisons-only under the value-then-index total order, so any mismatch is
// a bug, not floating-point noise.

#include <cstdint>

namespace GInt {
namespace gpu {

struct DeviceInfo {
    char name[256];
    int sm_count;
    long long total_mem_bytes;
    int cc_major;
    int cc_minor;
};

// Fills 'out' with the properties of CUDA device 0.
// Returns false if no usable CUDA device is present.
bool QueryDevice(DeviceInfo& out);

// ---------------------------------------------------------------------------
// 2D path (Stage 1)
// ---------------------------------------------------------------------------

// Local topology of the 3x3 cell neighborhood of a vertex in the doubled grid.
// Slot s = (dx+1) + 3*(dy+1) for dx,dy in {-1,0,1} - the same enumeration as
// TopologicalRegularGrid2D::m_adjacent_cell_offsets, so slot 4 is the vertex.
// The tables MUST be generated from the real mesh's iterators (see the
// validation harness) so the device kernel provably replicates the CPU
// iteration orders; they are never hand-typed.
struct Dgrad2DTables {
    uint8_t dim[9];                            // cell dimension per slot
    uint8_t fac_count[9]; uint8_t fac[9][2];   // facets inside the 3x3, FacetsIterator order
    uint8_t cof_count[9]; uint8_t cof[9][4];   // cofacets inside the 3x3, CofacetsIterator order
    uint8_t dir_px, dir_mx, dir_py, dir_my;    // GradBitfield pair codes for +x,-x,+y,-y
};

struct Dgrad2DTimings {
    float h2d_ms;
    float kernel_ms;
    float d2h_ms;
};

// Computes the Robins lower-star discrete gradient of an X*Y float image
// (row-major, x-fastest) for INTERIOR vertices only (grid coords in
// [1,X-2] x [1,Y-2]). max_labels/min_labels are the CPU-precomputed
// per-cell bytes from RegularGridMaxMinVertexLabeling2D ((2X-1)*(2Y-1)
// bytes each). out_grad must be (2X-1)*(2Y-1) bytes, zero-initialized by
// the caller; cells owned by interior vertices are written with the exact
// GradBitfield byte the CPU path produces, all other cells are untouched.
// interior_only must be true in Stage 1. timings is optional.
bool ComputeDiscreteGradient2D(const float* values, int64_t X, int64_t Y,
                               const uint8_t* max_labels,
                               const uint8_t* min_labels,
                               const Dgrad2DTables& tables,
                               uint8_t* out_grad,
                               bool interior_only,
                               Dgrad2DTimings* timings);

// Computes the discrete gradient of a scalar grid of X*Y*Z floats
// (row-major, x-fastest) into out_grad, sized (2X-1)*(2Y-1)*(2Z-1) bytes.
//
// interior_only: fill only cells owned by interior vertices (Stage 1);
// Stage 3 removes the restriction via boundary validity masking.
//
// Stage 1 additionally takes the CPU-precomputed max/min vertex label arrays;
// Stage 2 (label fusion) drops them, at which point this signature loses the
// label parameters. max_labels/min_labels: 1 byte per cell as produced by
// RegularGridMaxMinVertexLabeling3D, or nullptr once fusion lands.
//
// Returns false on failure (no device, allocation failure, or - until
// Stage 1 lands - not yet implemented).
bool ComputeDiscreteGradient3D(const float* values,
                               int64_t X, int64_t Y, int64_t Z,
                               const uint8_t* max_labels,
                               const uint8_t* min_labels,
                               uint8_t* out_grad,
                               bool interior_only);

} // namespace gpu
} // namespace GInt

#endif
