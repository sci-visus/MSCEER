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

// Stage 1 A/B variant: Robins lower-star gradient of an X*Y float image
// (row-major, x-fastest) for INTERIOR vertices only (grid coords in
// [1,X-2] x [1,Y-2]), reading the CPU-precomputed per-cell label bytes from
// RegularGridMaxMinVertexLabeling2D ((2X-1)*(2Y-1) bytes each). out_grad
// must be (2X-1)*(2Y-1) bytes, zero-initialized by the caller; only cells
// owned by interior vertices are written. interior_only must be true and
// the label pointers non-null. timings is optional.
bool ComputeDiscreteGradient2D(const float* values, int64_t X, int64_t Y,
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
bool ComputeDiscreteGradient2DBatched(const float* values, int64_t X, int64_t Y,
                                      int64_t nslices,
                                      const Dgrad2DTables& tables,
                                      uint8_t* out_grad, int block_shape,
                                      Dgrad2DTimings* timings);

// Single-image convenience wrapper for the batched full-domain fused path.
bool ComputeDiscreteGradient2DFused(const float* values, int64_t X, int64_t Y,
                                    const Dgrad2DTables& tables,
                                    uint8_t* out_grad,
                                    Dgrad2DTimings* timings);

// Flow-terminal labeling of the 2D discrete gradient (Stage 4 downstream).
// Exploits the offset-doubling structure of the discrete flow: a non-critical
// vertex v is paired with edge e = v + delta and the V-path continues at the
// edge's other endpoint v + 2*delta, so the vertex flow's successor is one
// byte read plus a doubled offset on the vertex lattice - no facet/cofacet
// enumeration, no dimension checks. Symmetrically, a non-critical quad q
// paired with edge e = q + delta is reached (on the reverse path toward its
// maximum) from the opposite cofacet q + 2*delta. Both flows are forests
// rooted at critical cells; terminal labels are computed by in-place pointer
// doubling (log-depth passes, race-tolerant, deterministic result).
//
// grad: the (2X-1)*(2Y-1) GradBitfield bytes. Only the pair bits are read,
// so this works before or after setAscendingManifoldDimensions.
// out_vertex_min: X*Y int32; terminal critical VERTEX grid number per vertex
//   (the base ascending-2-manifold / basin identity of each pixel).
// out_quad_max: (X-1)*(Y-1) int32; terminal critical QUAD lattice index
//   (qx + qy*(X-1)) per quad (the base descending-2-manifold identity).
// Either output may be null to skip that labeling. timings->kernel_ms covers
// init + all doubling passes.
bool Label2DFlowTerminals(const uint8_t* grad, int64_t X, int64_t Y,
                          const Dgrad2DTables& tables,
                          int32_t* out_vertex_min, int32_t* out_quad_max,
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

// ---------------------------------------------------------------------------
// Persistent 2D label context (Stage 6)
// ---------------------------------------------------------------------------

// Every entry point above is alloc -> upload -> compute -> download -> free,
// which is the right shape for a one-shot computation but the wrong one for
// interactive relabeling: the per-pixel base labeling is the large buffer and
// it does not change when the persistence threshold moves. This is the first
// persistent handle - it keeps the base labeling resident on the device so a
// relabel uploads only the (small) remap and downloads only the result.
//
// A LabelCtx2D owns, per direction (ascending / descending), a device copy of
// the base labeling, a device remap buffer, and a device output buffer. All
// three are allocated lazily on first use of that direction and reused across
// calls. The handle is NOT thread-safe; use one per owning object.
struct LabelCtx2D;

// Creates a context for an X*Y image. Returns nullptr if there is no usable
// CUDA device (or on any allocation failure), so the caller can fall back to
// its CPU path without having to unwind a partially built handle.
LabelCtx2D* CreateLabelCtx2D(int64_t X, int64_t Y);

// Frees every device buffer the context owns. Safe on nullptr.
void DestroyLabelCtx2D(LabelCtx2D* ctx);

// Uploads a host COMPACT per-pixel base labeling (X*Y int32, -1 = unlabeled)
// for one direction. Allocates the device buffer on first call; re-uploading
// over an existing labeling is allowed and invalidates any device pointer
// previously returned for that direction.
bool CtxSetBaseLabels(LabelCtx2D* ctx, bool ascending, const int32_t* base_compact);

// out[i] = base[i] < 0 ? -1 : remap[base[i]], with a compact id outside
// [0, m) also yielding -1 (defensive; the CPU reference and any CPU fallback
// must apply the identical rule). remap holds m int32 in ANY caller-defined
// id space - node ids, class ids, colors - and a -1 in remap propagates to
// every pixel of that region. m == 0 is legal and produces an all -1 image.
//
// Uploads remap, runs a gather kernel, and downloads X*Y int32 into
// out_labels_host. Returns false (having written nothing) if the base
// labeling for that direction has not been uploaded, or on any CUDA error.
bool CtxRelabel(LabelCtx2D* ctx, bool ascending, const int32_t* remap, int32_t m,
                int32_t* out_labels_host, Dgrad2DTimings* timings);

// As CtxRelabel but without the download: on success *out_dev_labels is the
// context's internal device int32[X*Y]. The pointer stays valid until the next
// CtxRelabel/CtxRelabelDevice/CtxSetBaseLabels on the SAME direction, or until
// the context is destroyed. The caller must not free it.
bool CtxRelabelDevice(LabelCtx2D* ctx, bool ascending, const int32_t* remap,
                      int32_t m, const void** out_dev_labels,
                      Dgrad2DTimings* timings);

// The device int32[X*Y] base labeling for a direction, or nullptr if it has
// not been uploaded. Valid until the next CtxSetBaseLabels on that direction
// or until destroy. The caller must not free it.
const void* CtxBaseLabelsDevice(const LabelCtx2D* ctx, bool ascending);

// Free / total device memory in bytes for the current device. Either output
// may be null. Returns false if there is no usable device. Exists so leak
// gates can be written in a plain host TU without including <cuda_runtime.h>.
bool QueryDeviceMemory(int64_t* free_bytes, int64_t* total_bytes);

} // namespace gpu
} // namespace GInt

#endif
