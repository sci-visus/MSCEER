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

#include "tile2d.h"
#include "tile3d.h"

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

// Stage A5: CONSTRAINED fused pairing. classes holds one restriction class
// byte per doubled-lattice cell ((2X-1)*(2Y-1)); each lower star is split
// into subsets keyed by class*maxDim() + boundaryValue and every subset is
// homotopy-expanded independently in ascending key order, exactly as the
// 4-argument MyRobinsNoalloc(mesh, maxv, restriction, grad) does
// (gi_modified_robins.h:1281-1340, collisions of that key included), so the
// output memcmp-matches the CPU's constrained pairing. An all-zero class
// array reproduces ComputeDiscreteGradient2DFused byte for byte.
bool ComputeDiscreteGradient2DConstrained(const float* values, const uint8_t* classes,
                                          int64_t X, int64_t Y,
                                          const Dgrad2DTables& tables,
                                          uint8_t* out_grad,
                                          Dgrad2DTimings* timings);

// Tiled variant: the class bytes come from a source callback (cell boxes).
bool ComputeDiscreteGradient2DConstrainedTiled(const TileSource2D& source,
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
// ComputeDiscreteGradient2DFused for every tile size and residency (the same
// kernel runs on "views"; see dgrad2d_kernel.cu). grid comes from
// MakeTileGrid2D(X, Y, T); T >= 1. stats is optional.
bool ComputeDiscreteGradient2DTiled(const TileSource2D& source,
                                    const TileGrid2D& grid,
                                    const Dgrad2DTables& tables,
                                    const ByteSink2D& sink,
                                    const TiledOptions2D& options,
                                    TileStats2D* stats);

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

// Stage A1: TILED flow-terminal labeling - the first client of the tiled
// SuccessorResolver (successor2d.cu). Per tile the successor field is
// resolved by in-tile pointer doubling with successors clamped at the tile
// box (a fixed log-depth number of passes, no host flag); cells whose flow
// leaves the tile become frontier records indexed by ring position. Phase 1
// streams the gradient and downloads only perimeter-sized ring data, the
// host resolves the frontier graph, phase 2 streams again and substitutes
// the resolved payloads on the way to the sink. Nothing proportional to the
// image is ever held. ascending=true labels the vertex lattice (X*Y words:
// terminal minimum), false the quad lattice ((X-1)*(Y-1) words: terminal
// maximum). In lattice-index mode the output is byte-identical to
// Label2DFlowTerminals; see ResolveOptions2D for compact ids.
bool Label2DFlowTerminalsTiled(const ByteSource2D& gradSource,
                               const TileGrid2D& grid,
                               const Dgrad2DTables& tables, bool ascending,
                               const WordSink2D& sink,
                               const ResolveOptions2D& options,
                               ResolveStats2D* stats);

// Stage A2: TILED numeric region labeling of the vertex lattice (basins of
// minima at a persistence): per owned vertex, decide (base minimum ->
// remap[base]; f(x) < sStar[remap[base(x)]] -> certain) or integrate a
// numeric hop (numeric2d.cuh) whose snapped successor feeds the same tiled
// resolver as Label2DFlowTerminalsTiled. Inputs per tile come from
// callbacks: function values (halo 2), gradient bytes (owned cell box),
// compact base labels (halo 1). remap[nBase] -> living id in [0, nLiving);
// sStar[nLiving] = lowest boundary saddle value of each living region (the
// certain-core threshold). Output: living id per vertex (X*Y words).
bool LabelNumericRegions2DTiled(const TileSource2D& values,
                                const ByteSource2D& gradSource,
                                const WordSource2D& baseLabels,
                                const TileGrid2D& grid,
                                const Dgrad2DTables& tables,
                                const int32_t* remap, int32_t nBase,
                                const float* sStar, int32_t nLiving,
                                const WordSink2D& sink,
                                const NumericOptions2D& options,
                                NumericStats2D* stats);

// Stage A3: TILED consistency pass on a vertex region labeling. Successor of
// a vertex = its lowest axis neighbour (under the (value, index) order) that
// carries the same region label, else its lowest axis neighbour overall; a
// vertex with no lower axis neighbour is a base minimum and a terminal with
// its own label. The output label is that of the chain's terminal, so every
// non-minimum vertex ends up with a lower AXIS neighbour in its own region -
// the property the constrained lower-star pairing needs to never manufacture
// a boundary minimum (the CPU's IsolatedCCRegionRemover role). Successors
// are +-1 steps, so the result is independent of the tile size. Inputs:
// values and region words (halo 1 each); options.compactIds is ignored.
bool ConsistencyPass2DTiled(const TileSource2D& values,
                            const WordSource2D& regions,
                            const TileGrid2D& grid,
                            const WordSink2D& sink,
                            const ResolveOptions2D& options,
                            ResolveStats2D* stats);

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

// ---------------------------------------------------------------------------
// 3D lower-star (Robins) gradient (dgrad3d_*.cu)
// ---------------------------------------------------------------------------

// Local topology of the 27 cells around a vertex of the doubled grid. Slot
// s = (dx+1) + 3*(dy+1) + 9*(dz+1), the AdjacentCellsIterator order of both
// TopologicalRegularGrid3D and Explicit3x3x3SmallRegularGrid, so slot 13 is
// the vertex. Vertex-window positions w use the same encoding for the grid
// offset of a vertex (the small RegularGrid3D index of the CPU sliding
// window), so w orders the 27 window vertices exactly as the global vertex
// index does. Generated from the real iterators (dgrad3d_tables.h), never
// hand-typed.
struct Dgrad3DTables {
    uint8_t dim[27];
    uint8_t fac_count[27]; uint8_t fac[27][6];   // facets inside the window, FacetsIterator order
    uint8_t cof_count[27]; uint8_t cof[27][6];   // cofacets inside the window, CofacetsIterator order
    uint8_t vert_count[27]; uint8_t vert[27][8]; // vertices of the cell as window positions
};

struct Dgrad3DTimings {
    float h2d_ms;      // summed per-batch event windows
    float kernel_ms;
    float d2h_ms;
    float wall_ms;     // end to end, host copies included
};

// Lower-star gradient of an X*Y*Z float volume (x fastest), memcmp-equal to
// SlidingWindowRobinsNoalloc::ComputePairing_sliding (== MyRobinsNoalloc::
// ComputePairing): one byte per cell of the (2X-1)(2Y-1)(2Z-1) doubled grid,
// every byte written. Membership and lowest-vertex queries are computed on the
// fly (no max/min label arrays), and each lower star is split by
// boundaryValue(cell) exactly as the CPU does (non-periodic grids). This
// in-core entry point is the tiled path with a single tile. Returns false on
// any CUDA failure so the caller can fall back to the CPU.
bool ComputeDiscreteGradient3D(const float* values, int64_t X, int64_t Y, int64_t Z,
                               const Dgrad3DTables& tables, uint8_t* out_grad,
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
// MyRobinsNoalloc<...,4,6> of ondemandaccurate), collisions of that key
// included, as the 4-argument MyRobinsNoalloc does. All-zero classes reproduce
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

// ---------------------------------------------------------------------------
// 3D discrete Morse theory rules (dmt3d_gradient.cu)
// ---------------------------------------------------------------------------

// GPU versions of the raw host drivers in include/gi_dmt_lowerstar_gradient.h
// (ComputeLowerStarGradient3D) and include/gi_dmt_dmf_gradient.h
// (ComputeDMFGradient3D). The device code runs the per-star / per-cell
// functions of include/dmt_core unchanged, so every output array is
// memcmp-identical to the CPU driver's (dmt3d_validate checks this). The
// enums and stats structs mirror the CPU ones field for field; they are
// redeclared here because the gi_* headers must not reach nvcc.
//
// ranks: X*Y*Z uint32, x fastest, pairwise distinct (only their order
// matters). Every cell output is (2X-1)(2Y-1)(2Z-1) bytes (or int64 for
// out_order), x fastest, in host memory. Optional outputs and stats/timings
// may be nullptr. Returns false on any CUDA failure (no device, allocation,
// launch), in which case the caller should run the CPU driver instead.

enum class DmtLowerStarRule { PLS = 0, PLS_BR = 1 };

struct DmtLowerStarStats {
    int64_t critical[4];
    int64_t orphans;
    int64_t unassigned;   // cells no star wrote; 0 on success
};

struct DmtDMFStats {
    int64_t guard1, guard2, alg1, alg2;
    int64_t critical[4];
};

struct Dmt3DTimings {
    float h2d_ms;
    float kernel_ms;   // every kernel of the call, including stats reductions
    float d2h_ms;
};

// PLS / PLS-BR, one thread per vertex. out_grad is required; out_step and
// out_part are the optional per-cell journal arrays. Returns false on a CUDA
// failure; unlike the CPU driver it returns true when cells were left
// unassigned (stats->unassigned reports them, 0 on every valid input).
bool ComputeLowerStarGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                                   DmtLowerStarRule rule, uint8_t* out_grad,
                                   uint8_t* out_step, uint8_t* out_part,
                                   DmtLowerStarStats* stats, Dmt3DTimings* timings);

// Shivashankar-Natarajan DMF: the top-facet precompute, then one kernel per
// algorithm per dimension d = 0..dmax-1 (Algorithm 1, then Algorithm 2 if
// second_pass), guard and claim counts by block reductions. out_stage: 1/2 =
// paired by Algorithm 1/2, 0 = critical. out_order: each cell's position in
// the (dimension, key) order, by per-dimension LSD radix sorts on the device.
bool ComputeDMFGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                             bool second_pass, uint8_t* out_grad, uint8_t* out_stage,
                             int64_t* out_order, DmtDMFStats* stats, Dmt3DTimings* timings);

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
