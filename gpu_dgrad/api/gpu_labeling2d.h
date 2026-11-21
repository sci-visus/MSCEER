#ifndef GPU_LABELING2D_H
#define GPU_LABELING2D_H

// 2D labeling on the GPU (labeling2d/): flow terminals, numeric regions, the consistency
// pass and the persistent label context. Plain C++.

#include <cstdint>
#include "gpu_tables.h"
#include "tile2d.h"

namespace GInt {
namespace gpu {

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

// ---------------------------------------------------------------------------
// Persistent 2D label context
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

} // namespace gpu
} // namespace GInt

#endif
