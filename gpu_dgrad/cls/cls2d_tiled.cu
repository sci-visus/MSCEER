// Stage A0: tiled driver for the 2D fused Robins kernel.
//
// tile = ownership unit (tile2d.h). Tiles are visited in row-major order and
// packed 'tilesPerLaunch' at a time into a TransferRing slot: the staged value
// boxes go back to back in the slot's pinned input buffer, the owned cell
// boxes back to back in the output buffer, and one launch with blockIdx.z
// over the batch computes them. The host loop overlaps source reads and sink
// writes with the device's H2D / kernel / D2H of other slots.
//
// Nothing here is allocated in proportion to the image: 'slots' batches of
// 'tilesPerLaunch' tiles is the whole footprint, host and device.

#include "cls2d_internal.cuh"
#include "common/gpu_transfer.h"
#include "tile2d.h"

#include <chrono>
#include <cstdio>
#include <vector>

namespace GInt {
namespace gpu {

namespace {

struct BatchPlan {
    int64_t firstTile = 0;
    int count = 0;
    size_t inBytes = 0, outBytes = 0;
    int maxCw = 0, maxCh = 0;
    std::vector<DeviceTileView2D> views;   // device pointers filled at enqueue
    std::vector<Box2D> staged, owned;
};

#define TILED_CUDA_CHECK(call)                                                  \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            return false;                                                       \
        }                                                                       \
    } while (0)

} // namespace

// Shared driver; classSource (Stage A5) is optional: when given, each tile
// also stages the restriction class bytes of the cells around its computed
// vertices and the kernel groups lower-star subsets by class*maxDim+bv.
static bool TiledImpl(const TileSource2D& source, const ByteSource2D* classSource,
                      const TileGrid2D& grid, const Dgrad2DTables& tables,
                      const ByteSink2D& sink, const TiledOptions2D& options,
                      TileStats2D* stats) {
    if (grid.X < 3 || grid.Y < 3 || grid.T < 1 || grid.numTiles < 1) return false;
    int tpl = options.tilesPerLaunch < 1 ? 1
            : options.tilesPerLaunch > 65535 ? 65535 : options.tilesPerLaunch;
    if ((int64_t)tpl > grid.numTiles) tpl = (int)grid.numTiles;
    int nslots = options.slots < 1 ? 1 : options.slots;
    const int64_t batches = (grid.numTiles + tpl - 1) / tpl;
    if ((int64_t)nslots > batches) nslots = (int)batches;   // residency never exceeds the tile count

    const auto t_wall0 = std::chrono::steady_clock::now();
    TileStats2D st;
    st.tiles = grid.numTiles;

    // capacities: the largest staged box is (T+4)^2 values (+ (2T+8)^2 class
    // bytes), the largest owned cell box (2T)^2 bytes, times the tiles per launch
    const size_t maxStaged = (size_t)(grid.T + 4) * (size_t)(grid.T + 4);
    const size_t maxCls = classSource ? (((size_t)(2 * grid.T + 8) * (size_t)(2 * grid.T + 8) + 3) & ~(size_t)3) : 0;
    const size_t maxOwned = (size_t)(2 * grid.T) * (size_t)(2 * grid.T);
    const size_t capIn = (maxStaged * sizeof(float) + maxCls) * (size_t)tpl;
    const size_t capOut = maxOwned * (size_t)tpl;

    TransferRing ring;
    if (!ring.Create(nslots, capIn, capOut)) return false;
    DeviceTileView2D* d_views = nullptr;
    TILED_CUDA_CHECK(cudaMalloc(&d_views, (size_t)nslots * tpl * sizeof(DeviceTileView2D)));
    if (!UploadTables2D(tables)) { cudaFree(d_views); return false; }
    st.hostPinnedBytes = (int64_t)ring.PinnedBytes();
    st.deviceBytesAllocated = (int64_t)ring.DeviceBytes() +
                              (int64_t)((size_t)nslots * tpl * sizeof(DeviceTileView2D));

    std::vector<BatchPlan> plans((size_t)nslots);

    // drain one slot: wait for its D2H, hand the owned boxes to the sink
    auto retire = [&](int slot) -> bool {
        BatchPlan& p = plans[(size_t)slot];
        if (p.count == 0) return true;
        float h2d = 0.f, kern = 0.f, d2h = 0.f;
        if (!ring.Retire(slot, &h2d, &kern, &d2h)) return false;
        st.h2d_ms += h2d;
        st.kernel_ms += kern;
        st.d2h_ms += d2h;
        const unsigned char* out = (const unsigned char*)ring.Slot(slot).h_out;
        size_t off = 0;
        for (int i = 0; i < p.count; i++) {
            const Box2D& o = p.owned[(size_t)i];
            if (!sink(o.x0, o.y0, o.w, o.h, out + off, o.w)) return false;
            off += (size_t)o.w * o.h;
        }
        p.count = 0;
        return true;
    };

    bool ok = true;
    int64_t next = 0;
    int slot = 0;
    while (next < grid.numTiles && ok) {
        if (!retire(slot)) { ok = false; break; }
        BatchPlan& p = plans[(size_t)slot];
        p.firstTile = next;
        p.count = (int)((grid.numTiles - next) < tpl ? (grid.numTiles - next) : tpl);
        p.views.resize((size_t)p.count);
        p.staged.resize((size_t)p.count);
        p.owned.resize((size_t)p.count);
        p.inBytes = p.outBytes = 0;
        p.maxCw = p.maxCh = 0;

        TransferSlot& s = ring.Slot(slot);
        unsigned char* h_in = (unsigned char*)s.h_in;
        const unsigned char* d_in = (const unsigned char*)s.d_in;
        size_t inOff = 0, outOff = 0;   // bytes
        for (int i = 0; i < p.count && ok; i++) {
            const TileDesc2D d = TileAt(grid, next + i);
            const Box2D sb = StagedVertexBox2D(grid, d);
            const Box2D cb = ComputedVertexBox2D(grid, d);
            const Box2D ob = OwnedCellBox2D(grid, d);
            p.staged[(size_t)i] = sb;
            p.owned[(size_t)i] = ob;
            if (!source(sb.x0, sb.y0, sb.w, sb.h, (float*)(h_in + inOff), sb.w)) { ok = false; break; }

            DeviceTileView2D& v = p.views[(size_t)i];
            v.vals = (const float*)(d_in + inOff);
            inOff += (size_t)sb.w * sb.h * sizeof(float);
            v.cls = nullptr;
            v.cpitch = 0; v.ccx0 = v.ccy0 = v.ccw = v.cch = 0;
            if (classSource) {
                // cells around the computed vertices: the last computed vertex
                // v = cx0+cw-1 owns cells up to 2v+1, so [2cx0-1, 2(cx0+cw)) clipped
                const int64_t x0 = 2 * cb.x0 - 1 < 0 ? 0 : 2 * cb.x0 - 1;
                const int64_t y0 = 2 * cb.y0 - 1 < 0 ? 0 : 2 * cb.y0 - 1;
                const int64_t x1 = 2 * (cb.x0 + cb.w) > 2 * grid.X - 1 ? 2 * grid.X - 1 : 2 * (cb.x0 + cb.w);
                const int64_t y1 = 2 * (cb.y0 + cb.h) > 2 * grid.Y - 1 ? 2 * grid.Y - 1 : 2 * (cb.y0 + cb.h);
                const int w = (int)(x1 - x0), h = (int)(y1 - y0);
                if (!(*classSource)(x0, y0, w, h, h_in + inOff, w)) { ok = false; break; }
                v.cls = d_in + inOff;
                v.cpitch = w; v.ccx0 = (int)x0; v.ccy0 = (int)y0; v.ccw = w; v.cch = h;
                inOff += ((size_t)w * h + 3) & ~(size_t)3;
            }
            v.pitch = sb.w;
            v.sx0 = (int)sb.x0; v.sy0 = (int)sb.y0; v.svx = sb.w; v.svy = sb.h;
            v.cx0 = (int)cb.x0; v.cy0 = (int)cb.y0; v.cw = cb.w; v.ch = cb.h;
            v.ox = (int)d.ox; v.oy = (int)d.oy; v.tw = d.tw; v.th = d.th;
            v.X = (int)grid.X; v.Y = (int)grid.Y;
            v.grad = (unsigned char*)s.d_out + outOff;
            v.gpitch = ob.w;
            v.gcx0 = (int)ob.x0; v.gcy0 = (int)ob.y0;
            if (cb.w > p.maxCw) p.maxCw = cb.w;
            if (cb.h > p.maxCh) p.maxCh = cb.h;
            outOff += (size_t)ob.w * ob.h;
        }
        if (!ok) break;
        p.inBytes = inOff;
        p.outBytes = outOff;

        // views for this slot live in their own region of d_views. They go
        // onto the H2D stream BEFORE the values, so the slot's ev_h2d (which
        // the compute stream waits on) covers both copies by stream order.
        // The previous batch of this slot has been retired above, so neither
        // its kernel nor its copies can still touch this region. The views
        // array is pageable; that copy is a few hundred bytes.
        DeviceTileView2D* dv = d_views + (size_t)slot * tpl;
        cudaError_t e = cudaMemcpyAsync(dv, p.views.data(),
                                        (size_t)p.count * sizeof(DeviceTileView2D),
                                        cudaMemcpyHostToDevice, ring.StreamH2D());
        if (e != cudaSuccess) { ok = false; break; }
        // (a pageable-source cudaMemcpyAsync returns only after the source has
        //  been staged, and p.views is not touched again before this slot is
        //  retired, so no synchronization is needed here)
        if (!ring.BeginUpload(slot, p.inBytes)) { ok = false; break; }
        if (!LaunchFused2DViews(dv, p.count, p.maxCw, p.maxCh, options.block_shape,
                                ring.StreamKernel())) { ok = false; break; }
        if (!ring.EndComputeAndDownload(slot, p.outBytes)) { ok = false; break; }
        st.launches++;
        st.bytesH2D += (int64_t)p.inBytes;
        st.bytesD2H += (int64_t)p.outBytes;

        size_t freeB = 0, totB = 0;
        if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) {
            if (st.deviceFreeMin < 0 || (int64_t)freeB < st.deviceFreeMin) st.deviceFreeMin = (int64_t)freeB;
            if ((int64_t)freeB > st.deviceFreeMax) st.deviceFreeMax = (int64_t)freeB;
        }

        next += p.count;
        slot = (slot + 1) % nslots;
    }
    // drain
    for (int i = 0; i < nslots && ok; i++) {
        if (!retire((slot + i) % nslots)) ok = false;
    }
    cudaFree(d_views);
    ring.Destroy();

    st.wall_ms = (float)std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - t_wall0).count();
    if (stats) *stats = st;
    return ok;
}

} // namespace gpu
} // namespace GInt

namespace GInt {
namespace gpu {

bool ComputeDiscreteGradient2DTiled(const TileSource2D& source,
                                    const TileGrid2D& grid,
                                    const Dgrad2DTables& tables,
                                    const ByteSink2D& sink,
                                    const TiledOptions2D& options,
                                    TileStats2D* stats) {
    return TiledImpl(source, nullptr, grid, tables, sink, options, stats);
}

bool ComputeDiscreteGradient2DConstrainedTiled(const TileSource2D& source,
                                               const ByteSource2D& classSource,
                                               const TileGrid2D& grid,
                                               const Dgrad2DTables& tables,
                                               const ByteSink2D& sink,
                                               const TiledOptions2D& options,
                                               TileStats2D* stats) {
    return TiledImpl(source, &classSource, grid, tables, sink, options, stats);
}

} // namespace gpu
} // namespace GInt
