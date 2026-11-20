// Tiled driver for the 3D lower-star kernel (dgrad3d_kernel.cu), and the
// in-core entry points, which are the same driver with one tile.
//
// Tiles are visited in tile-id order and packed 'tilesPerLaunch' at a time
// into a TransferRing slot: staged value boxes (and class boxes) back to back
// in the slot's pinned input buffer, owned cell boxes back to back in its
// output buffer, one launch per slot. While the device runs H2D / kernel /
// D2H of one slot on three streams, the host fills the next slot from the
// source and drains a finished one into the sink. Nothing is allocated in
// proportion to the volume: 'slots' x 'tilesPerLaunch' tiles is the whole
// footprint, host and device.
//
// All offsets and boxes are int64 here; the device views carry only int32
// tile-local pitches and per-axis coordinates (checked below).

#include "dgrad3d_internal.cuh"
#include "gpu_transfer.h"
#include "tile3d.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

namespace GInt {
namespace gpu {

namespace {

#define TILED3_CUDA_CHECK(call)                                                 \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            return false;                                                       \
        }                                                                       \
    } while (0)

struct BatchPlan3 {
    int count = 0;
    size_t inBytes = 0, outBytes = 0;
    int maxCw = 0, maxCh = 0, maxCd = 0;
    std::vector<DeviceTileView3D> views;
    std::vector<Box3D> owned;
};

double MsSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
        .count();
}

const int64_t kAxisLimit = (int64_t)1 << 30;
const int64_t kInt32Limit = ((int64_t)1 << 31) - 1;

inline size_t Align16(size_t b) { return (b + 15) & ~(size_t)15; }

int64_t AxisCap(int T, int64_t N, int extra) {
    const int64_t a = (int64_t)T + extra;
    return a < N ? a : N;
}

} // namespace

static bool Tiled3DImpl(const TileSource3D& source, const ByteSource3D* classSource,
                        int keyClamp, const TileGrid3D& grid, const Dgrad3DTables& tables,
                        const ByteSink3D& sink, const TiledOptions3D& options,
                        TileStats3D* stats) {
    if (grid.X < 1 || grid.Y < 1 || grid.Z < 1 || grid.numTiles < 1) return false;
    if (grid.X >= kAxisLimit || grid.Y >= kAxisLimit || grid.Z >= kAxisLimit) {
        std::fprintf(stderr, "dgrad3d: axis too large for 32-bit tile-local indexing\n");
        return false;
    }
    if (classSource && (keyClamp < 0 || keyClamp > 255)) return false;

    // per-tile capacities (the largest tile: interior, T+3 staged per axis)
    const int64_t sx = AxisCap(grid.Tx, grid.X, 3), sy = AxisCap(grid.Ty, grid.Y, 3),
                  sz = AxisCap(grid.Tz, grid.Z, 3);
    const int64_t ocx = std::min<int64_t>(2 * (int64_t)grid.Tx, 2 * grid.X - 1),
                  ocy = std::min<int64_t>(2 * (int64_t)grid.Ty, 2 * grid.Y - 1),
                  ocz = std::min<int64_t>(2 * (int64_t)grid.Tz, 2 * grid.Z - 1);
    // class box: the computed box is up to T+1 wide, its stars span 2(T+1)+1 cells
    const int64_t kcx = std::min<int64_t>(2 * (int64_t)grid.Tx + 3, 2 * grid.X - 1),
                  kcy = std::min<int64_t>(2 * (int64_t)grid.Ty + 3, 2 * grid.Y - 1),
                  kcz = std::min<int64_t>(2 * (int64_t)grid.Tz + 3, 2 * grid.Z - 1);
    const int64_t stagedMax = sx * sy * sz, ownedMax = ocx * ocy * ocz,
                  clsMax = classSource ? kcx * kcy * kcz : 0;
    if (stagedMax > kInt32Limit || ownedMax > kInt32Limit || clsMax > kInt32Limit) {
        std::fprintf(stderr, "dgrad3d: tile too large for 32-bit tile-local indexing\n");
        return false;
    }

    int tpl = options.tilesPerLaunch < 1 ? 1 : options.tilesPerLaunch;
    if ((int64_t)tpl > grid.numTiles) tpl = (int)grid.numTiles;
    int nslots = options.slots < 1 ? 1 : options.slots;
    const int64_t batches = (grid.numTiles + tpl - 1) / tpl;
    if ((int64_t)nslots > batches) nslots = (int)batches;

    const auto t_wall0 = std::chrono::steady_clock::now();
    TileStats3D st;
    st.tiles = grid.numTiles;

    const size_t perTileIn = Align16((size_t)stagedMax * sizeof(float)) + Align16((size_t)clsMax);
    const size_t perTileOut = Align16((size_t)ownedMax);
    const size_t capIn = perTileIn * (size_t)tpl;
    const size_t capOut = perTileOut * (size_t)tpl;

    if (!UploadTables3D(tables)) {
        std::fprintf(stderr, "dgrad3d: tables rejected\n");
        return false;
    }
    TransferRing ring;
    if (!ring.Create(nslots, capIn, capOut)) return false;
    DeviceTileView3D* d_views = nullptr;
    TILED3_CUDA_CHECK(cudaMalloc(&d_views, (size_t)nslots * tpl * sizeof(DeviceTileView3D)));
    st.hostPinnedBytes = (int64_t)ring.PinnedBytes();
    st.deviceBytesAllocated =
        (int64_t)ring.DeviceBytes() + (int64_t)((size_t)nslots * tpl * sizeof(DeviceTileView3D));

    unsigned char dircode[8];
    DirCodes3D(grid.X, grid.Y, dircode);

    std::vector<BatchPlan3> plans((size_t)nslots);

    auto retire = [&](int slot) -> bool {
        BatchPlan3& p = plans[(size_t)slot];
        if (p.count == 0) return true;
        float h2d = 0.f, kern = 0.f, d2h = 0.f;
        if (!ring.Retire(slot, &h2d, &kern, &d2h)) return false;
        st.h2d_ms += h2d;
        st.kernel_ms += kern;
        st.d2h_ms += d2h;
        const auto t0 = std::chrono::steady_clock::now();
        const unsigned char* out = (const unsigned char*)ring.Slot(slot).h_out;
        size_t off = 0;
        for (int i = 0; i < p.count; i++) {
            const Box3D& o = p.owned[(size_t)i];
            if (!sink(o, out + off, o.w, (int64_t)o.w * o.h)) return false;
            off += Align16((size_t)o.Count());
        }
        st.sink_ms += (float)MsSince(t0);
        p.count = 0;
        return true;
    };

    bool ok = true;
    int64_t next = 0;
    int slot = 0;
    while (next < grid.numTiles && ok) {
        if (!retire(slot)) { ok = false; break; }
        BatchPlan3& p = plans[(size_t)slot];
        p.count = (int)std::min<int64_t>(grid.numTiles - next, tpl);
        p.views.resize((size_t)p.count);
        p.owned.resize((size_t)p.count);
        p.maxCw = p.maxCh = p.maxCd = 0;

        TransferSlot& s = ring.Slot(slot);
        unsigned char* h_in = (unsigned char*)s.h_in;
        const unsigned char* d_in = (const unsigned char*)s.d_in;
        size_t inOff = 0, outOff = 0;
        const auto t_src = std::chrono::steady_clock::now();
        for (int i = 0; i < p.count && ok; i++) {
            const TileDesc3D d = TileAt3D(grid, next + i);
            const Box3D sb = StagedVertexBox3D(grid, d);
            const Box3D cb = ComputedVertexBox3D(grid, d);
            const Box3D ob = OwnedCellBox3D(grid, d);
            p.owned[(size_t)i] = ob;
            if (sb.Count() > stagedMax || ob.Count() > ownedMax ||
                inOff + perTileIn > capIn || outOff + perTileOut > capOut) {
                std::fprintf(stderr, "dgrad3d: internal error, tile exceeds its buffer\n");
                ok = false;
                break;
            }
            if (!source(sb, (float*)(h_in + inOff), sb.w, (int64_t)sb.w * sb.h)) { ok = false; break; }

            DeviceTileView3D& v = p.views[(size_t)i];
            std::memset(&v, 0, sizeof(v));
            v.vals = (const float*)(d_in + inOff);
            inOff += Align16((size_t)sb.Count() * sizeof(float));
            v.vpy = sb.w;
            v.vpz = sb.w * sb.h;
            v.sx0 = (int)sb.x0; v.sy0 = (int)sb.y0; v.sz0 = (int)sb.z0;
            v.svx = sb.w; v.svy = sb.h; v.svz = sb.d;
            v.cx0 = (int)cb.x0; v.cy0 = (int)cb.y0; v.cz0 = (int)cb.z0;
            v.cw = cb.w; v.ch = cb.h; v.cd = cb.d;
            v.ox = (int)d.ox; v.oy = (int)d.oy; v.oz = (int)d.oz;
            v.tw = d.tw; v.th = d.th; v.td = d.td;
            v.X = (int)grid.X; v.Y = (int)grid.Y; v.Z = (int)grid.Z;
            v.grad = (unsigned char*)s.d_out + outOff;
            v.gpy = ob.w;
            v.gpz = ob.w * ob.h;
            v.gcx0 = (int)ob.x0; v.gcy0 = (int)ob.y0; v.gcz0 = (int)ob.z0;
            outOff += Align16((size_t)ob.Count());
            v.cls = nullptr;
            v.keyClamp = keyClamp;
            if (classSource) {
                const Box3D kb = ClassCellBox3D(grid, d);
                if (kb.Count() > clsMax) {
                    std::fprintf(stderr, "dgrad3d: internal error, class box exceeds its buffer\n");
                    ok = false;
                    break;
                }
                if (!(*classSource)(kb, h_in + inOff, kb.w, (int64_t)kb.w * kb.h)) { ok = false; break; }
                v.cls = d_in + inOff;
                v.kpy = kb.w;
                v.kpz = kb.w * kb.h;
                v.kcx0 = (int)kb.x0; v.kcy0 = (int)kb.y0; v.kcz0 = (int)kb.z0;
                inOff += Align16((size_t)kb.Count());
            }
            std::memcpy(v.dircode, dircode, sizeof(dircode));
            p.maxCw = std::max(p.maxCw, cb.w);
            p.maxCh = std::max(p.maxCh, cb.h);
            p.maxCd = std::max(p.maxCd, cb.d);
        }
        st.source_ms += (float)MsSince(t_src);
        if (!ok) break;
        p.inBytes = inOff;
        p.outBytes = outOff;

        // the views go onto the H2D stream before the values, so the event the
        // compute stream waits on covers both; this slot's previous batch was
        // retired above, so nothing still reads this region of d_views
        DeviceTileView3D* dv = d_views + (size_t)slot * tpl;
        if (cudaMemcpyAsync(dv, p.views.data(), (size_t)p.count * sizeof(DeviceTileView3D),
                            cudaMemcpyHostToDevice, ring.StreamH2D()) != cudaSuccess) {
            ok = false;
            break;
        }
        if (!ring.BeginUpload(slot, p.inBytes)) { ok = false; break; }
        if (!LaunchDgrad3DViews(dv, p.count, p.maxCw, p.maxCh, p.maxCd, options.block_shape,
                                ring.StreamKernel())) {
            ok = false;
            break;
        }
        if (!ring.EndComputeAndDownload(slot, p.outBytes)) { ok = false; break; }
        st.launches++;
        st.bytesH2D += (int64_t)p.inBytes;
        st.bytesD2H += (int64_t)p.outBytes;

        size_t freeB = 0, totB = 0;
        if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) {
            if (st.deviceFreeMin < 0 || (int64_t)freeB < st.deviceFreeMin)
                st.deviceFreeMin = (int64_t)freeB;
            if ((int64_t)freeB > st.deviceFreeMax) st.deviceFreeMax = (int64_t)freeB;
        }
        next += p.count;
        slot = (slot + 1) % nslots;
    }
    for (int i = 0; i < nslots && ok; i++)
        if (!retire((slot + i) % nslots)) ok = false;
    if (!ok) cudaDeviceSynchronize();   // let in-flight work finish before freeing
    cudaFree(d_views);
    ring.Destroy();

    st.wall_ms = (float)MsSince(t_wall0);
    if (stats) *stats = st;
    return ok;
}

namespace {

// Runs fn(z) for z in [0, n) on up to 'threads' std::threads (the in-core
// adapters copy whole volumes between pageable and pinned memory).
template <class F>
void ParallelSlices(int n, int threads, F fn) {
    if (threads <= 1 || n <= 1) {
        for (int z = 0; z < n; z++) fn(z);
        return;
    }
    if (threads > n) threads = n;
    std::atomic<int> nextz(0);
    std::vector<std::thread> pool;
    pool.reserve((size_t)threads);
    for (int t = 0; t < threads; t++)
        pool.emplace_back([&]() {
            for (int z = nextz++; z < n; z = nextz++) fn(z);
        });
    for (auto& th : pool) th.join();
}

int CopyThreads(int64_t bytes, int requested) {
    if (requested > 0) return requested;
    if (bytes < ((int64_t)8 << 20)) return 1;
    const unsigned hc = std::thread::hardware_concurrency();
    return hc == 0 ? 4 : (int)std::min<unsigned>(hc, 16u);
}

bool InCore3D(const float* values, const uint8_t* classes, int keyClamp, int64_t X, int64_t Y,
              int64_t Z, const Dgrad3DTables& tables, uint8_t* out_grad,
              Dgrad3DTimings* timings) {
    if (!values || !out_grad || X < 1 || Y < 1 || Z < 1) return false;
    if (X > INT32_MAX || Y > INT32_MAX || Z > INT32_MAX) return false;
    const TileGrid3D grid = MakeTileGrid3D(X, Y, Z, (int)X, (int)Y, (int)Z);
    TiledOptions3D opt;
    opt.slots = 1;
    opt.tilesPerLaunch = 1;
    const int64_t W = 2 * X - 1, H = 2 * Y - 1;
    TileSource3D src = [&](const Box3D& b, float* dst, int64_t pitch, int64_t spitch) {
        ParallelSlices(b.d, CopyThreads(b.Count() * 4, opt.hostThreads), [&](int z) {
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch,
                            values + ((b.z0 + z) * Y + (b.y0 + y)) * X + b.x0,
                            sizeof(float) * (size_t)b.w);
        });
        return true;
    };
    ByteSource3D csrc = [&](const Box3D& b, uint8_t* dst, int64_t pitch, int64_t spitch) {
        ParallelSlices(b.d, CopyThreads(b.Count(), opt.hostThreads), [&](int z) {
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch,
                            classes + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0, (size_t)b.w);
        });
        return true;
    };
    ByteSink3D sink = [&](const Box3D& b, const uint8_t* srcb, int64_t pitch, int64_t spitch) {
        ParallelSlices(b.d, CopyThreads(b.Count(), opt.hostThreads), [&](int z) {
            for (int y = 0; y < b.h; y++)
                std::memcpy(out_grad + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0,
                            srcb + z * spitch + y * pitch, (size_t)b.w);
        });
        return true;
    };
    TileStats3D st;
    const bool ok = Tiled3DImpl(src, classes ? &csrc : nullptr, keyClamp, grid, tables, sink,
                                opt, &st);
    if (timings) {
        timings->h2d_ms = st.h2d_ms;
        timings->kernel_ms = st.kernel_ms;
        timings->d2h_ms = st.d2h_ms;
        timings->wall_ms = st.wall_ms;
    }
    return ok;
}

} // namespace

bool ComputeDiscreteGradient3DTiled(const TileSource3D& source, const TileGrid3D& grid,
                                    const Dgrad3DTables& tables, const ByteSink3D& sink,
                                    const TiledOptions3D& options, TileStats3D* stats) {
    return Tiled3DImpl(source, nullptr, 0, grid, tables, sink, options, stats);
}

bool ComputeDiscreteGradient3DConstrainedTiled(const TileSource3D& source,
                                               const ByteSource3D& classSource,
                                               const TileGrid3D& grid,
                                               const Dgrad3DTables& tables, int keyClamp,
                                               const ByteSink3D& sink,
                                               const TiledOptions3D& options,
                                               TileStats3D* stats) {
    return Tiled3DImpl(source, &classSource, keyClamp, grid, tables, sink, options, stats);
}

bool ComputeDiscreteGradient3D(const float* values, int64_t X, int64_t Y, int64_t Z,
                               const Dgrad3DTables& tables, uint8_t* out_grad,
                               Dgrad3DTimings* timings) {
    return InCore3D(values, nullptr, 0, X, Y, Z, tables, out_grad, timings);
}

bool ComputeDiscreteGradient3DConstrained(const float* values, const uint8_t* classes,
                                          int64_t X, int64_t Y, int64_t Z,
                                          const Dgrad3DTables& tables, int keyClamp,
                                          uint8_t* out_grad, Dgrad3DTimings* timings) {
    if (!classes) return false;
    return InCore3D(values, classes, keyClamp, X, Y, Z, tables, out_grad, timings);
}

} // namespace gpu
} // namespace GInt
