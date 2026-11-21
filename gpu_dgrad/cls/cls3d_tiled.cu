// Tiled driver for the 3D lower-star kernel (cls3d_kernel.cu), and the
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

#include "cls3d_internal.cuh"
#include "common/gpu_transfer.h"
#include "tile3d.h"

#include <algorithm>
#include <atomic>
#include <functional>
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
    ViewGroups3D groups;   // views split for the interior / boundary kernels
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

// packed byte -> the two gradient bytes it holds (little endian: low nibble first)
struct UnpackLut {
    uint16_t v[256];
    UnpackLut() {
        for (int b = 0; b < 256; b++)
            v[b] = (uint16_t)(UnpackNibble((unsigned)b & 15u) | (UnpackNibble((unsigned)b >> 4) << 8));
    }
};
const UnpackLut kUnpackLut;

// Expands n cells starting at nibble index nib0 of a packed buffer.
inline void UnpackRun(const unsigned char* packed, int64_t nib0, int64_t n, unsigned char* dst) {
    int64_t i = nib0, k = 0;
    if ((i & 1) && k < n) { dst[k++] = UnpackNibble(packed[i >> 1] >> 4); i++; }
    const unsigned char* pp = packed + (i >> 1);
    const int64_t pairs = (n - k) / 2;
    int64_t j = 0;
    for (; j + 4 <= pairs; j += 4) {   // 4 packed bytes -> 8 gradient bytes
        const uint64_t w = (uint64_t)kUnpackLut.v[pp[j]] | ((uint64_t)kUnpackLut.v[pp[j + 1]] << 16) |
                           ((uint64_t)kUnpackLut.v[pp[j + 2]] << 32) | ((uint64_t)kUnpackLut.v[pp[j + 3]] << 48);
        std::memcpy(dst + k + 2 * j, &w, 8);
    }
    for (; j < pairs; j++) {
        const uint16_t w = kUnpackLut.v[pp[j]];
        std::memcpy(dst + k + 2 * j, &w, 2);
    }
    k += 2 * pairs;
    i += 2 * pairs;
    if (k < n) dst[k] = UnpackNibble(packed[i >> 1] & 15u);
}

template <class F> void ParallelSlices(int n, int threads, F fn);
int CopyThreads(int64_t bytes, int requested);

// A sink of packed owned boxes: cell (x, y, z) of the box is nibble
// nib0 + (z * h + y) * w + x of `packed`.
typedef std::function<bool(const Box3D& box, const unsigned char* packed, int64_t nib0)> PackedSink3D;

} // namespace

static bool Tiled3DImpl(const TileSource3D& source, const ByteSource3D* classSource,
                        int keyClamp, const TileGrid3D& grid, const Dgrad3DTables& tables,
                        const ByteSink3D& sink, const TiledOptions3D& options,
                        TileStats3D* stats, const PackedSink3D* packedSink = nullptr) {
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
    const bool pack = options.packOutput;
    TransferRing ring;
    if (!ring.Create(nslots, capIn, capOut)) return false;
    DeviceTileView3D* d_views = nullptr;
    // up to 7 views per tile after the interior / face split
    TILED3_CUDA_CHECK(cudaMalloc(&d_views, (size_t)nslots * tpl * 7 * sizeof(DeviceTileView3D)));
    // packed copies of the slots' outputs (capOut is a multiple of 16)
    std::vector<unsigned char*> d_pack((size_t)nslots, nullptr);
    if (pack)
        for (int i = 0; i < nslots; i++)
            if (cudaMalloc(&d_pack[(size_t)i], capOut / 2) != cudaSuccess) {
                for (auto q : d_pack) if (q) cudaFree(q);
                cudaFree(d_views);
                return false;
            }
    std::vector<unsigned char> unpacked;   // host expansion for a plain sink
    if (pack && !packedSink) unpacked.resize((size_t)ownedMax);
    st.hostPinnedBytes = (int64_t)ring.PinnedBytes();
    st.deviceBytesAllocated =
        (int64_t)ring.DeviceBytes() + (int64_t)((size_t)nslots * tpl * 7 * sizeof(DeviceTileView3D)) +
        (pack ? (int64_t)(capOut / 2) * nslots : 0);

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
            if (!pack) {
                if (!sink(o, out + off, o.w, (int64_t)o.w * o.h)) return false;
            } else if (packedSink) {
                if (!(*packedSink)(o, out, (int64_t)off)) return false;
            } else {
                const int64_t n = o.Count();
                const int nt = CopyThreads(n, options.hostThreads);
                const int64_t chunk = ((n + nt - 1) / nt + 1) & ~(int64_t)1;
                ParallelSlices(nt, nt, [&](int c) {
                    const int64_t a = (int64_t)c * chunk;
                    const int64_t e = std::min<int64_t>(n, a + chunk);
                    if (a < e) UnpackRun(out, (int64_t)off + a, e - a, unpacked.data() + a);
                });
                if (!sink(o, unpacked.data(), o.w, (int64_t)o.w * o.h)) return false;
            }
            off += Align16((size_t)o.Count());
        }
        st.sink_ms += (float)MsSince(t0);
        p.count = 0;
        return true;
    };

    int pendingSlot = -1;
    size_t pendingBytes = 0;
    bool ok = true;
    int64_t next = 0;
    int slot = 0;
    while (next < grid.numTiles && ok) {
        if (pendingSlot == slot) {   // single slot: its download must be queued before reuse
            if (!ring.Download(pendingSlot, pendingBytes, pack ? d_pack[(size_t)pendingSlot] : nullptr)) {
                ok = false;
                break;
            }
            pendingSlot = -1;
        }
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
        DeviceTileView3D* dv = d_views + (size_t)slot * tpl * 7;
        p.groups = SplitViews3D(p.views.data(), p.count);   // interior / domain-face boxes
        if (cudaMemcpyAsync(dv, p.groups.views.data(), p.groups.views.size() * sizeof(DeviceTileView3D),
                            cudaMemcpyHostToDevice, ring.StreamH2D()) != cudaSuccess) {
            ok = false;
            break;
        }
        if (!ring.BeginUpload(slot, p.inBytes)) { ok = false; break; }
        if (!LaunchDgrad3DGroups(dv, p.groups, options.block_shape, ring.StreamKernel())) {
            ok = false;
            break;
        }
        // download deferred until the next slot's upload is queued (copies
        // run in submission order on a single copy engine)
        if (pack && !LaunchPackNibbles((const unsigned char*)s.d_out, d_pack[(size_t)slot], p.outBytes / 2,
                                       ring.StreamKernel())) { ok = false; break; }
        if (!ring.MarkComputed(slot)) { ok = false; break; }
        if (pendingSlot >= 0 && !ring.Download(pendingSlot, pendingBytes,
                                               pack ? d_pack[(size_t)pendingSlot] : nullptr)) { ok = false; break; }
        pendingSlot = slot;
        pendingBytes = pack ? p.outBytes / 2 : p.outBytes;
        st.launches++;
        st.bytesH2D += (int64_t)p.inBytes;
        st.bytesD2H += (int64_t)pendingBytes;

        size_t freeB = 0, totB = 0;
        if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) {
            if (st.deviceFreeMin < 0 || (int64_t)freeB < st.deviceFreeMin)
                st.deviceFreeMin = (int64_t)freeB;
            if ((int64_t)freeB > st.deviceFreeMax) st.deviceFreeMax = (int64_t)freeB;
        }
        next += p.count;
        slot = (slot + 1) % nslots;
    }
    if (ok && pendingSlot >= 0 &&
        !ring.Download(pendingSlot, pendingBytes, pack ? d_pack[(size_t)pendingSlot] : nullptr))
        ok = false;
    for (int i = 0; i < nslots && ok; i++)
        if (!retire((slot + i) % nslots)) ok = false;
    if (!ok) cudaDeviceSynchronize();   // let in-flight work finish before freeing
    for (auto q : d_pack) if (q) cudaFree(q);
    cudaFree(d_views);
    ring.Destroy();

    st.wall_ms = (float)MsSince(t_wall0);
    if (stats) *stats = st;
    return ok;
}

namespace {

// Runs fn(z) for z in [0, n) on up to 'threads' std::threads (the in-core
// adapters copy whole volumes between pageable and pinned memory).
// Runs fn(z) for z in [0, n) on up to 'threads' OpenMP threads (the host
// copies of the in-core adapters, once per slab). OpenMP's idle threads spin
// briefly between regions; a condition-variable pool and per-call
// std::threads both cost more than the copies' parallelism gained at ~1000
// short copies per 1024^3 run.
template <class F>
void ParallelSlices(int n, int threads, F fn) {
    if (threads > n) threads = n;
#if defined(_OPENMP)
    if (threads > 1) {
#pragma omp parallel for num_threads(threads) schedule(static)
        for (int z = 0; z < n; z++) fn(z);
        return;
    }
#endif
    for (int z = 0; z < n; z++) fn(z);
}

int CopyThreads(int64_t bytes, int requested) {
    if (requested > 0) return requested;
    if (bytes < ((int64_t)256 << 10)) return 1;   // below 256 KB one thread is faster
    const unsigned hc = std::thread::hardware_concurrency();
    return hc == 0 ? 4 : (int)std::min<unsigned>(hc, 16u);
}

// ~64 MB of cells per slab: enough to amortize a launch and two DMAs per
// slab, small enough that filling and draining the pipeline is short.
int AutoSlabDepth(int64_t X, int64_t Y, int64_t Z, int requested, bool /*direct*/) {
    if (requested > 0) return (int)std::min<int64_t>(requested, Z);
    const int64_t plane = (2 * X - 1) * (2 * Y - 1);         // cells per cell-plane
    int64_t tz = ((int64_t)64 << 20) / (2 * plane);
    if (tz < 1) tz = 1;
    return (int)std::min<int64_t>(tz, Z);
}

// Copy mode: z-slabs through the pinned TransferRing, host copies on threads.
bool InCoreCopy3D(const float* values, const uint8_t* classes, int keyClamp, int64_t X, int64_t Y,
                  int64_t Z, const Dgrad3DTables& tables, uint8_t* out_grad, int tz, int slots,
                  bool pack, TileStats3D* stats) {
    const TileGrid3D grid = MakeTileGrid3D(X, Y, Z, (int)X, (int)Y, tz);
    TiledOptions3D opt;
    opt.slots = slots;
    opt.tilesPerLaunch = 1;
    opt.packOutput = pack;
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
    // the packed download is expanded straight into the caller's array
    PackedSink3D psink = [&](const Box3D& b, const unsigned char* packed, int64_t nib0) {
        ParallelSlices(b.d, CopyThreads(b.Count(), opt.hostThreads), [&](int z) {
            for (int y = 0; y < b.h; y++)
                UnpackRun(packed, nib0 + ((int64_t)z * b.h + y) * b.w, b.w,
                          out_grad + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0);
        });
        return true;
    };
    return Tiled3DImpl(src, classes ? &csrc : nullptr, keyClamp, grid, tables, sink, opt, stats,
                       opt.packOutput ? &psink : nullptr);
}

// Direct mode: full-XY z-slabs are contiguous in the caller's (page-locked)
// arrays, so each slab is one DMA in and one DMA out, straight between the
// caller's memory and per-slot device buffers. Slot reuse is ordered on the
// device by events (the upload into a slot waits for the kernel that last
// read it, the kernel waits for the download that last drained it), so the
// host only enqueues and synchronizes once at the end. This GPU has a single
// copy engine: uploads and downloads share it, the kernel overlaps both.
bool InCoreDirect3D(const float* values, const uint8_t* classes, int keyClamp, int64_t X, int64_t Y,
                    int64_t Z, const Dgrad3DTables& tables, uint8_t* out_grad, int tz, int nslots,
                    TileStats3D* stats) {
    const auto t_wall0 = std::chrono::steady_clock::now();
    if (X >= kAxisLimit || Y >= kAxisLimit || Z >= kAxisLimit) return false;
    if (classes && (keyClamp < 0 || keyClamp > 255)) return false;
    const TileGrid3D grid = MakeTileGrid3D(X, Y, Z, (int)X, (int)Y, tz);
    const int64_t W = 2 * X - 1, H = 2 * Y - 1;
    const int64_t stagedMax = X * Y * std::min<int64_t>(grid.Tz + 3, Z);
    const int64_t ownedMax = W * H * std::min<int64_t>(2 * (int64_t)grid.Tz, 2 * Z - 1);
    const int64_t clsMax = classes ? W * H * std::min<int64_t>(2 * (int64_t)grid.Tz + 3, 2 * Z - 1) : 0;
    if (stagedMax > kInt32Limit || ownedMax > kInt32Limit || clsMax > kInt32Limit) return false;
    const int64_t ntiles = grid.numTiles;
    if (nslots < 1) nslots = 1;
    if ((int64_t)nslots > ntiles) nslots = (int)ntiles;
    if (!UploadTables3D(tables)) return false;

    const size_t inBytes = Align16((size_t)stagedMax * sizeof(float)) + Align16((size_t)clsMax);
    const size_t outBytes = Align16((size_t)ownedMax);
    std::vector<unsigned char*> d_in((size_t)nslots, nullptr), d_out((size_t)nslots, nullptr);
    DeviceTileView3D* d_views = nullptr;
    cudaStream_t sh = nullptr, sk = nullptr, sd = nullptr;
    std::vector<cudaEvent_t> ev;   // per tile: h2d start/end, kernel start/end, d2h start/end
    bool ok = true;
    auto chk = [&](cudaError_t e) {
        if (e != cudaSuccess && ok) {
            std::fprintf(stderr, "dgrad3d direct: %s\n", cudaGetErrorString(e));
            ok = false;
        }
    };
    for (int i = 0; i < nslots && ok; i++) {
        chk(cudaMalloc(&d_in[(size_t)i], inBytes));
        chk(cudaMalloc(&d_out[(size_t)i], outBytes));
    }
    chk(cudaMalloc(&d_views, (size_t)ntiles * 7 * sizeof(DeviceTileView3D)));
    chk(cudaStreamCreateWithFlags(&sh, cudaStreamNonBlocking));
    chk(cudaStreamCreateWithFlags(&sk, cudaStreamNonBlocking));
    chk(cudaStreamCreateWithFlags(&sd, cudaStreamNonBlocking));
    ev.assign((size_t)ntiles * 6, nullptr);
    for (auto& e : ev)
        if (ok) chk(cudaEventCreate(&e));

    TileStats3D st;
    st.tiles = ntiles;
    st.deviceBytesAllocated = (int64_t)nslots * (int64_t)(inBytes + outBytes) +
                              ntiles * 7 * (int64_t)sizeof(DeviceTileView3D);
    unsigned char dircode[8];
    DirCodes3D(X, Y, dircode);

    // every view up front (each points into its slot's buffers), one upload
    std::vector<DeviceTileView3D> views((size_t)ntiles);
    std::vector<Box3D> sbox((size_t)ntiles), obox((size_t)ntiles), kbox((size_t)ntiles);
    int maxCd = 1;
    for (int64_t t = 0; t < ntiles && ok; t++) {
        const int slot = (int)(t % nslots);
        const TileDesc3D d = TileAt3D(grid, t);
        const Box3D sb = StagedVertexBox3D(grid, d), cb = ComputedVertexBox3D(grid, d),
                    ob = OwnedCellBox3D(grid, d);
        sbox[(size_t)t] = sb;
        obox[(size_t)t] = ob;
        DeviceTileView3D& v = views[(size_t)t];
        std::memset(&v, 0, sizeof(v));
        v.vals = (const float*)d_in[(size_t)slot];
        v.vpy = sb.w; v.vpz = sb.w * sb.h;
        v.sx0 = (int)sb.x0; v.sy0 = (int)sb.y0; v.sz0 = (int)sb.z0;
        v.svx = sb.w; v.svy = sb.h; v.svz = sb.d;
        v.cx0 = (int)cb.x0; v.cy0 = (int)cb.y0; v.cz0 = (int)cb.z0;
        v.cw = cb.w; v.ch = cb.h; v.cd = cb.d;
        v.ox = (int)d.ox; v.oy = (int)d.oy; v.oz = (int)d.oz;
        v.tw = d.tw; v.th = d.th; v.td = d.td;
        v.X = (int)X; v.Y = (int)Y; v.Z = (int)Z;
        v.grad = d_out[(size_t)slot];
        v.gpy = ob.w; v.gpz = ob.w * ob.h;
        v.gcx0 = (int)ob.x0; v.gcy0 = (int)ob.y0; v.gcz0 = (int)ob.z0;
        v.keyClamp = keyClamp;
        if (classes) {
            const Box3D kb = ClassCellBox3D(grid, d);
            kbox[(size_t)t] = kb;
            v.cls = d_in[(size_t)slot] + Align16((size_t)stagedMax * sizeof(float));
            v.kpy = kb.w; v.kpz = kb.w * kb.h;
            v.kcx0 = (int)kb.x0; v.kcy0 = (int)kb.y0; v.kcz0 = (int)kb.z0;
        }
        std::memcpy(v.dircode, dircode, sizeof(dircode));
        maxCd = std::max(maxCd, cb.d);
    }
    // interior / domain-face split per slab, 7 view slots per slab
    std::vector<ViewGroups3D> groups((size_t)ntiles);
    std::vector<DeviceTileView3D> flat((size_t)ntiles * 7);
    for (int64_t t = 0; t < ntiles; t++) {
        groups[(size_t)t] = SplitViews3D(&views[(size_t)t], 1);
        for (size_t k = 0; k < groups[(size_t)t].views.size(); k++) flat[(size_t)t * 7 + k] = groups[(size_t)t].views[k];
    }
    if (ok) chk(cudaMemcpy(d_views, flat.data(), flat.size() * sizeof(DeviceTileView3D), cudaMemcpyHostToDevice));

    // The copy engine runs copies in submission order (this GPU has one), so
    // a download queued right behind its own upload would hold every later
    // upload until its kernel finished. Uploads are therefore issued
    // `ahead` slabs early and each download after the next upload: the engine
    // sees H2D(t+ahead), D2H(t), ... and moves one slab out while the next
    // comes in, both overlapping the kernel in between.
    const int64_t ahead = nslots > 1 ? nslots - 1 : 1;
    auto issueUpload = [&](int64_t t) {
        const int slot = (int)(t % nslots);
        cudaEvent_t* e = &ev[(size_t)t * 6];
        if (t >= nslots) chk(cudaStreamWaitEvent(sh, ev[(size_t)(t - nslots) * 6 + 3], 0));   // slot read
        const Box3D& sb = sbox[(size_t)t];
        chk(cudaEventRecord(e[0], sh));
        // full XY planes: the staged box starts at plane sb.z0 and is contiguous
        chk(cudaMemcpyAsync(d_in[(size_t)slot], values + sb.z0 * X * Y, (size_t)sb.Count() * sizeof(float),
                            cudaMemcpyHostToDevice, sh));
        if (classes) {
            const Box3D& kb = kbox[(size_t)t];
            chk(cudaMemcpyAsync(d_in[(size_t)slot] + Align16((size_t)stagedMax * sizeof(float)),
                                classes + kb.z0 * W * H, (size_t)kb.Count(), cudaMemcpyHostToDevice, sh));
        }
        chk(cudaEventRecord(e[1], sh));
        st.bytesH2D += sb.Count() * 4 + (classes ? kbox[(size_t)t].Count() : 0);
    };
    auto issueKernel = [&](int64_t t) {
        cudaEvent_t* e = &ev[(size_t)t * 6];
        if (t >= nslots) chk(cudaStreamWaitEvent(sk, ev[(size_t)(t - nslots) * 6 + 5], 0));   // slot drained
        chk(cudaStreamWaitEvent(sk, e[1], 0));
        chk(cudaEventRecord(e[2], sk));
        if (ok && !LaunchDgrad3DGroups(d_views + t * 7, groups[(size_t)t], 0, sk)) ok = false;
        chk(cudaEventRecord(e[3], sk));
        st.launches++;
    };
    auto issueDownload = [&](int64_t t) {
        const int slot = (int)(t % nslots);
        cudaEvent_t* e = &ev[(size_t)t * 6];
        const Box3D& ob = obox[(size_t)t];
        chk(cudaStreamWaitEvent(sd, e[3], 0));
        chk(cudaEventRecord(e[4], sd));
        chk(cudaMemcpyAsync(out_grad + ob.z0 * W * H, d_out[(size_t)slot], (size_t)ob.Count(),
                            cudaMemcpyDeviceToHost, sd));
        chk(cudaEventRecord(e[5], sd));
        st.bytesD2H += ob.Count();
    };
    // prologue: the first `ahead` uploads; then per slab: its kernel, the
    // upload `ahead` slabs later (its slot was last read by kernel t + ahead -
    // nslots < t, already queued), and its own download
    for (int64_t t = 0; t < std::min<int64_t>(ahead, ntiles) && ok; t++) issueUpload(t);
    for (int64_t t = 0; t < ntiles && ok; t++) {
        issueKernel(t);
        if (t + ahead < ntiles) issueUpload(t + ahead);
        issueDownload(t);
    }
    chk(cudaDeviceSynchronize());
    if (ok)
        for (int64_t t = 0; t < ntiles; t++) {
            const cudaEvent_t* e = &ev[(size_t)t * 6];
            float a1 = 0, a2 = 0, a3 = 0;
            cudaEventElapsedTime(&a1, e[0], e[1]);
            cudaEventElapsedTime(&a2, e[2], e[3]);
            cudaEventElapsedTime(&a3, e[4], e[5]);
            st.h2d_ms += a1; st.kernel_ms += a2; st.d2h_ms += a3;
        }
    size_t freeB = 0, totB = 0;
    if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) st.deviceFreeMin = st.deviceFreeMax = (int64_t)freeB;
    for (auto& e : ev) if (e) cudaEventDestroy(e);
    if (sh) cudaStreamDestroy(sh);
    if (sk) cudaStreamDestroy(sk);
    if (sd) cudaStreamDestroy(sd);
    for (auto p : d_in) if (p) cudaFree(p);
    for (auto p : d_out) if (p) cudaFree(p);
    if (d_views) cudaFree(d_views);
    st.wall_ms = (float)MsSince(t_wall0);
    if (stats) *stats = st;
    (void)maxCd;
    return ok;
}

bool PinnedRange(const void* ptr, size_t bytes) {
    if (!ptr || bytes == 0) return false;
    const unsigned char* p = (const unsigned char*)ptr;
    const void* ends[2] = {p, p + bytes - 1};
    for (const void* q : ends) {
        cudaPointerAttributes at;
        if (cudaPointerGetAttributes(&at, q) != cudaSuccess) { cudaGetLastError(); return false; }
        if (at.type != cudaMemoryTypeHost) return false;
    }
    return true;
}

bool InCore3D(const float* values, const uint8_t* classes, int keyClamp, int64_t X, int64_t Y,
              int64_t Z, const Dgrad3DTables& tables, uint8_t* out_grad,
              const InCoreOptions3D& opt, Dgrad3DTimings* timings) {
    if (!values || !out_grad || X < 1 || Y < 1 || Z < 1) return false;
    if (X > INT32_MAX || Y > INT32_MAX || Z > INT32_MAX) return false;
    const auto t0 = std::chrono::steady_clock::now();
    const size_t vb = (size_t)(X * Y * Z) * sizeof(float);
    const size_t gb = (size_t)((2 * X - 1) * (2 * Y - 1) * (2 * Z - 1));
    bool direct = false, registered = false;
    if (opt.hostMode == 0) {
        direct = PinnedRange(values, vb) && PinnedRange(out_grad, gb) && (!classes || PinnedRange(classes, gb));
    } else if (opt.hostMode == 2) {
        direct = true;
        registered = true;
        bool ok = cudaHostRegister((void*)values, vb, cudaHostRegisterDefault) == cudaSuccess;
        ok = ok && cudaHostRegister(out_grad, gb, cudaHostRegisterDefault) == cudaSuccess;
        ok = ok && (!classes || cudaHostRegister((void*)classes, gb, cudaHostRegisterDefault) == cudaSuccess);
        if (!ok) {
            cudaGetLastError();
            cudaHostUnregister((void*)values);
            cudaHostUnregister(out_grad);
            if (classes) cudaHostUnregister((void*)classes);
            cudaGetLastError();
            direct = registered = false;   // fall back to copy mode
        }
    }
    TileStats3D st;
    bool ok;
    const int tz = AutoSlabDepth(X, Y, Z, opt.slabDepth, direct);
    if (direct)
        ok = InCoreDirect3D(values, classes, keyClamp, X, Y, Z, tables, out_grad, tz, opt.slots, &st);
    else
        ok = InCoreCopy3D(values, classes, keyClamp, X, Y, Z, tables, out_grad, tz, opt.slots,
                          opt.packOutput, &st);
    if (registered) {
        cudaHostUnregister((void*)values);
        cudaHostUnregister(out_grad);
        if (classes) cudaHostUnregister((void*)classes);
    }
    if (timings) {
        timings->h2d_ms = st.h2d_ms;
        timings->kernel_ms = st.kernel_ms;
        timings->d2h_ms = st.d2h_ms;
        timings->wall_ms = (float)MsSince(t0);
    }
    return ok;
}

} // namespace

bool HostParallelFor(int n, int threads, const std::function<void(int)>& f) {
    ParallelSlices(n, threads, f);
    return true;
}

int HostCopyThreads(int64_t bytes) { return CopyThreads(bytes, 0); }

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
    return InCore3D(values, nullptr, 0, X, Y, Z, tables, out_grad, InCoreOptions3D(), timings);
}

bool ComputeDiscreteGradient3DConstrained(const float* values, const uint8_t* classes,
                                          int64_t X, int64_t Y, int64_t Z,
                                          const Dgrad3DTables& tables, int keyClamp,
                                          uint8_t* out_grad, Dgrad3DTimings* timings) {
    if (!classes) return false;
    return InCore3D(values, classes, keyClamp, X, Y, Z, tables, out_grad, InCoreOptions3D(), timings);
}

bool ComputeDiscreteGradient3DEx(const float* values, const uint8_t* classes, int keyClamp,
                                 int64_t X, int64_t Y, int64_t Z, const Dgrad3DTables& tables,
                                 uint8_t* out_grad, const InCoreOptions3D& options,
                                 Dgrad3DTimings* timings) {
    return InCore3D(values, classes, keyClamp, X, Y, Z, tables, out_grad, options, timings);
}

bool PinHostMemory(void* ptr, size_t bytes) {
    if (!ptr || bytes == 0) return false;
    if (cudaHostRegister(ptr, bytes, cudaHostRegisterDefault) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    return true;
}

bool UnpinHostMemory(void* ptr) {
    if (!ptr) return false;
    if (cudaHostUnregister(ptr) != cudaSuccess) {
        cudaGetLastError();
        return false;
    }
    return true;
}

bool IsHostMemoryPinned(const void* ptr, size_t bytes) { return PinnedRange(ptr, bytes); }

} // namespace gpu
} // namespace GInt
