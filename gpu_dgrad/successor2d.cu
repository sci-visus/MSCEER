// Stage A1: tiled SuccessorResolver + its first client, the tiled flow-
// terminal labeling.
//
// Per tile (lattice box w x h, local index i = x + y*w), a 32-bit word per
// lattice point: tag 00 = local successor index, 10 = TERMINAL(v), 11 =
// FRONTIER(r) with r the ring position of the cell whose successor leaves the
// box (exit target recorded in exitGlobal[r]). In-tile pointer doubling
// W[i] = W[W[i]] runs a FIXED ceil(log2(n))+1 passes - a chain inside a box
// of n cells is at most n long, so that is a proof of convergence, and it
// removes the per-pass host flag of the in-core path entirely. Terminal
// cells get a scan-order rank (CUB exclusive scan) for the compact-id mode.
//
// Two streaming phases over the gradient tiles (see frontier2d.h): phase 1
// downloads only ring words + exit targets + terminal count per tile, the
// host resolves the frontier graph, phase 2 recomputes the doubling (cheap)
// and substitutes the resolved payloads while writing the sink. The device
// footprint is slots x tilesPerLaunch tiles; the host keeps ring records
// only (perimeter-sized).

#include "dgrad2d_internal.cuh"
#include "frontier2d.h"
#include "gpu_transfer.h"
#include "tile2d.h"

#include <cub/cub.cuh>

#include <chrono>
#include <cstdio>
#include <vector>

namespace GInt {
namespace gpu {

namespace {

// device copy of the ring enumeration (tile2d.h)
__host__ __device__ __forceinline__ int ring_index(int x, int y, int w, int h) {
    if (y == 0) return x;
    if (y == h - 1) return w + x;
    if (x == 0) return 2 * w + (y - 1);
    if (x == w - 1) return 2 * w + (h - 2) + (y - 1);
    return -1;
}
__host__ __device__ __forceinline__ void ring_to_local(int r, int w, int h, int& x, int& y) {
    if (r < w) { x = r; y = 0; return; }
    if (h > 1 && r < 2 * w) { x = r - w; y = h - 1; return; }
    const int k = r - 2 * w;
    if (w == 1 || k < h - 2) { x = 0; y = 1 + k; return; }
    x = w - 1;
    y = 1 + k - (h - 2);
}

struct Step2D {
    int dx[8], dy[8];   // per 3-bit pair code; code 7 = terminal
};

// One tile's work in a batch. Device pointers.
struct TileWork {
    const unsigned char* grad;  // owned cell box, origin (2ox, 2oy)
    int gpitch;
    int w, h, n;                // lattice box
    long long ox, oy;           // global lattice origin
    int NX;                     // global lattice width
    int quad;                   // 0: vertex lattice, 1: quad lattice
    int R;                      // ring cells
    int32_t* W;                 // [n] words
    int32_t* flag;              // [n] is-terminal
    int32_t* rank;              // [n] exclusive scan of flag
    int32_t* count;             // [1] terminals (packed output)
    int32_t* ringWord;          // [R] (packed output)
    long long* exitGlobal;      // [R] (packed output)
    // phase 2 only: resolved payload per ring slot, compact-id base, output
    const long long* resolvedSlot;  // [R]
    long long base;
    int compact;
    int32_t* out;               // [n] final int32 words (packed output)
};

__constant__ Step2D c_step;

__global__ void init_words_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int x = i % t.w, y = i / t.w;
    const int cx = t.quad ? 2 * x + 1 : 2 * x;
    const int cy = t.quad ? 2 * y + 1 : 2 * y;
    const unsigned code = (t.grad[cy * t.gpitch + cx] >> 2) & 7u;
    const int r = ring_index(x, y, t.w, t.h);
    if (code == 7u) {
        t.W[i] = (int32_t)((kWordTagTerminal << kWordTagShift) | (uint32_t)i);
        t.flag[i] = 1;
        if (r >= 0) t.exitGlobal[r] = -1;
        return;
    }
    t.flag[i] = 0;
    const int nx = x + c_step.dx[code], ny = y + c_step.dy[code];
    if (nx >= 0 && nx < t.w && ny >= 0 && ny < t.h) {
        t.W[i] = nx + ny * t.w;
        if (r >= 0) t.exitGlobal[r] = -1;
    } else {
        // r >= 0 is guaranteed: a cell stepping out of the box is on the ring
        t.W[i] = (int32_t)((kWordTagFrontier << kWordTagShift) | (uint32_t)r);
        t.exitGlobal[r] = (t.oy + ny) * (long long)t.NX + (t.ox + nx);
    }
}

__global__ void jump_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int32_t w = t.W[i];
    if ((((uint32_t)w) >> kWordTagShift) == kWordTagSucc) t.W[i] = t.W[w];
}

// compact mode: TERMINAL(local) -> TERMINAL(rank)
__global__ void rank_words_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int32_t w = t.W[i];
    if ((((uint32_t)w) >> kWordTagShift) == kWordTagTerminal) {
        const uint32_t v = ((uint32_t)w) & kWordValueMask;
        t.W[i] = (int32_t)((kWordTagTerminal << kWordTagShift) | (uint32_t)t.rank[v]);
    }
}

__global__ void extract_ring_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r == 0) t.count[0] = t.n > 0 ? t.rank[t.n - 1] + t.flag[t.n - 1] : 0;
    if (r >= t.R) return;
    int x, y;
    ring_to_local(r, t.w, t.h, x, y);
    t.ringWord[r] = t.W[x + y * t.w];
}

// phase 2: substitute the final payload per lattice point (device-side; the
// host only copies the box to the sink). -1 marks a payload that does not
// fit an int32, which the host treats as failure.
__global__ void output_words_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int32_t w = t.W[i];
    const uint32_t tag = ((uint32_t)w) >> kWordTagShift;
    const uint32_t v = ((uint32_t)w) & kWordValueMask;
    long long p = -1;
    if (tag == kWordTagTerminal) {
        if (t.compact) p = t.base + (long long)v;
        else p = (t.oy + (long long)(v / (uint32_t)t.w)) * (long long)t.NX + (t.ox + (long long)(v % (uint32_t)t.w));
    } else if (tag == kWordTagFrontier) {
        p = t.resolvedSlot[v];
    }
    t.out[i] = (p < 0 || p > 0x7FFFFFFFll) ? (int32_t)-1 : (int32_t)p;
}

#define SUCC_CUDA_CHECK(call)                                                   \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            return false;                                                       \
        }                                                                       \
    } while (0)

static int ceil_log2_plus1(int n) {
    int p = 1;
    while ((1 << (p - 1)) < n) p++;
    return p + 1;
}

// packed phase-1 output layout per tile (Rmax ring slots reserved):
// [exitGlobal int64 x Rmax][ringWord int32 x Rmax][count int32], stride
// rounded up to 8 bytes so the int64 array stays aligned tile after tile.
struct PackedRing {
    static size_t Bytes(int Rmax) {
        const size_t b = (size_t)Rmax * sizeof(long long) + (size_t)Rmax * sizeof(int32_t) + sizeof(int32_t);
        return (b + 7) & ~(size_t)7;
    }
    static long long* Exits(unsigned char* base, int) { return (long long*)base; }
    static int32_t* Words(unsigned char* base, int Rmax) {
        return (int32_t*)(base + (size_t)Rmax * sizeof(long long));
    }
    static int32_t* Count(unsigned char* base, int Rmax) {
        return (int32_t*)(base + (size_t)Rmax * sizeof(long long) + (size_t)Rmax * sizeof(int32_t));
    }
};

struct Phase {
    bool ringOnly;   // phase 1: download packed ring data; phase 2: full words
};

struct Engine {
    const TileGrid2D& grid;
    int64_t NX, NY;
    bool quad;
    int nslots, tpl, Rmax;
    size_t maxN;
    TransferRing ring;
    TileWork* d_works = nullptr;
    int32_t* d_W = nullptr;
    int32_t* d_flag = nullptr;
    int32_t* d_rank = nullptr;
    void* d_temp = nullptr;
    size_t tempBytes = 0;
    long long* d_resolved = nullptr;   // [nslots * tpl * Rmax]
    int passes = 0;
    int64_t launches = 0;
    const FrontierGraph2D* fgRun = nullptr;   // phase 2 only
    bool compactRun = false;

    Engine(const TileGrid2D& g, bool q, int slots, int tilesPerLaunch)
        : grid(g), NX(q ? g.X - 1 : g.X), NY(q ? g.Y - 1 : g.Y), quad(q),
          nslots(slots), tpl(tilesPerLaunch) {
        Rmax = RingCount2D(g.T, g.T);
        maxN = (size_t)g.T * (size_t)g.T;
    }
    ~Engine() {
        if (d_works) cudaFree(d_works);
        if (d_W) cudaFree(d_W);
        if (d_flag) cudaFree(d_flag);
        if (d_rank) cudaFree(d_rank);
        if (d_temp) cudaFree(d_temp);
        if (d_resolved) cudaFree(d_resolved);
    }

    bool Create() {
        const size_t capIn = (size_t)(2 * grid.T) * (size_t)(2 * grid.T) * (size_t)tpl;
        const size_t capOut1 = PackedRing::Bytes(Rmax) * (size_t)tpl;
        const size_t capOut2 = maxN * sizeof(int32_t) * (size_t)tpl;
        const size_t capOut = capOut1 > capOut2 ? capOut1 : capOut2;
        if (!ring.Create(nslots, capIn, capOut)) return false;
        const size_t perSlot = maxN * (size_t)tpl;
        SUCC_CUDA_CHECK(cudaMalloc(&d_works, (size_t)nslots * tpl * sizeof(TileWork)));
        SUCC_CUDA_CHECK(cudaMalloc(&d_W, (size_t)nslots * perSlot * sizeof(int32_t)));
        SUCC_CUDA_CHECK(cudaMalloc(&d_flag, (size_t)nslots * perSlot * sizeof(int32_t)));
        SUCC_CUDA_CHECK(cudaMalloc(&d_rank, (size_t)nslots * perSlot * sizeof(int32_t)));
        SUCC_CUDA_CHECK(cub::DeviceScan::ExclusiveSum(nullptr, tempBytes, d_flag, d_rank,
                                                      (int)maxN));
        tempBytes = tempBytes * (size_t)nslots + 256;
        SUCC_CUDA_CHECK(cudaMalloc(&d_temp, tempBytes));
        SUCC_CUDA_CHECK(cudaMalloc(&d_resolved, (size_t)nslots * tpl * (size_t)Rmax * sizeof(long long)));
        passes = ceil_log2_plus1((int)maxN);
        return true;
    }
    int64_t DeviceBytes() const {
        return (int64_t)ring.DeviceBytes() + (int64_t)((size_t)nslots * tpl * sizeof(TileWork)) +
               (int64_t)(3 * (size_t)nslots * maxN * (size_t)tpl * sizeof(int32_t)) +
               (int64_t)tempBytes + (int64_t)((size_t)nslots * tpl * (size_t)Rmax * sizeof(long long));
    }

    // Streams every tile through the device; onRetire(slot, firstTile, count,
    // works) consumes the slot's pinned output.
    template <class Retire>
    bool Run(const ByteSource2D& source, const Phase& phase, bool compact,
             const Retire& onRetire) {
        std::vector<std::vector<TileWork>> works((size_t)nslots);
        std::vector<int64_t> firstTile((size_t)nslots, -1);
        std::vector<int> counts((size_t)nslots, 0);
        std::vector<size_t> outBytes((size_t)nslots, 0);

        auto retire = [&](int slot) -> bool {
            if (counts[(size_t)slot] == 0) return true;
            if (!ring.Retire(slot, nullptr, nullptr, nullptr)) return false;
            if (!onRetire(slot, firstTile[(size_t)slot], counts[(size_t)slot],
                          works[(size_t)slot], (unsigned char*)ring.Slot(slot).h_out))
                return false;
            counts[(size_t)slot] = 0;
            return true;
        };

        int64_t next = 0;
        int slot = 0;
        bool ok = true;
        while (next < grid.numTiles && ok) {
            if (!retire(slot)) { ok = false; break; }
            const int cnt = (int)((grid.numTiles - next) < tpl ? (grid.numTiles - next) : tpl);
            std::vector<TileWork>& wk = works[(size_t)slot];
            wk.assign((size_t)cnt, TileWork{});
            TransferSlot& s = ring.Slot(slot);
            unsigned char* h_in = (unsigned char*)s.h_in;
            size_t inOff = 0, outOff2 = 0;
            int maxN_b = 0, maxR_b = 0;
            const size_t perSlot = maxN * (size_t)tpl;
            for (int i = 0; i < cnt && ok; i++) {
                const TileDesc2D d = TileAt(grid, next + i);
                const Box2D ob = OwnedCellBox2D(grid, d);
                const Box2D lb = LatticeBox2D(d, NX, NY);
                if (!source(ob.x0, ob.y0, ob.w, ob.h, h_in + inOff, ob.w)) { ok = false; break; }
                TileWork& t = wk[(size_t)i];
                t.grad = (const unsigned char*)s.d_in + inOff;
                t.gpitch = ob.w;
                t.w = lb.w; t.h = lb.h; t.n = lb.w * lb.h;
                t.ox = lb.x0; t.oy = lb.y0;
                t.NX = (int)NX;
                t.quad = quad ? 1 : 0;
                t.R = RingCount2D(lb.w, lb.h);
                const size_t base = (size_t)slot * perSlot + (size_t)i * maxN;
                t.W = d_W + base;
                t.flag = d_flag + base;
                t.rank = d_rank + base;
                unsigned char* packed = (unsigned char*)s.d_out + (size_t)i * PackedRing::Bytes(Rmax);
                t.count = PackedRing::Count(packed, Rmax);
                t.ringWord = PackedRing::Words(packed, Rmax);
                t.exitGlobal = PackedRing::Exits(packed, Rmax);
                t.resolvedSlot = d_resolved + ((size_t)slot * tpl + (size_t)i) * (size_t)Rmax;
                t.base = 0;
                t.compact = compact ? 1 : 0;
                t.out = nullptr;
                if (!phase.ringOnly) {
                    const TileRingRecord2D& tr = fgRun->tiles[(size_t)(next + i)];
                    t.base = tr.base;
                    t.out = (int32_t*)((unsigned char*)s.d_out + outOff2);
                    outOff2 += (size_t)t.n * sizeof(int32_t);
                    if (tr.R > 0)
                        SUCC_CUDA_CHECK(cudaMemcpyAsync((void*)t.resolvedSlot, tr.resolved.data(),
                                                        (size_t)tr.R * sizeof(long long),
                                                        cudaMemcpyHostToDevice, ring.StreamH2D()));
                }
                if (t.n > maxN_b) maxN_b = t.n;
                if (t.R > maxR_b) maxR_b = t.R;
                inOff += (size_t)ob.w * ob.h;
            }
            if (!ok) break;
            firstTile[(size_t)slot] = next;
            counts[(size_t)slot] = cnt;

            TileWork* dw = d_works + (size_t)slot * tpl;
            SUCC_CUDA_CHECK(cudaMemcpyAsync(dw, wk.data(), (size_t)cnt * sizeof(TileWork),
                                            cudaMemcpyHostToDevice, ring.StreamH2D()));
            if (!ring.BeginUpload(slot, inOff)) { ok = false; break; }
            cudaStream_t ks = ring.StreamKernel();
            if (maxN_b > 0) {
                const int TPB = 256;
                dim3 grid1((unsigned)((maxN_b + TPB - 1) / TPB), (unsigned)cnt);
                init_words_kernel<<<grid1, TPB, 0, ks>>>(dw);
                SUCC_CUDA_CHECK(cudaGetLastError());
                const int np = ceil_log2_plus1(maxN_b);
                for (int p = 0; p < np; p++) jump_kernel<<<grid1, TPB, 0, ks>>>(dw);
                SUCC_CUDA_CHECK(cudaGetLastError());
                // terminal ranks (per tile; CUB temp partitioned per slot)
                size_t tb = tempBytes / (size_t)nslots;
                void* temp = (unsigned char*)d_temp + (size_t)slot * tb;
                for (int i = 0; i < cnt; i++) {
                    const TileWork& t = wk[(size_t)i];
                    if (t.n <= 0) continue;
                    SUCC_CUDA_CHECK(cub::DeviceScan::ExclusiveSum(temp, tb, t.flag, t.rank,
                                                                  t.n, ks));
                }
                if (compact) {
                    rank_words_kernel<<<grid1, TPB, 0, ks>>>(dw);
                    SUCC_CUDA_CHECK(cudaGetLastError());
                }
                dim3 gridR((unsigned)((maxR_b + TPB - 1) / TPB > 0 ? (maxR_b + TPB - 1) / TPB : 1),
                           (unsigned)cnt);
                extract_ring_kernel<<<gridR, TPB, 0, ks>>>(dw);
                SUCC_CUDA_CHECK(cudaGetLastError());
            }
            size_t ob = 0;
            if (phase.ringOnly) {
                ob = (size_t)cnt * PackedRing::Bytes(Rmax);
            } else {
                // phase 2: final words per tile, back to back in d_out
                if (maxN_b > 0) {
                    const int TPB = 256;
                    dim3 grid1((unsigned)((maxN_b + TPB - 1) / TPB), (unsigned)cnt);
                    output_words_kernel<<<grid1, TPB, 0, ks>>>(dw);
                    SUCC_CUDA_CHECK(cudaGetLastError());
                }
                ob = outOff2;
            }
            outBytes[(size_t)slot] = ob;
            if (!ring.EndComputeAndDownload(slot, ob)) { ok = false; break; }
            launches++;
            next += cnt;
            slot = (slot + 1) % nslots;
        }
        for (int i = 0; i < nslots && ok; i++)
            if (!retire((slot + i) % nslots)) ok = false;
        return ok;
    }
};

} // namespace

bool Label2DFlowTerminalsTiled(const ByteSource2D& gradSource,
                               const TileGrid2D& grid,
                               const Dgrad2DTables& tables, bool ascending,
                               const WordSink2D& sink,
                               const ResolveOptions2D& options,
                               ResolveStats2D* stats) {
    if (grid.X < 2 || grid.Y < 2 || grid.T < 1 || grid.numTiles < 1) return false;
    const auto t0 = std::chrono::steady_clock::now();
    ResolveStats2D st;
    st.tiles = grid.numTiles;

    Step2D step{};
    step.dx[tables.dir_px] = 1;  step.dy[tables.dir_px] = 0;
    step.dx[tables.dir_mx] = -1; step.dy[tables.dir_mx] = 0;
    step.dx[tables.dir_py] = 0;  step.dy[tables.dir_py] = 1;
    step.dx[tables.dir_my] = 0;  step.dy[tables.dir_my] = -1;
    SUCC_CUDA_CHECK(cudaMemcpyToSymbol(c_step, &step, sizeof(step)));

    Engine eng(grid, !ascending, options.slots < 1 ? 1 : options.slots,
               options.tilesPerLaunch < 1 ? 1 : options.tilesPerLaunch);
    if (!eng.Create()) return false;
    st.passes = eng.passes;
    st.deviceBytesAllocated = eng.DeviceBytes();
    st.hostPinnedBytes = (int64_t)eng.ring.PinnedBytes();

    // ---- phase 1: ring records ----
    FrontierGraph2D fg;
    fg.NX = eng.NX;
    fg.NY = eng.NY;
    fg.grid = grid;
    fg.tiles.resize((size_t)grid.numTiles);
    const int Rmax = eng.Rmax;
    Phase p1{true};
    bool ok = eng.Run(gradSource, p1, options.compactIds,
        [&](int, int64_t first, int cnt, const std::vector<TileWork>& wk, unsigned char* h_out) {
            for (int i = 0; i < cnt; i++) {
                const TileWork& t = wk[(size_t)i];
                TileRingRecord2D& tr = fg.tiles[(size_t)(first + i)];
                unsigned char* packed = h_out + (size_t)i * PackedRing::Bytes(Rmax);
                tr.w = t.w; tr.h = t.h; tr.ox = t.ox; tr.oy = t.oy; tr.R = t.R;
                tr.count = t.n > 0 ? *PackedRing::Count(packed, Rmax) : 0;
                tr.word.assign(PackedRing::Words(packed, Rmax), PackedRing::Words(packed, Rmax) + t.R);
                tr.exitGlobal.assign(PackedRing::Exits(packed, Rmax),
                                     PackedRing::Exits(packed, Rmax) + t.R);
            }
            return true;
        });
    st.phase1_ms = (float)std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - t0).count();
    if (!ok) return false;

    // ---- stitch ----
    const auto t1 = std::chrono::steady_clock::now();
    int64_t base = 0;
    for (size_t t = 0; t < fg.tiles.size(); t++) {
        fg.tiles[t].base = base;
        base += fg.tiles[t].count;
        st.hostStitchBytes += (int64_t)fg.tiles[t].R * (int64_t)(sizeof(int32_t) + 2 * sizeof(int64_t));
    }
    st.terminals = base;
    if (!ResolveFrontier2D(fg, options.compactIds, &st.stitchSteps)) {
        std::fprintf(stderr, "Label2DFlowTerminalsTiled: malformed frontier graph\n");
        return false;
    }
    for (size_t t = 0; t < fg.tiles.size(); t++) {
        int64_t n = 0;
        for (int r = 0; r < fg.tiles[t].R; r++) if (fg.tiles[t].exitGlobal[(size_t)r] >= 0) n++;
        st.frontierSlots += n;
        if (n > st.frontierSlotsMaxTile) st.frontierSlotsMaxTile = n;
    }
    st.stitch_ms = (float)std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - t1).count();

    // ---- phase 2: output ----
    const auto t2 = std::chrono::steady_clock::now();
    eng.fgRun = &fg;
    eng.compactRun = options.compactIds;
    Phase p2{false};
    ok = eng.Run(gradSource, p2, options.compactIds,
        [&](int, int64_t, int cnt, const std::vector<TileWork>& wk, unsigned char* h_out) {
            size_t off = 0;
            for (int i = 0; i < cnt; i++) {
                const TileWork& t = wk[(size_t)i];
                if (t.n <= 0) continue;
                const int32_t* W = (const int32_t*)(h_out + off);
                off += (size_t)t.n * sizeof(int32_t);
                if (!sink(t.ox, t.oy, t.w, t.h, W, t.w)) return false;
            }
            return true;
        });
    st.phase2_ms = (float)std::chrono::duration<double, std::milli>(
                       std::chrono::steady_clock::now() - t2).count();
    st.launches = eng.launches;
    st.wall_ms = (float)std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - t0).count();
    if (stats) *stats = st;
    return ok;
}

} // namespace gpu
} // namespace GInt
