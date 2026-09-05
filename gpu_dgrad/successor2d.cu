// Stage A1/A2: tiled SuccessorResolver and its clients - the tiled flow-
// terminal labeling (A1) and the numeric region labeling (A2).
//
// Per tile (lattice box w x h, local index i = x + y*w), a 32-bit word per
// lattice point: tag 00 = local successor index, 10 = TERMINAL(v), 11 =
// FRONTIER(r) with r the position on the tile's OUTER ring (the one-cell
// layer around the box) of the exit target - a successor is never more than
// one step outside the box, so the outer ring enumerates every target and
// the target is implied by (tile, r) (frontier2d.h). In-tile pointer
// doubling W[i] = W[W[i]] runs a FIXED ceil(log2(n))+1 passes - a chain
// inside a box of n cells is at most n long, so that is a proof of
// convergence, and it removes the per-pass host flag of the in-core path
// entirely. Terminal cells get a scan-order rank (CUB exclusive scan) for the
// compact-id payload mode.
//
// Two streaming phases over the inputs: phase 1 downloads only each tile's
// inner-ring root words and terminal count, the host resolves the frontier
// graph, phase 2 recomputes the successor field + doubling (cheap; nothing
// proportional to the image is ever held) and a device gather kernel
// substitutes the resolved payloads straight into the output tile.

#include "dgrad2d_internal.cuh"
#include "frontier2d.h"
#include "gpu_transfer.h"
#include "numeric2d.cuh"
#include "tile2d.h"

#include <cub/cub.cuh>

#include <chrono>
#include <cstdio>
#include <functional>
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

enum Kind { kKindFlow = 0, kKindNumeric = 1 };

struct Step2D {
    int dx[8], dy[8];   // per 3-bit pair code; code 7 = terminal
};

// One tile's work in a batch. Device pointers.
struct TileWork {
    int w, h, n;                // lattice box
    long long ox, oy;           // global lattice origin
    int NX;                     // global lattice width
    int R, Rout;                // inner ring cells, outer ring slots
    int payloadMode;            // PayloadMode
    int32_t* W;                 // [n] words
    int32_t* flag;              // [n] is-terminal
    int32_t* rank;              // [n] exclusive scan of flag
    int32_t* count;             // [1] terminals (packed output)
    int32_t* ringWord;          // [R] (packed output)
    const long long* resolvedSlot;  // [Rout] phase 2
    long long base;             // compact-id base, phase 2
    int32_t* out;               // [n] final words, phase 2
    // flow kind
    const unsigned char* grad;  // owned cell box, origin (2ox, 2oy)
    int gpitch;
    int quad;                   // 0: vertex lattice, 1: quad lattice
    // numeric kind
    numeric2d::TileIn nin;
    int* events;
};

__constant__ Step2D c_step;
__constant__ numeric2d::Params c_np;

// word for a successor at global lattice point (gx, gy)
__device__ __forceinline__ int32_t succ_word(const TileWork& t, long long gx, long long gy) {
    const long long lx = gx - t.ox, ly = gy - t.oy;
    if (lx >= 0 && lx < t.w && ly >= 0 && ly < t.h) return (int32_t)(lx + ly * t.w);
    const int r = ring_index((int)lx + 1, (int)ly + 1, t.w + 2, t.h + 2);
    return (int32_t)((kWordTagFrontier << kWordTagShift) | (uint32_t)r);
}

__global__ void init_flow_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int x = i % t.w, y = i / t.w;
    const int cx = t.quad ? 2 * x + 1 : 2 * x;
    const int cy = t.quad ? 2 * y + 1 : 2 * y;
    const unsigned code = (t.grad[cy * t.gpitch + cx] >> 2) & 7u;
    if (code == 7u) {
        t.W[i] = (int32_t)((kWordTagTerminal << kWordTagShift) | (uint32_t)i);
        t.flag[i] = 1;
        return;
    }
    t.flag[i] = 0;
    t.W[i] = succ_word(t, t.ox + x + c_step.dx[code], t.oy + y + c_step.dy[code]);
}

__global__ void init_numeric_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int x = i % t.w, y = i / t.w;
    int living = -1, ev = 0;
    const long long s = numeric2d::hop(t.nin, c_np, (int)(t.ox + x), (int)(t.oy + y), living, ev);
    atomicAdd(t.events + ev, 1);
    if (s < 0) {
        const uint32_t v = living < 0 ? 0u : (uint32_t)living;
        t.W[i] = (int32_t)((kWordTagTerminal << kWordTagShift) | v);
        t.flag[i] = 1;
        return;
    }
    t.flag[i] = 0;
    t.W[i] = succ_word(t, s % (long long)t.NX, s / (long long)t.NX);
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
// fit an int32.
__global__ void output_words_kernel(const TileWork* __restrict__ works) {
    const TileWork t = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= t.n) return;
    const int32_t w = t.W[i];
    const uint32_t tag = ((uint32_t)w) >> kWordTagShift;
    const uint32_t v = ((uint32_t)w) & kWordValueMask;
    long long p = -1;
    if (tag == kWordTagTerminal) {
        if (t.payloadMode == kPayloadDirect) p = (long long)v;
        else if (t.payloadMode == kPayloadCompactId) p = t.base + (long long)v;
        else p = (t.oy + (long long)(v / (uint32_t)t.w)) * (long long)t.NX +
                 (t.ox + (long long)(v % (uint32_t)t.w));
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

// packed phase-1 output layout per tile: [ringWord int32 x Rmax][count int32]
struct PackedRing {
    static size_t Bytes(int Rmax) { return ((size_t)Rmax + 1) * sizeof(int32_t); }
    static int32_t* Words(unsigned char* base) { return (int32_t*)base; }
    static int32_t* Count(unsigned char* base, int Rmax) { return (int32_t*)base + Rmax; }
};

enum PhaseMode { kPhaseRing = 0, kPhaseOutput = 1, kPhaseRawWords = 2 };

// Fills the slot's pinned input for one tile and the TileWork's input
// pointers; advances inOff by the bytes consumed (multiple of 4).
typedef std::function<bool(int64_t tileId, const TileDesc2D& d, unsigned char* h_in,
                           const unsigned char* d_in, size_t& inOff, TileWork& t)>
    PackFn;

struct Engine {
    const TileGrid2D& grid;
    int64_t NX, NY;
    Kind kind;
    int nslots, tpl, Rmax, RoutMax;
    size_t maxN;
    size_t capInPerTile;
    TransferRing ring;
    TileWork* d_works = nullptr;
    int32_t* d_W = nullptr;
    int32_t* d_flag = nullptr;
    int32_t* d_rank = nullptr;
    void* d_temp = nullptr;
    size_t tempBytes = 0;
    long long* d_resolved = nullptr;   // [nslots * tpl * RoutMax]
    int* d_events = nullptr;           // [16]
    int passes = 0;
    int64_t launches = 0;
    PayloadMode payloadMode = kPayloadLatticeIndex;
    const FrontierGraph2D* fgRun = nullptr;   // phase 2 only

    Engine(const TileGrid2D& g, Kind k, bool quadLattice, int slots, int tilesPerLaunch,
           size_t capInPerTile_)
        : grid(g), NX(quadLattice ? g.X - 1 : g.X), NY(quadLattice ? g.Y - 1 : g.Y), kind(k),
          nslots(slots), tpl(tilesPerLaunch), capInPerTile(capInPerTile_) {
        Rmax = RingCount2D(g.T, g.T);
        RoutMax = RingCount2D(g.T + 2, g.T + 2);
        maxN = (size_t)g.T * (size_t)g.T;
    }
    ~Engine() {
        if (d_works) cudaFree(d_works);
        if (d_W) cudaFree(d_W);
        if (d_flag) cudaFree(d_flag);
        if (d_rank) cudaFree(d_rank);
        if (d_temp) cudaFree(d_temp);
        if (d_resolved) cudaFree(d_resolved);
        if (d_events) cudaFree(d_events);
    }

    bool Create() {
        const size_t capIn = capInPerTile * (size_t)tpl;
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
        SUCC_CUDA_CHECK(cudaMalloc(&d_resolved,
                                   (size_t)nslots * tpl * (size_t)RoutMax * sizeof(long long)));
        SUCC_CUDA_CHECK(cudaMalloc(&d_events, 16 * sizeof(int)));
        SUCC_CUDA_CHECK(cudaMemset(d_events, 0, 16 * sizeof(int)));
        passes = ceil_log2_plus1((int)maxN);
        return true;
    }
    int64_t DeviceBytes() const {
        return (int64_t)ring.DeviceBytes() + (int64_t)((size_t)nslots * tpl * sizeof(TileWork)) +
               (int64_t)(3 * (size_t)nslots * maxN * (size_t)tpl * sizeof(int32_t)) +
               (int64_t)tempBytes +
               (int64_t)((size_t)nslots * tpl * (size_t)RoutMax * sizeof(long long));
    }

    // Streams every tile through the device; onRetire(slot, firstTile, count,
    // works, h_out) consumes the slot's pinned output.
    template <class Retire>
    bool Run(const PackFn& pack, PhaseMode phase, const Retire& onRetire) {
        std::vector<std::vector<TileWork>> works((size_t)nslots);
        std::vector<int64_t> firstTile((size_t)nslots, -1);
        std::vector<int> counts((size_t)nslots, 0);

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
        const size_t perSlot = maxN * (size_t)tpl;
        while (next < grid.numTiles && ok) {
            if (!retire(slot)) { ok = false; break; }
            const int cnt = (int)((grid.numTiles - next) < tpl ? (grid.numTiles - next) : tpl);
            std::vector<TileWork>& wk = works[(size_t)slot];
            wk.assign((size_t)cnt, TileWork{});
            TransferSlot& s = ring.Slot(slot);
            size_t inOff = 0, outOff2 = 0;
            int maxN_b = 0, maxR_b = 0;
            for (int i = 0; i < cnt && ok; i++) {
                const int64_t tid = next + i;
                const TileDesc2D d = TileAt(grid, tid);
                const Box2D lb = LatticeBox2D(d, NX, NY);
                TileWork& t = wk[(size_t)i];
                t.w = lb.w; t.h = lb.h; t.n = lb.w * lb.h;
                t.ox = lb.x0; t.oy = lb.y0;
                t.NX = (int)NX;
                t.R = RingCount2D(lb.w, lb.h);
                t.Rout = RingCount2D(lb.w + 2, lb.h + 2);
                t.payloadMode = (int)payloadMode;
                if (!pack(tid, d, (unsigned char*)s.h_in, (const unsigned char*)s.d_in, inOff, t)) {
                    ok = false;
                    break;
                }
                const size_t base = (size_t)slot * perSlot + (size_t)i * maxN;
                t.W = d_W + base;
                t.flag = d_flag + base;
                t.rank = d_rank + base;
                unsigned char* packed = (unsigned char*)s.d_out + (size_t)i * PackedRing::Bytes(Rmax);
                t.count = PackedRing::Count(packed, Rmax);
                t.ringWord = PackedRing::Words(packed);
                t.resolvedSlot = d_resolved + ((size_t)slot * tpl + (size_t)i) * (size_t)RoutMax;
                t.base = 0;
                t.out = nullptr;
                t.events = d_events;
                if (phase == kPhaseOutput) {
                    const TileRingRecord2D& tr = fgRun->tiles[(size_t)tid];
                    t.base = tr.base;
                    t.out = (int32_t*)((unsigned char*)s.d_out + outOff2);
                    outOff2 += (size_t)t.n * sizeof(int32_t);
                    if (tr.Rout > 0)
                        SUCC_CUDA_CHECK(cudaMemcpyAsync((void*)t.resolvedSlot, tr.resolved.data(),
                                                        (size_t)tr.Rout * sizeof(long long),
                                                        cudaMemcpyHostToDevice, ring.StreamH2D()));
                } else if (phase == kPhaseRawWords) {
                    t.out = (int32_t*)((unsigned char*)s.d_out + outOff2);
                    outOff2 += (size_t)t.n * sizeof(int32_t);
                }
                if (t.n > maxN_b) maxN_b = t.n;
                if (t.R > maxR_b) maxR_b = t.R;
            }
            if (!ok) break;
            firstTile[(size_t)slot] = next;
            counts[(size_t)slot] = cnt;

            TileWork* dw = d_works + (size_t)slot * tpl;
            SUCC_CUDA_CHECK(cudaMemcpyAsync(dw, wk.data(), (size_t)cnt * sizeof(TileWork),
                                            cudaMemcpyHostToDevice, ring.StreamH2D()));
            if (!ring.BeginUpload(slot, inOff)) { ok = false; break; }
            cudaStream_t ks = ring.StreamKernel();
            const int TPB = 256;
            dim3 grid1((unsigned)(maxN_b > 0 ? (maxN_b + TPB - 1) / TPB : 1), (unsigned)cnt);
            if (maxN_b > 0) {
                if (kind == kKindFlow) init_flow_kernel<<<grid1, TPB, 0, ks>>>(dw);
                else init_numeric_kernel<<<grid1, TPB, 0, ks>>>(dw);
                SUCC_CUDA_CHECK(cudaGetLastError());
            }
            size_t ob = 0;
            if (phase == kPhaseRawWords) {
                for (int i = 0; i < cnt; i++) {
                    const TileWork& t = wk[(size_t)i];
                    if (t.n > 0)
                        SUCC_CUDA_CHECK(cudaMemcpyAsync(t.out, t.W, (size_t)t.n * sizeof(int32_t),
                                                        cudaMemcpyDeviceToDevice, ks));
                }
                ob = outOff2;
            } else {
                if (maxN_b > 0) {
                    const int np = ceil_log2_plus1(maxN_b);
                    for (int p = 0; p < np; p++) jump_kernel<<<grid1, TPB, 0, ks>>>(dw);
                    SUCC_CUDA_CHECK(cudaGetLastError());
                    size_t tb = tempBytes / (size_t)nslots;
                    void* temp = (unsigned char*)d_temp + (size_t)slot * tb;
                    for (int i = 0; i < cnt; i++) {
                        const TileWork& t = wk[(size_t)i];
                        if (t.n <= 0) continue;
                        SUCC_CUDA_CHECK(cub::DeviceScan::ExclusiveSum(temp, tb, t.flag, t.rank,
                                                                      t.n, ks));
                    }
                    if (payloadMode == kPayloadCompactId) {
                        rank_words_kernel<<<grid1, TPB, 0, ks>>>(dw);
                        SUCC_CUDA_CHECK(cudaGetLastError());
                    }
                }
                if (phase == kPhaseRing) {
                    dim3 gridR((unsigned)(maxR_b > 0 ? (maxR_b + TPB - 1) / TPB : 1), (unsigned)cnt);
                    extract_ring_kernel<<<gridR, TPB, 0, ks>>>(dw);
                    SUCC_CUDA_CHECK(cudaGetLastError());
                    ob = (size_t)cnt * PackedRing::Bytes(Rmax);
                } else {
                    if (maxN_b > 0) {
                        output_words_kernel<<<grid1, TPB, 0, ks>>>(dw);
                        SUCC_CUDA_CHECK(cudaGetLastError());
                    }
                    ob = outOff2;
                }
            }
            if (!ring.EndComputeAndDownload(slot, ob)) { ok = false; break; }
            launches++;
            next += cnt;
            slot = (slot + 1) % nslots;
        }
        for (int i = 0; i < nslots && ok; i++)
            if (!retire((slot + i) % nslots)) ok = false;
        return ok;
    }

    // The whole resolver: phase 1 -> stitch -> phase 2 (or raw words).
    bool Resolve(const PackFn& pack, bool rawWords, const WordSink2D& sink, ResolveStats2D& st) {
        const auto t0 = std::chrono::steady_clock::now();
        st.tiles = grid.numTiles;
        st.passes = passes;
        st.deviceBytesAllocated = DeviceBytes();
        st.hostPinnedBytes = (int64_t)ring.PinnedBytes();

        FrontierGraph2D fg;
        fg.NX = NX;
        fg.NY = NY;
        fg.grid = grid;
        fg.tiles.resize((size_t)grid.numTiles);
        bool ok = true;
        if (!rawWords) {
            const int Rmax_ = Rmax;
            ok = Run(pack, kPhaseRing,
                [&](int, int64_t first, int cnt, const std::vector<TileWork>& wk, unsigned char* h_out) {
                    for (int i = 0; i < cnt; i++) {
                        const TileWork& t = wk[(size_t)i];
                        TileRingRecord2D& tr = fg.tiles[(size_t)(first + i)];
                        unsigned char* packed = h_out + (size_t)i * PackedRing::Bytes(Rmax_);
                        tr.w = t.w; tr.h = t.h; tr.ox = t.ox; tr.oy = t.oy;
                        tr.R = t.R; tr.Rout = t.Rout;
                        tr.count = t.n > 0 ? *PackedRing::Count(packed, Rmax_) : 0;
                        tr.word.assign(PackedRing::Words(packed), PackedRing::Words(packed) + t.R);
                    }
                    return true;
                });
            st.phase1_ms = (float)std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - t0).count();
            if (!ok) return false;

            const auto t1 = std::chrono::steady_clock::now();
            int64_t base = 0;
            for (size_t t = 0; t < fg.tiles.size(); t++) {
                fg.tiles[t].base = base;
                base += fg.tiles[t].count;
                st.hostStitchBytes += (int64_t)fg.tiles[t].R * (int64_t)sizeof(int32_t) +
                                      (int64_t)fg.tiles[t].Rout * (int64_t)sizeof(int64_t);
            }
            st.terminals = base;
            if (!ResolveFrontier2D(fg, payloadMode, &st.stitchSteps)) {
                std::fprintf(stderr, "SuccessorResolver: malformed frontier graph\n");
                return false;
            }
            for (size_t t = 0; t < fg.tiles.size(); t++) {
                int64_t n = 0;
                for (int r = 0; r < fg.tiles[t].Rout; r++)
                    if (fg.tiles[t].resolved[(size_t)r] >= 0) n++;
                st.frontierSlots += n;
                if (n > st.frontierSlotsMaxTile) st.frontierSlotsMaxTile = n;
            }
            st.stitch_ms = (float)std::chrono::duration<double, std::milli>(
                               std::chrono::steady_clock::now() - t1).count();
        }

        const auto t2 = std::chrono::steady_clock::now();
        fgRun = &fg;
        SUCC_CUDA_CHECK(cudaMemset(d_events, 0, 16 * sizeof(int)));
        ok = Run(pack, rawWords ? kPhaseRawWords : kPhaseOutput,
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
        st.launches = launches;
        st.wall_ms = (float)std::chrono::duration<double, std::milli>(
                         std::chrono::steady_clock::now() - t0).count();
        return ok;
    }
};

static bool UploadStep(const Dgrad2DTables& tables, Step2D& step) {
    step = Step2D{};
    step.dx[tables.dir_px] = 1;  step.dy[tables.dir_px] = 0;
    step.dx[tables.dir_mx] = -1; step.dy[tables.dir_mx] = 0;
    step.dx[tables.dir_py] = 0;  step.dy[tables.dir_py] = 1;
    step.dx[tables.dir_my] = 0;  step.dy[tables.dir_my] = -1;
    SUCC_CUDA_CHECK(cudaMemcpyToSymbol(c_step, &step, sizeof(step)));
    return true;
}

} // namespace

bool Label2DFlowTerminalsTiled(const ByteSource2D& gradSource,
                               const TileGrid2D& grid,
                               const Dgrad2DTables& tables, bool ascending,
                               const WordSink2D& sink,
                               const ResolveOptions2D& options,
                               ResolveStats2D* stats) {
    if (grid.X < 2 || grid.Y < 2 || grid.T < 1 || grid.numTiles < 1) return false;
    Step2D step;
    if (!UploadStep(tables, step)) return false;
    const size_t capIn = (size_t)(2 * grid.T) * (size_t)(2 * grid.T);
    Engine eng(grid, kKindFlow, !ascending, options.slots < 1 ? 1 : options.slots,
               options.tilesPerLaunch < 1 ? 1 : options.tilesPerLaunch, capIn);
    eng.payloadMode = options.compactIds ? kPayloadCompactId : kPayloadLatticeIndex;
    if (!eng.Create()) return false;
    const bool quad = !ascending;
    PackFn pack = [&](int64_t, const TileDesc2D& d, unsigned char* h_in, const unsigned char* d_in,
                      size_t& inOff, TileWork& t) {
        const Box2D ob = OwnedCellBox2D(grid, d);
        if (!gradSource(ob.x0, ob.y0, ob.w, ob.h, h_in + inOff, ob.w)) return false;
        t.grad = d_in + inOff;
        t.gpitch = ob.w;
        t.quad = quad ? 1 : 0;
        inOff += ((size_t)ob.w * ob.h + 3) & ~(size_t)3;
        return true;
    };
    ResolveStats2D st;
    const bool ok = eng.Resolve(pack, false, sink, st);
    if (stats) *stats = st;
    return ok;
}

bool LabelNumericRegions2DTiled(const TileSource2D& values,
                                const ByteSource2D& gradSource,
                                const WordSource2D& baseLabels,
                                const TileGrid2D& grid,
                                const Dgrad2DTables& tables,
                                const int32_t* remap, int32_t nBase,
                                const float* sStar, int32_t nLiving,
                                const WordSink2D& sink,
                                const NumericOptions2D& options,
                                NumericStats2D* stats) {
    if (grid.X < 3 || grid.Y < 3 || grid.T < 1 || grid.numTiles < 1) return false;
    if (nBase < 0 || nLiving < 0) return false;
    Step2D step;
    if (!UploadStep(tables, step)) return false;

    int32_t* d_remap = nullptr;
    float* d_sstar = nullptr;
    if (nBase > 0) {
        SUCC_CUDA_CHECK(cudaMalloc(&d_remap, (size_t)nBase * sizeof(int32_t)));
        SUCC_CUDA_CHECK(cudaMemcpy(d_remap, remap, (size_t)nBase * sizeof(int32_t),
                                   cudaMemcpyHostToDevice));
    }
    if (nLiving > 0) {
        SUCC_CUDA_CHECK(cudaMalloc(&d_sstar, (size_t)nLiving * sizeof(float)));
        SUCC_CUDA_CHECK(cudaMemcpy(d_sstar, sStar, (size_t)nLiving * sizeof(float),
                                   cudaMemcpyHostToDevice));
    }
    numeric2d::Params np{};
    np.X = (int)grid.X;
    np.Y = (int)grid.Y;
    np.hopCells = options.hopCells < 1 ? 1 : options.hopCells;
    np.iterCap = options.iterCap < 1 ? 1 : options.iterCap;
    np.errThresh = options.errorThreshold;
    np.gradThresh = options.gradThreshold;
    np.remap = d_remap;
    np.sStar = d_sstar;
    np.nBase = nBase;
    np.nLiving = nLiving;
    for (int c = 0; c < 8; c++) { np.dx[c] = step.dx[c]; np.dy[c] = step.dy[c]; }
    cudaError_t e = cudaMemcpyToSymbol(c_np, &np, sizeof(np));
    if (e != cudaSuccess) { cudaFree(d_remap); cudaFree(d_sstar); return false; }

    // per-tile input: labels (owned+1 halo, int32), values (owned+2 halo,
    // float), gradient bytes (owned cell box), each 4-byte aligned
    const size_t capIn = (size_t)(grid.T + 2) * (grid.T + 2) * sizeof(int32_t) +
                         (size_t)(grid.T + 4) * (grid.T + 4) * sizeof(float) +
                         (((size_t)(2 * grid.T) * (2 * grid.T) + 3) & ~(size_t)3);
    Engine eng(grid, kKindNumeric, false, options.slots < 1 ? 1 : options.slots,
               options.tilesPerLaunch < 1 ? 1 : options.tilesPerLaunch, capIn);
    eng.payloadMode = kPayloadDirect;
    if (!eng.Create()) { cudaFree(d_remap); cudaFree(d_sstar); return false; }

    PackFn pack = [&](int64_t, const TileDesc2D& d, unsigned char* h_in, const unsigned char* d_in,
                      size_t& inOff, TileWork& t) {
        numeric2d::TileIn& n = t.nin;
        n.ox = (int)d.ox; n.oy = (int)d.oy; n.tw = d.tw; n.th = d.th;
        // labels, halo 1
        {
            const int64_t x0 = d.ox - 1 < 0 ? 0 : d.ox - 1, y0 = d.oy - 1 < 0 ? 0 : d.oy - 1;
            const int64_t x1 = d.ox + d.tw + 1 > grid.X ? grid.X : d.ox + d.tw + 1;
            const int64_t y1 = d.oy + d.th + 1 > grid.Y ? grid.Y : d.oy + d.th + 1;
            const int w = (int)(x1 - x0), h = (int)(y1 - y0);
            if (!baseLabels(x0, y0, w, h, (int32_t*)(h_in + inOff), w)) return false;
            n.base = (const int32_t*)(d_in + inOff);
            n.bpitch = w; n.bx0 = (int)x0; n.by0 = (int)y0; n.bw = w; n.bh = h;
            inOff += (size_t)w * h * sizeof(int32_t);
        }
        // values, halo 2
        {
            const int64_t x0 = d.ox - 2 < 0 ? 0 : d.ox - 2, y0 = d.oy - 2 < 0 ? 0 : d.oy - 2;
            const int64_t x1 = d.ox + d.tw + 2 > grid.X ? grid.X : d.ox + d.tw + 2;
            const int64_t y1 = d.oy + d.th + 2 > grid.Y ? grid.Y : d.oy + d.th + 2;
            const int w = (int)(x1 - x0), h = (int)(y1 - y0);
            if (!values(x0, y0, w, h, (float*)(h_in + inOff), w)) return false;
            n.vals = (const float*)(d_in + inOff);
            n.vpitch = w; n.sx0 = (int)x0; n.sy0 = (int)y0; n.svx = w; n.svy = h;
            inOff += (size_t)w * h * sizeof(float);
        }
        // gradient bytes, owned cell box
        {
            const Box2D ob = OwnedCellBox2D(grid, d);
            if (!gradSource(ob.x0, ob.y0, ob.w, ob.h, h_in + inOff, ob.w)) return false;
            n.grad = d_in + inOff;
            n.gpitch = ob.w;
            inOff += ((size_t)ob.w * ob.h + 3) & ~(size_t)3;
        }
        return true;
    };
    NumericStats2D st;
    bool ok = eng.Resolve(pack, options.debugRawWords, sink, st.resolve);
    if (ok) {
        int h_ev[16];
        e = cudaMemcpy(h_ev, eng.d_events, sizeof(h_ev), cudaMemcpyDeviceToHost);
        if (e != cudaSuccess) ok = false;
        else for (int i = 0; i < 16; i++) st.events[i] = h_ev[i];
    }
    st.vertices = grid.X * grid.Y;
    cudaFree(d_remap);
    cudaFree(d_sstar);
    if (stats) *stats = st;
    return ok;
}

} // namespace gpu
} // namespace GInt
