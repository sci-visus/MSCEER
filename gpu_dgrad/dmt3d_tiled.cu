// Tiled (blocked) GPU drivers for the 3D DMT rules: PLS / PLS-BR and the
// multi-pass Shivashankar-Natarajan DMF. Host API in dgrad_gpu_api.h.
//
// Tiling follows the measured halos of the blocking notes:
//   - ownership: cell c belongs to the tile owning vertex floor(c/2) on every
//     axis (tile3d.h, the same convention as the 2D tiles and the 3D Robins
//     path), so the owned cells touch the vertices [o, o+t] per axis;
//   - staged vertices = [o - lo, o + t + 1 + hi): lo/hi = 1/1 for PLS, PLS-BR
//     and first-pass DMF, 2/1 for DMF with the second pass (a low-face vertex's
//     second pass looks at the first-pass choice of its -1 neighbour, which
//     compares keys reaching vertex o-2).
// Each tile is a cubical complex of its own (Cubical3DT<int>, the staged box)
// and runs every pass of its rule on local data: DMF top facets d = 1..dmax,
// algorithm 1 for d = 0..dmax-1, algorithm 2 for d = 0..dmax-1, as separate
// kernels on one stream (the kernel boundary is the barrier between passes).
// Cells outside the owned box are recomputed redundantly and never leave the
// device; counters (critical cells, guards, claims, orphans, unassigned) are
// counted on owned cells only.
//
// Indexing: int64 on the host (tile origins, boxes); the kernels see int32
// tile-local coordinates and pitches only, and carry coordinates instead of
// ids (3D launches over one parity class of cells at a time), so there are no
// id -> coordinate divisions.

#include "dgrad_gpu_api.h"
#include "gpu_transfer.h"
#include "tile3d.h"

#include "dmt_core/dmt_hd.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_star27_tables.h"
#include "dmt_core/dmt_lowerstar_pq.h"
#include "dmt_core/dmt_dmf_cell.h"
#include "dmt_core/dmt_dmf_keycache.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <vector>

namespace GInt {
namespace gpu {

namespace {

#define DMTT_CHECK(call)                                                               \
    do {                                                                               \
        cudaError_t e_ = (call);                                                       \
        if (e_ != cudaSuccess) {                                                       \
            std::fprintf(stderr, "dmt3d_tiled %s:%d %s: %s\n", __FILE__, __LINE__, #call, \
                         cudaGetErrorString(e_));                                      \
            return false;                                                              \
        }                                                                              \
    } while (0)

typedef dmt::Cubical3DT<int> LGrid;   // tile-local cubical complex

enum : int {
    kCrit0 = 0,          // 0..3 critical cells per dimension (owned)
    kUnassigned = 4,
    kOrphans = 5,
    kGuard1 = 6, kAlg1 = 7, kGuard2 = 8, kAlg2 = 9,
    kNumCnt = 10
};

// Everything a kernel needs about one tile, in local coordinates: cell q of
// the staged complex has global doubled coordinates q + 2*s0.
struct TileGeom {
    LGrid L;              // staged box as a complex
    int sx0, sy0, sz0;    // global vertex coords of the staged origin
    int X, Y, Z;          // volume size (vertices)
    int olo[3], ohi[3];   // owned cell box, local cell coords, [lo, hi)
    int clo[3], chi[3];   // computed vertex box (lower star), local vertex coords
    int opy, opz;         // pitches of the owned output box
};

__device__ __forceinline__ bool OwnedLocal(const TileGeom& t, int x, int y, int z) {
    return x >= t.olo[0] && x < t.ohi[0] && y >= t.olo[1] && y < t.ohi[1] && z >= t.olo[2] &&
           z < t.ohi[2];
}

__device__ __forceinline__ void BlockCount(unsigned* sh, int slot, unsigned v) {
    if (v) atomicAdd(&sh[slot], v);
}

// One parity class of cells: q = 2*t + par on each axis, t in [0, ext).
struct ParClass {
    int par[3];
    int ext[3];
};

ParClass MakeClass(const LGrid& L, int p) {
    ParClass c;
    for (int i = 0; i < 3; i++) {
        c.par[i] = (p >> i) & 1;
        c.ext[i] = c.par[i] ? L.n[i] - 1 : L.n[i];
    }
    return c;
}

const dim3 kBlk(32, 4, 2);

bool GridFor(const int ext[3], dim3& g) {
    const long long gx = (ext[0] + kBlk.x - 1) / kBlk.x, gy = (ext[1] + kBlk.y - 1) / kBlk.y,
                    gz = (ext[2] + kBlk.z - 1) / kBlk.z;
    if (gx < 1 || gy < 1 || gz < 1 || gy > 65535 || gz > 65535) return false;
    g = dim3((unsigned)gx, (unsigned)gy, (unsigned)gz);
    return true;
}

__device__ __forceinline__ bool ClassCellAt(const ParClass& K, int& x, int& y, int& z) {
    const int tx = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    const int ty = (int)(blockIdx.y * blockDim.y + threadIdx.y);
    const int tz = (int)(blockIdx.z * blockDim.z + threadIdx.z);
    if (tx >= K.ext[0] || ty >= K.ext[1] || tz >= K.ext[2]) return false;
    x = 2 * tx + K.par[0];
    y = 2 * ty + K.par[1];
    z = 2 * tz + K.par[2];
    return true;
}

// ---------------------------------------------------------------------------
// DMF kernels
// ---------------------------------------------------------------------------

__global__ void DMFTopTileKernel(TileGeom t, ParClass K, const uint32_t* __restrict__ ranks,
                                 uint8_t* top) {
    int x, y, z;
    if (!ClassCellAt(K, x, y, z)) return;
    uint8_t g0 = 7, h2 = 7;
    dmt::DMFTopFacetsFromTopQ(t.L, ranks, top, x, y, z, g0, h2);
    top[t.L.Id(x, y, z)] = (uint8_t)(g0 | (h2 << 4));
}

__global__ void DMFPassTileKernel(TileGeom t, ParClass K, const uint32_t* __restrict__ ranks,
                                  const uint8_t* __restrict__ top, int alg, uint8_t* grad,
                                  uint8_t* stage, unsigned long long* counters) {
    __shared__ unsigned sh[2];
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) { sh[0] = 0; sh[1] = 0; }
    __syncthreads();
    int x, y, z;
    if (ClassCellAt(K, x, y, z)) {
        const int a = t.L.Id(x, y, z);
        if (!(alg == 2 && grad[a] != 0)) {
            const int dir = dmt::DMFLowestSelectingCofacetFromTopQ(t.L, ranks, top, x, y, z, alg - 1);
            if (dir >= 0) {
                const int b = t.L.Neighbor(a, dir);
                const bool own = OwnedLocal(t, x, y, z);
                if (grad[a] != 0 || grad[b] != 0) {
                    if (own) atomicAdd(&sh[0], 1u);
                } else {
                    grad[a] = dmt::MakeGradByte((uint8_t)dir);
                    grad[b] = dmt::MakeGradByte((uint8_t)(dir ^ 1));
                    if (stage) { stage[a] = (uint8_t)alg; stage[b] = (uint8_t)alg; }
                    if (own) atomicAdd(&sh[1], 1u);
                }
            }
        }
    }
    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
        if (sh[0]) atomicAdd(&counters[alg == 1 ? kGuard1 : kGuard2], (unsigned long long)sh[0]);
        if (sh[1]) atomicAdd(&counters[alg == 1 ? kAlg1 : kAlg2], (unsigned long long)sh[1]);
    }
}

// Owned cells -> output box: unpaired (0) cells become critical; critical
// cells counted per dimension.
__global__ void DMFExtractKernel(TileGeom t, const uint8_t* __restrict__ grad,
                                 const uint8_t* __restrict__ stage, uint8_t* __restrict__ outGrad,
                                 uint8_t* __restrict__ outStage, unsigned long long* counters) {
    __shared__ unsigned sh[4];
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    if (lt < 4) sh[lt] = 0;
    __syncthreads();
    const int ox = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    const int oy = (int)(blockIdx.y * blockDim.y + threadIdx.y);
    const int oz = (int)(blockIdx.z * blockDim.z + threadIdx.z);
    const int w = t.ohi[0] - t.olo[0], h = t.ohi[1] - t.olo[1], d = t.ohi[2] - t.olo[2];
    if (ox < w && oy < h && oz < d) {
        const int x = t.olo[0] + ox, y = t.olo[1] + oy, z = t.olo[2] + oz;
        const int c = t.L.Id(x, y, z);
        uint8_t b = grad[c];
        if (b == 0) {
            b = dmt::MakeGradByte(dmt::kPairCritical);
            atomicAdd(&sh[(x & 1) + (y & 1) + (z & 1)], 1u);
        }
        const int o = oz * t.opz + oy * t.opy + ox;
        outGrad[o] = b;
        if (outStage) outStage[o] = stage[c];
    }
    __syncthreads();
    if (lt < 4 && sh[lt]) atomicAdd(&counters[kCrit0 + lt], (unsigned long long)sh[lt]);
}

// ---------------------------------------------------------------------------
// Lower-star kernel (PLS / PLS-BR)
// ---------------------------------------------------------------------------


// The star tables are read at thread-dependent slots (facet / cofacet masks of
// the cell being collapsed); from kernel-parameter space that serializes, so
// every lower-star kernel copies them to shared memory first.
__device__ __forceinline__ void LoadStarTablesShared(const dmt::Star27Tables& T, dmt::Star27Tables& sT) {
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    const int nt = blockDim.x * blockDim.y * blockDim.z;
    const unsigned char* src = reinterpret_cast<const unsigned char*>(&T);
    unsigned char* dst = reinterpret_cast<unsigned char*>(&sT);
    for (int i = lt; i < (int)sizeof(dmt::Star27Tables); i += nt) dst[i] = src[i];
    __syncthreads();
}

// 64 registers (4 blocks of 256 per SM) measured fastest: 20 ms vs 25 ms at 84
// registers for 256^3 PLS.
__global__ void __launch_bounds__(256, 4) LowerStarTileKernel(dmt::Star27Tables Tp, TileGeom t, const uint32_t* __restrict__ ranks,
                                    int boundary_restricted, uint8_t* __restrict__ outGrad,
                                    unsigned long long* counters) {
    __shared__ dmt::Star27Tables T;
    LoadStarTablesShared(Tp, T);
    // local vertex of the computed box
    const int vx = t.clo[0] + (int)(blockIdx.x * blockDim.x + threadIdx.x);
    const int vy = t.clo[1] + (int)(blockIdx.y * blockDim.y + threadIdx.y);
    const int vz = t.clo[2] + (int)(blockIdx.z * blockDim.z + threadIdx.z);
    if (vx >= t.chi[0] || vy >= t.chi[1] || vz >= t.chi[2]) return;
    // the staged box holds every in-domain neighbour of a computed vertex, and
    // a neighbour outside the staged box is outside the volume
    uint32_t key[27], vvalid = 0;
#pragma unroll
    for (int dz = -1; dz <= 1; dz++)
#pragma unroll
        for (int dy = -1; dy <= 1; dy++)
#pragma unroll
            for (int dx = -1; dx <= 1; dx++) {
                const int j = dmt::StarSlot(dx, dy, dz);
                const int x = vx + dx, y = vy + dy, z = vz + dz;
                if (x < 0 || y < 0 || z < 0 || x >= t.L.n[0] || y >= t.L.n[1] || z >= t.L.n[2]) {
                    key[j] = 0;
                    continue;
                }
                key[j] = ranks[x + t.L.n[0] * (y + t.L.n[1] * z)];
                vvalid |= 1u << j;
            }
    uint32_t bm = 0;
    if (boundary_restricted) {
        // StarBoundaryMask with global extents
        const int g[3] = {2 * (vx + t.sx0), 2 * (vy + t.sy0), 2 * (vz + t.sz0)};
        const int top[3] = {2 * t.X - 2, 2 * t.Y - 2, 2 * t.Z - 2};
#pragma unroll
        for (int s = 0; s < 27; s++) {
            const int ox = s % 3 - 1, oy = (s / 3) % 3 - 1, oz = s / 9 - 1;
            const int q0 = g[0] + ox, q1 = g[1] + oy, q2 = g[2] + oz;
            if (q0 == 0 || q0 == top[0] || q1 == 0 || q1 == top[1] || q2 == 0 || q2 == top[2]) bm |= 1u << s;
        }
    }
    uint32_t L;
    dmt::StarCodeRecorder rec;
    dmt::LowerStarCodes(T, key, vvalid, boundary_restricted != 0, bm, L, rec);
    // the star belongs to the vertex: count its orphans if the tile owns it
    const int cvx = 2 * vx, cvy = 2 * vy, cvz = 2 * vz;
    if (rec.n_orphans && OwnedLocal(t, cvx, cvy, cvz))
        atomicAdd(&counters[kOrphans], (unsigned long long)rec.n_orphans);
    uint32_t m = rec.done & L;   // the cells PLSStar would write
    while (m) {
        const int s = __ffs((int)m) - 1;
        m &= m - 1;
        const int x = cvx + s % 3 - 1, y = cvy + (s / 3) % 3 - 1, z = cvz + s / 9 - 1;
        if (!OwnedLocal(t, x, y, z)) continue;
        outGrad[(z - t.olo[2]) * t.opz + (y - t.olo[1]) * t.opy + (x - t.olo[0])] =
            dmt::MakeGradByte((uint8_t)rec.Code(s));
    }
}

// Census of an owned output box: unassigned bytes and critical cells per
// dimension (cell dimension from the global parity = local parity).
__global__ void LowerStarCensusKernel(TileGeom t, const uint8_t* __restrict__ outGrad,
                                      unsigned long long* counters) {
    __shared__ unsigned sh[5];
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    if (lt < 5) sh[lt] = 0;
    __syncthreads();
    const int ox = (int)(blockIdx.x * blockDim.x + threadIdx.x);
    const int oy = (int)(blockIdx.y * blockDim.y + threadIdx.y);
    const int oz = (int)(blockIdx.z * blockDim.z + threadIdx.z);
    if (ox < t.ohi[0] - t.olo[0] && oy < t.ohi[1] - t.olo[1] && oz < t.ohi[2] - t.olo[2]) {
        const uint8_t b = outGrad[oz * t.opz + oy * t.opy + ox];
        const int dim = ((t.olo[0] + ox) & 1) + ((t.olo[1] + oy) & 1) + ((t.olo[2] + oz) & 1);
        if (!dmt::GradAssigned(b)) atomicAdd(&sh[4], 1u);
        else if (dmt::GradCritical(b)) atomicAdd(&sh[dim], 1u);
    }
    __syncthreads();
    if (lt < 5 && sh[lt]) atomicAdd(&counters[lt == 4 ? kUnassigned : kCrit0 + lt], (unsigned long long)sh[lt]);
}

double MsSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

inline size_t Align16(size_t b) { return (b + 15) & ~(size_t)15; }

struct DevMem {
    void* p = nullptr;
    ~DevMem() { if (p) cudaFree(p); }
    cudaError_t alloc(size_t n) { return cudaMalloc(&p, n ? n : 1); }
};

enum class Rule { PLS, PLS_BR, DMF };

struct TilePlan {
    TileGeom g;
    Box3D owned;
};

// Staged box for halo (lo, hi) beyond [o, o+t].
Box3D StagedBox(const TileGrid3D& G, const TileDesc3D& d, int lo, int hi) {
    Box3D b;
    auto ax = [](int64_t o, int t, int64_t N, int lo, int hi, int64_t& x0, int& w) {
        int64_t a = o - lo, e = o + t + 1 + hi;
        if (a < 0) a = 0;
        if (e > N) e = N;
        x0 = a;
        w = (int)(e - a);
    };
    ax(d.ox, d.tw, G.X, lo, hi, b.x0, b.w);
    ax(d.oy, d.th, G.Y, lo, hi, b.y0, b.h);
    ax(d.oz, d.td, G.Z, lo, hi, b.z0, b.d);
    return b;
}

TilePlan PlanTile(const TileGrid3D& G, const TileDesc3D& d, const Box3D& sb) {
    TilePlan p;
    TileGeom& g = p.g;
    std::memset(&g, 0, sizeof(g));
    g.L = LGrid::Make(sb.w, sb.h, sb.d);
    g.sx0 = (int)sb.x0; g.sy0 = (int)sb.y0; g.sz0 = (int)sb.z0;
    g.X = (int)G.X; g.Y = (int)G.Y; g.Z = (int)G.Z;
    p.owned = OwnedCellBox3D(G, d);
    const int64_t s0[3] = {sb.x0, sb.y0, sb.z0};
    const int64_t ob0[3] = {p.owned.x0, p.owned.y0, p.owned.z0};
    const int ow[3] = {p.owned.w, p.owned.h, p.owned.d};
    const Box3D cb = ComputedVertexBox3D(G, d);
    const int64_t cb0[3] = {cb.x0, cb.y0, cb.z0};
    const int cw[3] = {cb.w, cb.h, cb.d};
    for (int i = 0; i < 3; i++) {
        g.olo[i] = (int)(ob0[i] - 2 * s0[i]);
        g.ohi[i] = g.olo[i] + ow[i];
        g.clo[i] = (int)(cb0[i] - s0[i]);
        g.chi[i] = g.clo[i] + cw[i];
    }
    g.opy = p.owned.w;
    g.opz = p.owned.w * p.owned.h;
    return p;
}

int MaxDimOf(const TileGrid3D& G) {
    return (G.X > 1) + (G.Y > 1) + (G.Z > 1);
}

bool TiledDmt(Rule rule, bool second_pass, const RankSource3D& source, const TileGrid3D& grid,
              const ByteSink3D& gradSink, const ByteSink3D* stageSink,
              const DmtTiledOptions3D& opt, unsigned long long cntOut[kNumCnt], TileStats3D* tstats) {
    if (grid.X < 1 || grid.Y < 1 || grid.Z < 1 || grid.numTiles < 1) return false;
    if (grid.X >= (1 << 30) || grid.Y >= (1 << 30) || grid.Z >= (1 << 30)) return false;
    const int defLo = (rule == Rule::DMF && second_pass) ? 2 : 1;
    const int lo = opt.haloLo < 0 ? defLo : opt.haloLo;
    const int hi = opt.haloHi < 0 ? 1 : opt.haloHi;

    // capacities from the largest tile shape
    auto cap = [&](int T, int64_t N) { return std::min<int64_t>((int64_t)T + 1 + lo + hi, N); };
    const int64_t sx = cap(grid.Tx, grid.X), sy = cap(grid.Ty, grid.Y), sz = cap(grid.Tz, grid.Z);
    const int64_t stagedMax = sx * sy * sz;
    const int64_t localCells = (2 * sx - 1) * (2 * sy - 1) * (2 * sz - 1);
    const int64_t ownedMax = std::min<int64_t>(2 * (int64_t)grid.Tx, 2 * grid.X - 1) *
                             std::min<int64_t>(2 * (int64_t)grid.Ty, 2 * grid.Y - 1) *
                             std::min<int64_t>(2 * (int64_t)grid.Tz, 2 * grid.Z - 1);
    if (localCells >= ((int64_t)1 << 31)) {
        std::fprintf(stderr, "dmt3d_tiled: tile too large for 32-bit tile-local indexing\n");
        return false;
    }
    const bool wantStage = rule == Rule::DMF && stageSink != nullptr;
    const size_t capIn = Align16((size_t)stagedMax * sizeof(uint32_t));
    const size_t outPlane = Align16((size_t)ownedMax);
    const size_t capOut = outPlane * (wantStage ? 2 : 1);

    int nslots = opt.slots < 1 ? 1 : opt.slots;
    if ((int64_t)nslots > grid.numTiles) nslots = (int)grid.numTiles;

    const auto t_wall0 = std::chrono::steady_clock::now();
    TileStats3D st;
    st.tiles = grid.numTiles;

    TransferRing ring;
    if (!ring.Create(nslots, capIn, capOut)) return false;
    // local complex scratch: kernels of all slots run on one stream, so one
    // set is enough
    DevMem d_top, d_grad, d_stage, d_cnt;
    if (rule == Rule::DMF) {
        DMTT_CHECK(d_top.alloc((size_t)localCells));
        DMTT_CHECK(d_grad.alloc((size_t)localCells));
        if (wantStage) DMTT_CHECK(d_stage.alloc((size_t)localCells));
    }
    DMTT_CHECK(d_cnt.alloc(sizeof(unsigned long long) * kNumCnt));
    DMTT_CHECK(cudaMemsetAsync(d_cnt.p, 0, sizeof(unsigned long long) * kNumCnt, ring.StreamKernel()));
    st.hostPinnedBytes = (int64_t)ring.PinnedBytes();
    st.deviceBytesAllocated = (int64_t)ring.DeviceBytes() +
                              (rule == Rule::DMF ? (int64_t)localCells * (wantStage ? 3 : 2) : 0) +
                              (int64_t)(sizeof(unsigned long long) * kNumCnt);

    dmt::Star27Tables T27;
    dmt::BuildStar27TablesCanonical(T27);
    const int dmax = MaxDimOf(grid);
    unsigned long long* cnt = (unsigned long long*)d_cnt.p;
    cudaStream_t ks = ring.StreamKernel();

    std::vector<TilePlan> plans((size_t)nslots);
    std::vector<int> busy((size_t)nslots, 0);

    auto retire = [&](int slot) -> bool {
        if (!busy[(size_t)slot]) return true;
        float h2d = 0.f, kern = 0.f, d2h = 0.f;
        if (!ring.Retire(slot, &h2d, &kern, &d2h)) return false;
        st.h2d_ms += h2d;
        st.kernel_ms += kern;
        st.d2h_ms += d2h;
        const auto t0 = std::chrono::steady_clock::now();
        const Box3D& o = plans[(size_t)slot].owned;
        const uint8_t* out = (const uint8_t*)ring.Slot(slot).h_out;
        if (!gradSink(o, out, o.w, (int64_t)o.w * o.h)) return false;
        if (wantStage && !(*stageSink)(o, out + outPlane, o.w, (int64_t)o.w * o.h)) return false;
        st.sink_ms += (float)MsSince(t0);
        busy[(size_t)slot] = 0;
        return true;
    };

    bool ok = true;
    int slot = 0;
    for (int64_t tile = 0; tile < grid.numTiles && ok; tile++) {
        if (!retire(slot)) { ok = false; break; }
        const TileDesc3D d = TileAt3D(grid, tile);
        const Box3D sb = StagedBox(grid, d, lo, hi);
        TilePlan& P = plans[(size_t)slot];
        P = PlanTile(grid, d, sb);
        if (sb.Count() > stagedMax || P.owned.Count() > ownedMax) {
            std::fprintf(stderr, "dmt3d_tiled: internal error, tile exceeds its buffer\n");
            ok = false;
            break;
        }
        TransferSlot& s = ring.Slot(slot);
        const auto t_src = std::chrono::steady_clock::now();
        if (!source(sb, (uint32_t*)s.h_in, sb.w, (int64_t)sb.w * sb.h)) { ok = false; break; }
        st.source_ms += (float)MsSince(t_src);
        if (!ring.BeginUpload(slot, (size_t)sb.Count() * sizeof(uint32_t))) { ok = false; break; }

        const uint32_t* ranks = (const uint32_t*)s.d_in;
        uint8_t* outG = (uint8_t*)s.d_out;
        uint8_t* outS = wantStage ? outG + outPlane : nullptr;
        const TileGeom& g = P.g;
        const int oext[3] = {P.owned.w, P.owned.h, P.owned.d};
        dim3 go;
        if (!GridFor(oext, go)) { ok = false; break; }
        if (rule == Rule::DMF) {
            const int64_t nloc = (int64_t)g.L.NumCells();
            if (cudaMemsetAsync(d_top.p, 0x77, (size_t)nloc, ks) != cudaSuccess ||
                cudaMemsetAsync(d_grad.p, 0, (size_t)nloc, ks) != cudaSuccess ||
                (wantStage && cudaMemsetAsync(d_stage.p, 0, (size_t)nloc, ks) != cudaSuccess)) {
                ok = false;
                break;
            }
            uint8_t* top = (uint8_t*)d_top.p;
            uint8_t* gr = (uint8_t*)d_grad.p;
            uint8_t* sg = wantStage ? (uint8_t*)d_stage.p : nullptr;
            auto forDim = [&](int dim, auto launch) -> bool {
                for (int p = 0; p < 8; p++) {
                    if (((p & 1) + ((p >> 1) & 1) + ((p >> 2) & 1)) != dim) continue;
                    const ParClass K = MakeClass(g.L, p);
                    if (K.ext[0] < 1 || K.ext[1] < 1 || K.ext[2] < 1) continue;
                    dim3 gk;
                    if (!GridFor(K.ext, gk)) return false;
                    launch(K, gk);
                    if (cudaGetLastError() != cudaSuccess) return false;
                }
                return true;
            };
            for (int dim = 1; dim <= dmax && ok; dim++)
                ok = forDim(dim, [&](const ParClass& K, dim3 gk) {
                    DMFTopTileKernel<<<gk, kBlk, 0, ks>>>(g, K, ranks, top);
                });
            for (int alg = 1; alg <= (second_pass ? 2 : 1) && ok; alg++)
                for (int dim = 0; dim < dmax && ok; dim++)
                    ok = forDim(dim, [&](const ParClass& K, dim3 gk) {
                        DMFPassTileKernel<<<gk, kBlk, 0, ks>>>(g, K, ranks, top, alg, gr, sg, cnt);
                    });
            if (!ok) break;
            DMFExtractKernel<<<go, kBlk, 0, ks>>>(g, gr, sg, outG, outS, cnt);
            if (cudaGetLastError() != cudaSuccess) { ok = false; break; }
        } else {
            if (cudaMemsetAsync(outG, 0, (size_t)P.owned.Count(), ks) != cudaSuccess) { ok = false; break; }
            const int cext[3] = {g.chi[0] - g.clo[0], g.chi[1] - g.clo[1], g.chi[2] - g.clo[2]};
            dim3 gv;
            if (!GridFor(cext, gv)) { ok = false; break; }
            LowerStarTileKernel<<<gv, kBlk, 0, ks>>>(T27, g, ranks, rule == Rule::PLS_BR ? 1 : 0, outG, cnt);
            LowerStarCensusKernel<<<go, kBlk, 0, ks>>>(g, outG, cnt);
            if (cudaGetLastError() != cudaSuccess) { ok = false; break; }
        }
        const size_t dbytes = wantStage ? outPlane + (size_t)P.owned.Count() : (size_t)P.owned.Count();
        if (!ring.EndComputeAndDownload(slot, dbytes)) { ok = false; break; }
        busy[(size_t)slot] = 1;
        st.launches++;
        st.bytesH2D += (int64_t)sb.Count() * 4;
        st.bytesD2H += (int64_t)P.owned.Count() * (wantStage ? 2 : 1);
        size_t freeB = 0, totB = 0;
        if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) {
            if (st.deviceFreeMin < 0 || (int64_t)freeB < st.deviceFreeMin) st.deviceFreeMin = (int64_t)freeB;
            if ((int64_t)freeB > st.deviceFreeMax) st.deviceFreeMax = (int64_t)freeB;
        }
        slot = (slot + 1) % nslots;
    }
    for (int i = 0; i < nslots && ok; i++)
        if (!retire((slot + i) % nslots)) ok = false;
    if (ok) {
        DMTT_CHECK(cudaMemcpy(cntOut, d_cnt.p, sizeof(unsigned long long) * kNumCnt, cudaMemcpyDeviceToHost));
    } else {
        cudaDeviceSynchronize();
    }
    ring.Destroy();
    st.wall_ms = (float)MsSince(t_wall0);
    if (tstats) *tstats = st;
    return ok;
}

} // namespace

bool ComputeLowerStarGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
                                     DmtLowerStarRule rule, const ByteSink3D& gradSink,
                                     const DmtTiledOptions3D& options, DmtLowerStarStats* stats,
                                     TileStats3D* tstats) {
    unsigned long long c[kNumCnt] = {};
    const bool ok = TiledDmt(rule == DmtLowerStarRule::PLS_BR ? Rule::PLS_BR : Rule::PLS, false, source,
                             grid, gradSink, nullptr, options, c, tstats);
    if (ok && stats) {
        for (int d = 0; d < 4; d++) stats->critical[d] = (int64_t)c[kCrit0 + d];
        stats->orphans = (int64_t)c[kOrphans];
        stats->unassigned = (int64_t)c[kUnassigned];
    }
    return ok;
}

bool ComputeDMFGradient3DTiled(const RankSource3D& source, const TileGrid3D& grid,
                               bool second_pass, const ByteSink3D& gradSink,
                               const ByteSink3D* stageSink, const DmtTiledOptions3D& options,
                               DmtDMFStats* stats, TileStats3D* tstats) {
    unsigned long long c[kNumCnt] = {};
    const bool ok = TiledDmt(Rule::DMF, second_pass, source, grid, gradSink, stageSink, options, c, tstats);
    if (ok && stats) {
        stats->guard1 = (int64_t)c[kGuard1];
        stats->guard2 = (int64_t)c[kGuard2];
        stats->alg1 = (int64_t)c[kAlg1];
        stats->alg2 = (int64_t)c[kAlg2];
        for (int d = 0; d < 4; d++) stats->critical[d] = (int64_t)c[kCrit0 + d];
    }
    return ok;
}

} // namespace gpu
} // namespace GInt
