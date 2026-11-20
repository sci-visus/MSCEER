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
#include "dgrad3d_internal.cuh"
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
    kWritten = 10,       // lower star: owned cells written (unassigned = owned - written)
    kNumCnt = 11
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
    int vs[3];            // vertex strides of the rank buffer: dense staged box
                          // {1, w, w*h}, or the full-XY slab {1, X, X*Y}
};

// Cell strides of the tile-local complex, vertex strides of the rank buffer.
__device__ __forceinline__ dmt::DMFStrides<int> TileStrides(const TileGeom& t) {
    dmt::DMFStrides<int> S = dmt::MakeDMFStrides(t.L);
    S.vs[0] = t.vs[0];
    S.vs[1] = t.vs[1];
    S.vs[2] = t.vs[2];
    return S;
}

__device__ __forceinline__ bool OwnedLocal(const TileGeom& t, int x, int y, int z) {
    return x >= t.olo[0] && x < t.ohi[0] && y >= t.olo[1] && y < t.ohi[1] && z >= t.olo[2] &&
           z < t.ohi[2];
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

// Kernels are templated on the parity pattern P of their cells (and the
// selector of the pass), so the key walks of dmt_dmf_keycache.h unroll into
// registers; one launch covers one parity class.
template <int P>
__global__ void DMFTopTileKernel(TileGeom t, ParClass K, const uint32_t* __restrict__ ranks,
                                 uint8_t* top) {
    int x, y, z;
    if (!ClassCellAt(K, x, y, z)) return;
    const dmt::DMFStrides<int> S = TileStrides(t);
    const int cid = t.L.Id(x, y, z);
    const int vb = (x >> 1) * t.vs[0] + (y >> 1) * t.vs[1] + (z >> 1) * t.vs[2];
    uint8_t g0 = 7, h2 = 7;
    dmt::DMFTopFacetsP<P>(S, ranks, top, cid, vb, g0, h2);
    top[cid] = (uint8_t)(g0 | (h2 << 4));
}

template <int P, int W>
__global__ void DMFPassTileKernel(TileGeom t, ParClass K, const uint32_t* __restrict__ ranks,
                                  const uint8_t* __restrict__ top, uint8_t* grad,
                                  uint8_t* stage, unsigned long long* counters) {
    constexpr int alg = W + 1;
    __shared__ unsigned sh[2];
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) { sh[0] = 0; sh[1] = 0; }
    __syncthreads();
    int x, y, z;
    bool guard = false, claim = false;
    if (ClassCellAt(K, x, y, z)) {
        const int a = t.L.Id(x, y, z);
        if (!(alg == 2 && grad[a] != 0)) {
            const dmt::DMFStrides<int> S = TileStrides(t);
            const int vb = (x >> 1) * t.vs[0] + (y >> 1) * t.vs[1] + (z >> 1) * t.vs[2];
            const int q[3] = {x, y, z};
            const int dir = dmt::DMFLowestSelectingP<P, W>(t.L, S, ranks, top, a, vb, q);
            if (dir >= 0) {
                const int st = dmt::DMFSelAxis(S.cs, dir >> 1);
                const int b = (dir & 1) ? a - st : a + st;
                const bool own = OwnedLocal(t, x, y, z);
                if (grad[a] != 0 || grad[b] != 0) {
                    guard = own;
                } else {
                    grad[a] = dmt::MakeGradByte((uint8_t)dir);
                    grad[b] = dmt::MakeGradByte((uint8_t)(dir ^ 1));
                    if (stage) { stage[a] = (uint8_t)alg; stage[b] = (uint8_t)alg; }
                    claim = own;
                }
            }
        }
    }
    if (guard) atomicAdd(&sh[0], 1u);
    if (claim) atomicAdd(&sh[1], 1u);
    __syncthreads();
    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0) {
        if (sh[0]) atomicAdd(&counters[alg == 1 ? kGuard1 : kGuard2], (unsigned long long)sh[0]);
        if (sh[1]) atomicAdd(&counters[alg == 1 ? kAlg1 : kAlg2], (unsigned long long)sh[1]);
    }
}

void LaunchDMFTop(int p, dim3 g, cudaStream_t s, const TileGeom& t, const ParClass& K,
                  const uint32_t* ranks, uint8_t* top) {
    switch (p) {
        case 1: DMFTopTileKernel<1><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 2: DMFTopTileKernel<2><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 3: DMFTopTileKernel<3><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 4: DMFTopTileKernel<4><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 5: DMFTopTileKernel<5><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 6: DMFTopTileKernel<6><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        case 7: DMFTopTileKernel<7><<<g, kBlk, 0, s>>>(t, K, ranks, top); break;
        default: break;
    }
}

template <int W>
void LaunchDMFPassW(int p, dim3 g, cudaStream_t s, const TileGeom& t, const ParClass& K,
                    const uint32_t* ranks, const uint8_t* top, uint8_t* gr, uint8_t* sg,
                    unsigned long long* cnt) {
    switch (p) {
        case 0: DMFPassTileKernel<0, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 1: DMFPassTileKernel<1, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 2: DMFPassTileKernel<2, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 3: DMFPassTileKernel<3, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 4: DMFPassTileKernel<4, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 5: DMFPassTileKernel<5, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        case 6: DMFPassTileKernel<6, W><<<g, kBlk, 0, s>>>(t, K, ranks, top, gr, sg, cnt); break;
        default: break;
    }
}

// Owned cells -> output box: unpaired (0) cells become critical; critical
// cells counted per dimension.
// Owned cells -> output box: unpaired (0) cells become critical; critical
// cells counted per dimension. Four consecutive cells along x per thread.
constexpr int kExtractX = 4;
__global__ void DMFExtractKernel(TileGeom t, const uint8_t* __restrict__ grad,
                                 const uint8_t* __restrict__ stage, uint8_t* __restrict__ outGrad,
                                 uint8_t* __restrict__ outStage, unsigned long long* counters) {
    __shared__ unsigned sh[4];
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    if (lt < 4) sh[lt] = 0;
    __syncthreads();
    const int ox0 = (int)(blockIdx.x * blockDim.x + threadIdx.x) * kExtractX;
    const int oy = (int)(blockIdx.y * blockDim.y + threadIdx.y);
    const int oz = (int)(blockIdx.z * blockDim.z + threadIdx.z);
    const int w = t.ohi[0] - t.olo[0], h = t.ohi[1] - t.olo[1], d = t.ohi[2] - t.olo[2];
    unsigned crit[4] = {0, 0, 0, 0};
    if (ox0 < w && oy < h && oz < d) {
        const int y = t.olo[1] + oy, z = t.olo[2] + oz;
        const int pyz = (y & 1) + (z & 1);
        const int c0 = t.L.Id(t.olo[0] + ox0, y, z);
        const int o0 = oz * t.opz + oy * t.opy + ox0;
        const int n = min(kExtractX, w - ox0);
#pragma unroll
        for (int k = 0; k < kExtractX; k++) {
            if (k >= n) break;
            uint8_t b = grad[c0 + k];
            if (b == 0) {
                b = dmt::MakeGradByte(dmt::kPairCritical);
                crit[pyz + ((t.olo[0] + ox0 + k) & 1)]++;
            }
            outGrad[o0 + k] = b;
            if (outStage) outStage[o0 + k] = stage[c0 + k];
        }
    }
#pragma unroll
    for (int k = 0; k < 4; k++)
        if (crit[k]) atomicAdd(&sh[k], crit[k]);
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
    constexpr int BX = 32, BY = 4;   // kBlk (32, 4, 2)
    __shared__ dmt::Star27Tables T;
    __shared__ unsigned scnt[6];             // crit0..3, written, orphans
    const int lt = threadIdx.x + BX * (threadIdx.y + BY * threadIdx.z);
    if (lt < 6) scnt[lt] = 0;
    // block origin in local vertex coords (may be -1)
    const int bx0 = t.clo[0] + (int)(blockIdx.x * BX) - 1;
    const int by0 = t.clo[1] + (int)(blockIdx.y * BY) - 1;
    const int bz0 = t.clo[2] + (int)(blockIdx.z * 2) - 1;
    LoadStarTablesShared(Tp, T);   // ends with __syncthreads
    // local vertex of the computed box
    const int vx = bx0 + 1 + (int)threadIdx.x;
    const int vy = by0 + 1 + (int)threadIdx.y;
    const int vz = bz0 + 1 + (int)threadIdx.z;
    uint32_t wmask = 0, cmask = 0;   // owned cells written / of them critical (slot masks)
    unsigned orph = 0;
    if (vx < t.chi[0] && vy < t.chi[1] && vz < t.chi[2]) {
        // the staged box holds every in-domain neighbour of a computed vertex,
        // and a neighbour outside the staged box is outside the volume
        uint32_t key[27], vvalid = 0;
#pragma unroll
        for (int dz = -1; dz <= 1; dz++)
#pragma unroll
            for (int dy = -1; dy <= 1; dy++)
#pragma unroll
                for (int dx = -1; dx <= 1; dx++) {
                    const int j = dmt::StarSlot(dx, dy, dz);
                    const int x = vx + dx, y = vy + dy, z = vz + dz;
                    const bool in = x >= 0 && y >= 0 && z >= 0 && x < t.L.n[0] && y < t.L.n[1] && z < t.L.n[2];
                    key[j] = in ? ranks[x * t.vs[0] + y * t.vs[1] + z * t.vs[2]] : 0u;
                    vvalid |= (uint32_t)in << j;
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
        if (rec.n_orphans && OwnedLocal(t, cvx, cvy, cvz)) orph = (unsigned)rec.n_orphans;
        uint32_t m = rec.done & L;   // the cells PLSStar would write
        while (m) {
            const int s = __ffs((int)m) - 1;
            m &= m - 1;
            const int x = cvx + s % 3 - 1, y = cvy + (s / 3) % 3 - 1, z = cvz + s / 9 - 1;
            if (!OwnedLocal(t, x, y, z)) continue;
            outGrad[(z - t.olo[2]) * t.opz + (y - t.olo[1]) * t.opy + (x - t.olo[0])] =
                dmt::MakeGradByte((uint8_t)rec.Code(s));
            wmask |= 1u << s;
        }
        cmask = wmask & rec.c0 & rec.c1 & rec.c2;   // pair code 7
    }
    // census fused here (every owned cell is written by exactly one star): a
    // cell's dimension is that of its star slot. Warp sums, then one shared
    // atomic per warp and counter, one global atomic per block and counter.
#pragma unroll
    for (int k = 0; k < 6; k++) {
        unsigned v = k < 4 ? (unsigned)__popc(cmask & T.dimmask[k]) : (k == 4 ? (unsigned)__popc(wmask) : orph);
        v += __shfl_xor_sync(0xffffffffu, v, 16);
        v += __shfl_xor_sync(0xffffffffu, v, 8);
        v += __shfl_xor_sync(0xffffffffu, v, 4);
        v += __shfl_xor_sync(0xffffffffu, v, 2);
        v += __shfl_xor_sync(0xffffffffu, v, 1);
        if ((lt & 31) == 0 && v) atomicAdd(&scnt[k], v);
    }
    __syncthreads();
    if (lt < 6 && scnt[lt]) {
        const int dst = lt < 4 ? kCrit0 + lt : (lt == 4 ? kWritten : kOrphans);
        atomicAdd(&counters[dst], (unsigned long long)scnt[lt]);
    }
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
    g.vs[0] = 1;
    g.vs[1] = sb.w;
    g.vs[2] = sb.w * sb.h;
    return p;
}

// Every kernel of one tile on stream ks: ranks / outG / outS are the tile's
// rank origin and owned-output origin (pitches in g). DMF uses the shared
// local-complex scratch (top, gr, sg); the lower-star rules expect outG to be
// zeroed by the caller for the census.
bool EnqueueTileWork(Rule rule, bool second_pass, int dmax, const dmt::Star27Tables& T27, const TileGeom& g,
                     const Box3D& owned, const uint32_t* ranks, uint8_t* outG, uint8_t* outS, uint8_t* top,
                     uint8_t* gr, uint8_t* sg, unsigned long long* cnt, cudaStream_t ks) {
    const int oext[3] = {owned.w, owned.h, owned.d};
    dim3 go;
    if (!GridFor(oext, go)) return false;
    if (rule == Rule::DMF) {
        const int64_t nloc = (int64_t)g.L.NumCells();
        // top[] needs no fill: entries of vertices are never read (key walks
        // stop at vertices; a selecting cofacet has dimension >= 1) and every
        // other entry is written by the top kernels before it is read
        if (cudaMemsetAsync(gr, 0, (size_t)nloc, ks) != cudaSuccess ||
            (sg && cudaMemsetAsync(sg, 0, (size_t)nloc, ks) != cudaSuccess))
            return false;
        bool ok = true;
        auto forDim = [&](int dim, auto launch) -> bool {
            for (int p = 0; p < 8; p++) {
                if (((p & 1) + ((p >> 1) & 1) + ((p >> 2) & 1)) != dim) continue;
                const ParClass K = MakeClass(g.L, p);
                if (K.ext[0] < 1 || K.ext[1] < 1 || K.ext[2] < 1) continue;
                dim3 gk;
                if (!GridFor(K.ext, gk)) return false;
                launch(p, K, gk);
                if (cudaGetLastError() != cudaSuccess) return false;
            }
            return true;
        };
        for (int dim = 1; dim <= dmax && ok; dim++)
            ok = forDim(dim, [&](int p, const ParClass& K, dim3 gk) { LaunchDMFTop(p, gk, ks, g, K, ranks, top); });
        for (int alg = 1; alg <= (second_pass ? 2 : 1) && ok; alg++)
            for (int dim = 0; dim < dmax && ok; dim++)
                ok = forDim(dim, [&](int p, const ParClass& K, dim3 gk) {
                    if (alg == 1) LaunchDMFPassW<0>(p, gk, ks, g, K, ranks, top, gr, sg, cnt);
                    else LaunchDMFPassW<1>(p, gk, ks, g, K, ranks, top, gr, sg, cnt);
                });
        if (!ok) return false;
        const int xext[3] = {(owned.w + kExtractX - 1) / kExtractX, owned.h, owned.d};
        dim3 gx;
        if (!GridFor(xext, gx)) return false;
        DMFExtractKernel<<<gx, kBlk, 0, ks>>>(g, gr, sg, outG, outS, cnt);
        return cudaGetLastError() == cudaSuccess;
    }
    const int cext[3] = {g.chi[0] - g.clo[0], g.chi[1] - g.clo[1], g.chi[2] - g.clo[2]};
    dim3 gv;
    if (!GridFor(cext, gv)) return false;
    LowerStarTileKernel<<<gv, kBlk, 0, ks>>>(T27, g, ranks, rule == Rule::PLS_BR ? 1 : 0, outG, cnt);
    return cudaGetLastError() == cudaSuccess;
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

    int pendingSlot = -1;
    size_t pendingBytes = 0;
    int64_t ownedTotal = 0;
    bool ok = true;
    int slot = 0;
    for (int64_t tile = 0; tile < grid.numTiles && ok; tile++) {
        if (pendingSlot == slot) {   // single slot: its download must be queued before reuse
            if (!ring.Download(pendingSlot, pendingBytes)) { ok = false; break; }
            pendingSlot = -1;
        }
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
        if (rule != Rule::DMF &&
            cudaMemsetAsync(outG, 0, (size_t)P.owned.Count(), ks) != cudaSuccess) { ok = false; break; }
        if (!EnqueueTileWork(rule, second_pass, dmax, T27, P.g, P.owned, ranks, outG, outS,
                             (uint8_t*)d_top.p, (uint8_t*)d_grad.p, wantStage ? (uint8_t*)d_stage.p : nullptr,
                             cnt, ks)) { ok = false; break; }
        const size_t dbytes = wantStage ? outPlane + (size_t)P.owned.Count() : (size_t)P.owned.Count();
        // download deferred until the next slot's upload is queued (copies
        // run in submission order on a single copy engine)
        if (!ring.MarkComputed(slot)) { ok = false; break; }
        if (pendingSlot >= 0 && !ring.Download(pendingSlot, pendingBytes)) { ok = false; break; }
        pendingSlot = slot;
        pendingBytes = dbytes;
        busy[(size_t)slot] = 1;
        st.launches++;
        st.bytesH2D += (int64_t)sb.Count() * 4;
        st.bytesD2H += (int64_t)P.owned.Count() * (wantStage ? 2 : 1);
        ownedTotal += P.owned.Count();
        size_t freeB = 0, totB = 0;
        if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) {
            if (st.deviceFreeMin < 0 || (int64_t)freeB < st.deviceFreeMin) st.deviceFreeMin = (int64_t)freeB;
            if ((int64_t)freeB > st.deviceFreeMax) st.deviceFreeMax = (int64_t)freeB;
        }
        slot = (slot + 1) % nslots;
    }
    if (ok && pendingSlot >= 0 && !ring.Download(pendingSlot, pendingBytes)) ok = false;
    for (int i = 0; i < nslots && ok; i++)
        if (!retire((slot + i) % nslots)) ok = false;
    if (ok) {
        DMTT_CHECK(cudaMemcpy(cntOut, d_cnt.p, sizeof(unsigned long long) * kNumCnt, cudaMemcpyDeviceToHost));
        if (rule != Rule::DMF) cntOut[kUnassigned] = (unsigned long long)ownedTotal - cntOut[kWritten];
    } else {
        cudaDeviceSynchronize();
    }
    ring.Destroy();
    st.wall_ms = (float)MsSince(t_wall0);
    if (tstats) *tstats = st;
    return ok;
}

bool PinnedRangeDmt(const void* ptr, size_t bytes) {
    if (!ptr || bytes == 0) return false;
    const unsigned char* q = (const unsigned char*)ptr;
    const void* ends[2] = {q, q + bytes - 1};
    for (const void* e : ends) {
        cudaPointerAttributes at;
        if (cudaPointerGetAttributes(&at, e) != cudaSuccess) { cudaGetLastError(); return false; }
        if (at.type != cudaMemoryTypeHost) return false;
    }
    return true;
}

// Direct streamed driver: the volume is cut into full-XY z-slabs; each slab's
// staged rank planes are one DMA from the caller's page-locked ranks into a
// slot buffer, the slab's owned cell planes one DMA from the slot back into
// the caller's gradient (and stage). Inside a slab the work still runs as
// (Tx, Ty) tiles of the same tile grid as the tiled driver, reading ranks
// through the slab's strides and writing owned cells at their slab offsets,
// so the kernels keep the small local complexes that make DMF fast (L2). Slot
// reuse is ordered on the device by events; uploads are issued `slots - 1`
// slabs ahead so that on a single copy engine the next slab comes in while
// the previous one goes out.
bool StreamedDmtDirect(Rule rule, bool second_pass, const uint32_t* ranksH, int64_t X, int64_t Y, int64_t Z,
                       int Tx, int Ty, int tz, uint8_t* gradH, uint8_t* stageH, int nslots,
                       unsigned long long cntOut[kNumCnt], TileStats3D* tstats) {
    const auto t_wall0 = std::chrono::steady_clock::now();
    const TileGrid3D grid = MakeTileGrid3D(X, Y, Z, Tx, Ty, tz);
    const int lo = (rule == Rule::DMF && second_pass) ? 2 : 1, hi = 1;
    const int64_t W = 2 * X - 1, H = 2 * Y - 1, plane = X * Y, cplane = W * H;
    const int nslab = grid.ntz;
    const int64_t perSlab = (int64_t)grid.ntx * grid.nty;
    // slab extents: staged vertex planes [zs0, zs1), owned cell planes [zo0, zo1)
    std::vector<int64_t> zs0(nslab), zs1(nslab), zo0(nslab), zo1(nslab);
    int64_t maxStagedPlanes = 0, maxOwnedPlanes = 0;
    for (int k = 0; k < nslab; k++) {
        const int64_t oz = (int64_t)k * grid.Tz, td = std::min<int64_t>(grid.Tz, Z - oz);
        zs0[k] = std::max<int64_t>(oz - lo, 0);
        zs1[k] = std::min<int64_t>(oz + td + 1 + hi, Z);
        zo0[k] = 2 * oz;
        zo1[k] = std::min<int64_t>(2 * (oz + td), 2 * Z - 1);
        maxStagedPlanes = std::max(maxStagedPlanes, zs1[k] - zs0[k]);
        maxOwnedPlanes = std::max(maxOwnedPlanes, zo1[k] - zo0[k]);
    }
    if (plane * maxStagedPlanes >= ((int64_t)1 << 31) || cplane * maxOwnedPlanes >= ((int64_t)1 << 31)) {
        std::fprintf(stderr, "dmt3d streamed: slab too large for 32-bit indexing, lower the slab depth\n");
        return false;
    }
    // largest tile-local complex (DMF scratch)
    auto cap = [&](int T, int64_t N) { return std::min<int64_t>((int64_t)T + 1 + lo + hi, N); };
    const int64_t sx = cap(grid.Tx, X), sy = cap(grid.Ty, Y), sz = cap(grid.Tz, Z);
    const int64_t localCells = (2 * sx - 1) * (2 * sy - 1) * (2 * sz - 1);
    if (localCells >= ((int64_t)1 << 31)) return false;
    const bool wantStage = rule == Rule::DMF && stageH != nullptr;
    if (nslots < 2) nslots = 2;
    if (nslots > nslab) nslots = std::max(1, nslab);

    const size_t inBytes = Align16((size_t)(plane * maxStagedPlanes) * sizeof(uint32_t));
    const size_t outPlaneB = Align16((size_t)(cplane * maxOwnedPlanes));
    const size_t outBytes = outPlaneB * (wantStage ? 2 : 1);
    std::vector<uint8_t*> d_in((size_t)nslots, nullptr), d_out((size_t)nslots, nullptr);
    DevMem d_top, d_grad, d_stage, d_cnt;
    cudaStream_t sh = nullptr, sk = nullptr, sd = nullptr;
    std::vector<cudaEvent_t> ev((size_t)nslab * 6, nullptr);
    bool ok = true;
    auto chk = [&](cudaError_t e) {
        if (e != cudaSuccess && ok) {
            std::fprintf(stderr, "dmt3d streamed: %s\n", cudaGetErrorString(e));
            ok = false;
        }
    };
    for (int i = 0; i < nslots && ok; i++) {
        chk(cudaMalloc(&d_in[(size_t)i], inBytes));
        chk(cudaMalloc(&d_out[(size_t)i], outBytes));
    }
    if (rule == Rule::DMF && ok) {
        chk(d_top.alloc((size_t)localCells));
        chk(d_grad.alloc((size_t)localCells));
        if (wantStage) chk(d_stage.alloc((size_t)localCells));
    }
    if (ok) chk(d_cnt.alloc(sizeof(unsigned long long) * kNumCnt));
    chk(cudaStreamCreateWithFlags(&sh, cudaStreamNonBlocking));
    chk(cudaStreamCreateWithFlags(&sk, cudaStreamNonBlocking));
    chk(cudaStreamCreateWithFlags(&sd, cudaStreamNonBlocking));
    for (auto& e : ev)
        if (ok) chk(cudaEventCreate(&e));
    if (ok) chk(cudaMemsetAsync(d_cnt.p, 0, sizeof(unsigned long long) * kNumCnt, sk));

    TileStats3D st;
    st.tiles = grid.numTiles;
    st.deviceBytesAllocated = (int64_t)nslots * (int64_t)(inBytes + outBytes) +
                              (rule == Rule::DMF ? localCells * (wantStage ? 3 : 2) : 0);
    dmt::Star27Tables T27;
    dmt::BuildStar27TablesCanonical(T27);
    const int dmax = MaxDimOf(grid);
    unsigned long long* cnt = (unsigned long long*)d_cnt.p;
    const int64_t ahead = nslots > 1 ? nslots - 1 : 1;   // one slot: upload k+1 after kernel k
    int64_t ownedTotal = 0;

    auto issueUpload = [&](int k) {
        const int slot = k % nslots;
        cudaEvent_t* e = &ev[(size_t)k * 6];
        if (k >= nslots) chk(cudaStreamWaitEvent(sh, ev[(size_t)(k - nslots) * 6 + 3], 0));   // slot read
        chk(cudaEventRecord(e[0], sh));
        chk(cudaMemcpyAsync(d_in[(size_t)slot], ranksH + zs0[k] * plane,
                            (size_t)((zs1[k] - zs0[k]) * plane) * sizeof(uint32_t), cudaMemcpyHostToDevice, sh));
        chk(cudaEventRecord(e[1], sh));
        st.bytesH2D += (zs1[k] - zs0[k]) * plane * 4;
    };
    auto issueWork = [&](int k) {
        const int slot = k % nslots;
        cudaEvent_t* e = &ev[(size_t)k * 6];
        if (k >= nslots) chk(cudaStreamWaitEvent(sk, ev[(size_t)(k - nslots) * 6 + 5], 0));   // slot drained
        chk(cudaStreamWaitEvent(sk, e[1], 0));
        chk(cudaEventRecord(e[2], sk));
        uint8_t* og = d_out[(size_t)slot];
        uint8_t* os = wantStage ? og + outPlaneB : nullptr;
        if (rule != Rule::DMF) chk(cudaMemsetAsync(og, 0, (size_t)((zo1[k] - zo0[k]) * cplane), sk));
        for (int64_t i = 0; i < perSlab && ok; i++) {
            const TileDesc3D d = TileAt3D(grid, (int64_t)k * perSlab + i);
            const Box3D sb = StagedBox(grid, d, lo, hi);
            TilePlan P = PlanTile(grid, d, sb);
            P.g.vs[0] = 1;
            P.g.vs[1] = (int)X;
            P.g.vs[2] = (int)plane;
            P.g.opy = (int)W;
            P.g.opz = (int)cplane;
            const uint32_t* r = (const uint32_t*)d_in[(size_t)slot] + ((sb.z0 - zs0[k]) * Y + sb.y0) * X + sb.x0;
            const int64_t ooff = ((P.owned.z0 - zo0[k]) * H + P.owned.y0) * W + P.owned.x0;
            ownedTotal += P.owned.Count();
            if (!EnqueueTileWork(rule, second_pass, dmax, T27, P.g, P.owned, r, og + ooff, os ? os + ooff : nullptr,
                                 (uint8_t*)d_top.p, (uint8_t*)d_grad.p, wantStage ? (uint8_t*)d_stage.p : nullptr,
                                 cnt, sk))
                ok = false;
            st.launches++;
        }
        chk(cudaEventRecord(e[3], sk));
    };
    auto issueDownload = [&](int k) {
        const int slot = k % nslots;
        cudaEvent_t* e = &ev[(size_t)k * 6];
        const size_t nb = (size_t)((zo1[k] - zo0[k]) * cplane);
        chk(cudaStreamWaitEvent(sd, e[3], 0));
        chk(cudaEventRecord(e[4], sd));
        chk(cudaMemcpyAsync(gradH + zo0[k] * cplane, d_out[(size_t)slot], nb, cudaMemcpyDeviceToHost, sd));
        if (wantStage)
            chk(cudaMemcpyAsync(stageH + zo0[k] * cplane, d_out[(size_t)slot] + outPlaneB, nb, cudaMemcpyDeviceToHost, sd));
        chk(cudaEventRecord(e[5], sd));
        st.bytesD2H += (int64_t)nb * (wantStage ? 2 : 1);
    };
    for (int k = 0; k < std::min<int64_t>(ahead, nslab) && ok; k++) issueUpload(k);
    for (int k = 0; k < nslab && ok; k++) {
        issueWork(k);
        if (k + ahead < nslab) issueUpload((int)(k + ahead));
        issueDownload(k);
    }
    chk(cudaDeviceSynchronize());
    if (ok) {
        for (int k = 0; k < nslab; k++) {
            const cudaEvent_t* e = &ev[(size_t)k * 6];
            float a1 = 0, a2 = 0, a3 = 0;
            cudaEventElapsedTime(&a1, e[0], e[1]);
            cudaEventElapsedTime(&a2, e[2], e[3]);
            cudaEventElapsedTime(&a3, e[4], e[5]);
            st.h2d_ms += a1; st.kernel_ms += a2; st.d2h_ms += a3;
        }
        chk(cudaMemcpy(cntOut, d_cnt.p, sizeof(unsigned long long) * kNumCnt, cudaMemcpyDeviceToHost));
        if (rule != Rule::DMF) cntOut[kUnassigned] = (unsigned long long)ownedTotal - cntOut[kWritten];
    }
    size_t freeB = 0, totB = 0;
    if (cudaMemGetInfo(&freeB, &totB) == cudaSuccess) st.deviceFreeMin = st.deviceFreeMax = (int64_t)freeB;
    for (auto& e : ev) if (e) cudaEventDestroy(e);
    if (sh) cudaStreamDestroy(sh);
    if (sk) cudaStreamDestroy(sk);
    if (sd) cudaStreamDestroy(sd);
    for (auto q : d_in) if (q) cudaFree(q);
    for (auto q : d_out) if (q) cudaFree(q);
    st.wall_ms = (float)MsSince(t_wall0);
    if (tstats) *tstats = st;
    return ok;
}

// Streamed entry: direct when every array is page-locked (or hostMode 2 pins
// them for the call), else the tiled ring driver over the host arrays.
bool StreamedDmt(Rule rule, bool second_pass, const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                 uint8_t* grad, uint8_t* stage, const DmtStreamOptions3D& o, unsigned long long c[kNumCnt],
                 TileStats3D* ts) {
    if (!ranks || !grad || X < 1 || Y < 1 || Z < 1) return false;
    if (X >= (1 << 30) || Y >= (1 << 30) || Z >= (1 << 30)) return false;
    const size_t rb = (size_t)(X * Y * Z) * sizeof(uint32_t);
    const size_t gb = (size_t)((2 * X - 1) * (2 * Y - 1) * (2 * Z - 1));
    const int Tx = o.tileX > 0 ? o.tileX : 128, Ty = o.tileY > 0 ? o.tileY : 128;
    const int tz = (int)std::min<int64_t>(o.slabDepth > 0 ? o.slabDepth : 64, Z);
    bool direct = false, registered = false;
    if (o.hostMode == 0) {
        direct = PinnedRangeDmt(ranks, rb) && PinnedRangeDmt(grad, gb) && (!stage || PinnedRangeDmt(stage, gb));
    } else if (o.hostMode == 2) {
        bool okr = cudaHostRegister((void*)ranks, rb, cudaHostRegisterDefault) == cudaSuccess;
        okr = okr && cudaHostRegister(grad, gb, cudaHostRegisterDefault) == cudaSuccess;
        okr = okr && (!stage || cudaHostRegister(stage, gb, cudaHostRegisterDefault) == cudaSuccess);
        if (okr) {
            direct = registered = true;
        } else {
            cudaGetLastError();
            cudaHostUnregister((void*)ranks);
            cudaHostUnregister(grad);
            if (stage) cudaHostUnregister(stage);
            cudaGetLastError();
        }
    }
    bool ok;
    if (direct) {
        ok = StreamedDmtDirect(rule, second_pass, ranks, X, Y, Z, Tx, Ty, tz, grad, stage, o.slots, c, ts);
    } else {
        // copy mode: the tiled ring over the host arrays
        const TileGrid3D tg = MakeTileGrid3D(X, Y, Z, Tx, Ty, tz);
        DmtTiledOptions3D to;
        to.slots = o.slots < 1 ? 3 : o.slots;
        // host copies on the library's worker pool (one z-plane of the box per task)
        const int64_t W = 2 * X - 1, H = 2 * Y - 1;
        const RankSource3D src = [&](const Box3D& b, uint32_t* dst, int64_t pitch, int64_t spitch) {
            return HostParallelFor(b.d, HostCopyThreads(b.Count() * 4), [&](int z) {
                for (int y = 0; y < b.h; y++)
                    std::memcpy(dst + z * spitch + y * pitch, ranks + ((b.z0 + z) * Y + (b.y0 + y)) * X + b.x0,
                                sizeof(uint32_t) * (size_t)b.w);
            });
        };
        auto sinkInto = [&](uint8_t* outArr) -> ByteSink3D {
            return [&, outArr](const Box3D& b, const uint8_t* sp, int64_t pitch, int64_t spitch) {
                return HostParallelFor(b.d, HostCopyThreads(b.Count()), [&](int z) {
                    for (int y = 0; y < b.h; y++)
                        std::memcpy(outArr + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0, sp + z * spitch + y * pitch,
                                    (size_t)b.w);
                });
            };
        };
        const ByteSink3D gs = sinkInto(grad);
        const ByteSink3D ss = stage ? sinkInto(stage) : ByteSink3D();
        ok = TiledDmt(rule, second_pass, src, tg, gs, stage ? &ss : nullptr, to, c, ts);
    }
    if (registered) {
        cudaHostUnregister((void*)ranks);
        cudaHostUnregister(grad);
        if (stage) cudaHostUnregister(stage);
    }
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

bool ComputeLowerStarGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                                        DmtLowerStarRule rule, uint8_t* out_grad,
                                        const DmtStreamOptions3D& options, DmtLowerStarStats* stats,
                                        TileStats3D* tstats) {
    unsigned long long c[kNumCnt] = {};
    const bool ok = StreamedDmt(rule == DmtLowerStarRule::PLS_BR ? Rule::PLS_BR : Rule::PLS, false, ranks, X, Y, Z,
                                out_grad, nullptr, options, c, tstats);
    if (ok && stats) {
        for (int d = 0; d < 4; d++) stats->critical[d] = (int64_t)c[kCrit0 + d];
        stats->orphans = (int64_t)c[kOrphans];
        stats->unassigned = (int64_t)c[kUnassigned];
    }
    return ok;
}

bool ComputeDMFGradient3DStreamed(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z, bool second_pass,
                                  uint8_t* out_grad, uint8_t* out_stage, const DmtStreamOptions3D& options,
                                  DmtDMFStats* stats, TileStats3D* tstats) {
    unsigned long long c[kNumCnt] = {};
    const bool ok = StreamedDmt(Rule::DMF, second_pass, ranks, X, Y, Z, out_grad, out_stage, options, c, tstats);
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
