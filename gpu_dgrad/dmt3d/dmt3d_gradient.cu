// 3D discrete Morse theory gradient rules on the GPU: PLS / PLS-BR
// (ComputeLowerStarGradient3DGPU) and Shivashankar-Natarajan DMF
// (ComputeDMFGradient3DGPU). Host API in dgrad_gpu_api.h.
//
// The device code calls the device-safe core of include/dmt_core unchanged
// (PLSStar / PLSBRStar / GatherStarKeys / StarBoundaryMask; for DMF the
// key-cache variants of dmt_dmf_keycache.h, which read cell keys back from the
// top-facet table); only the launch structure lives here.
// Nothing from the gi_* header stack is included.
//
//   Lower star: one thread per vertex, 3D launch (x fastest). Every cell lies
//   in the lower star of exactly one vertex, so threads write disjoint cells.
//   DMF: one top-facet kernel per dimension d = 1..dmax (each reads the keys
//   of dimension d-1 from the bytes the previous one wrote), then per algorithm one kernel per
//   dimension d = 0..dmax-1 over the cells of dimension d. A cofacet's highest
//   and second-highest facets are unique, so within one pass no two threads
//   claim the same cofacet and every cell a thread reads in the pass is
//   written either by itself or by an earlier kernel: the pass is race-free
//   and equals the CPU per-dimension pass (and hence the sequential order).
//
// Counts (critical cells, orphans, guards, claims) are block reductions in
// shared memory followed by one global atomicAdd per block.

#include "dgrad_gpu_api.h"

#include "dmt_core/dmt_hd.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_star27_tables.h"
#include "dmt_core/dmt_pls.h"
#include "dmt_core/dmt_pls_br.h"
#include "dmt_core/dmt_dmf.h"
#include "dmt_core/dmt_dmf_keycache.h"

#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include <climits>
#include <cstdio>
#include <cstring>

namespace GInt {
namespace gpu {

namespace {

#define DMT3D_CHECK(call)                                                              \
    do {                                                                               \
        cudaError_t e_ = (call);                                                       \
        if (e_ != cudaSuccess) {                                                       \
            std::fprintf(stderr, "dmt3d %s:%d %s: %s\n", __FILE__, __LINE__, #call,  \
                         cudaGetErrorString(e_));                                      \
            return false;                                                              \
        }                                                                              \
    } while (0)

// Counter slots shared by both drivers.
enum : int {
    kCntCrit0 = 0,      // 0..3 critical cells per dimension
    kCntUnassigned = 4,
    kCntOrphans = 5,
    kCntGuard1 = 6, kCntAlg1 = 7, kCntGuard2 = 8, kCntAlg2 = 9,
    kNumCounters = 10
};

struct DevBuf {
    void* p = nullptr;
    DevBuf() = default;
    DevBuf(const DevBuf&) = delete;
    DevBuf& operator=(const DevBuf&) = delete;
    ~DevBuf() { if (p) cudaFree(p); }
    cudaError_t alloc(size_t bytes) { return cudaMalloc(&p, bytes ? bytes : 1); }
    template <class T> T* as() const { return static_cast<T*>(p); }
};

struct EventSet {
    cudaEvent_t e[6] = {};
    bool ok = true;
    EventSet() { for (auto& x : e) ok = ok && cudaEventCreate(&x) == cudaSuccess; }
    ~EventSet() { for (auto& x : e) if (x) cudaEventDestroy(x); }
    float ms(int a, int b) const {
        float t = 0.f;
        cudaEventElapsedTime(&t, e[a], e[b]);
        return t;
    }
};

const dim3 kBlock3(32, 4, 2);

bool Grid3(int64_t nx, int64_t ny, int64_t nz, dim3& g) {
    const int64_t gx = (nx + kBlock3.x - 1) / kBlock3.x;
    const int64_t gy = (ny + kBlock3.y - 1) / kBlock3.y;
    const int64_t gz = (nz + kBlock3.z - 1) / kBlock3.z;
    if (gx > INT_MAX || gy > 65535 || gz > 65535) return false;
    g = dim3((unsigned)gx, (unsigned)gy, (unsigned)gz);
    return true;
}

int MaxDim(const dmt::Cubical3D& G) {
    int d = 0;
    for (int i = 0; i < 3; i++)
        if (G.n[i] > 1) d++;
    return d;
}

// Cells of one dimension, as up to three parity classes (the odd axes of the
// doubled coordinates), concatenated: thread t < start[n] maps to one cell.
struct DimClasses {
    int n;
    int pat[3];          // bit i set = coordinate i odd
    int64_t ext[3][3];   // cells per axis in the class
    int64_t start[4];    // prefix offsets; start[n] = cells of this dimension
};

DimClasses MakeDimClasses(const dmt::Cubical3D& G, int d) {
    DimClasses K;
    std::memset(&K, 0, sizeof(K));
    for (int p = 0; p < 8; p++) {
        if (((p & 1) + ((p >> 1) & 1) + ((p >> 2) & 1)) != d) continue;
        int64_t e[3], cnt = 1;
        for (int i = 0; i < 3; i++) {
            e[i] = ((p >> i) & 1) ? G.n[i] - 1 : G.n[i];
            cnt *= e[i];
        }
        if (cnt <= 0) continue;
        K.pat[K.n] = p;
        for (int i = 0; i < 3; i++) K.ext[K.n][i] = e[i];
        K.start[K.n + 1] = K.start[K.n] + cnt;
        K.n++;
    }
    return K;
}

__device__ __forceinline__ int64_t ClassCell(const dmt::Cubical3D& G, const DimClasses& K, int64_t t) {
    int k = 0;
    while (k + 1 < K.n && t >= K.start[k + 1]) k++;
    t -= K.start[k];
    const int64_t t0 = t % K.ext[k][0];
    const int64_t r = t / K.ext[k][0];
    const int64_t t1 = r % K.ext[k][1];
    const int64_t t2 = r / K.ext[k][1];
    const int p = K.pat[k];
    return G.Id(2 * t0 + (p & 1), 2 * t1 + ((p >> 1) & 1), 2 * t2 + ((p >> 2) & 1));
}

// ---------------------------------------------------------------------------
// Lower star (PLS / PLS-BR)
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

// kJournal: also write the step / part journals (StarResult path); without
// them the decisions are recorded as pair-code bit planes in registers.
template <bool kJournal>
__global__ void __launch_bounds__(256, 4) LowerStarKernel(dmt::Star27Tables Tp, dmt::Cubical3D G,
                                const uint32_t* __restrict__ ranks, int boundary_restricted,
                                uint8_t* __restrict__ grad, uint8_t* __restrict__ step,
                                uint8_t* __restrict__ part, unsigned long long* counters) {
    __shared__ dmt::Star27Tables T;
    LoadStarTablesShared(Tp, T);
    const int64_t x = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const int64_t y = (int64_t)blockIdx.y * blockDim.y + threadIdx.y;
    const int64_t z = (int64_t)blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= G.n[0] || y >= G.n[1] || z >= G.n[2]) return;
    uint32_t key[27], vvalid;
    dmt::GatherStarKeys(G, ranks, x, y, z, key, vvalid);
    if (kJournal) {
        dmt::StarResult R;
        if (boundary_restricted)
            dmt::PLSBRStar(T, key, vvalid, dmt::StarBoundaryMask(T, G, x, y, z), R);
        else
            dmt::PLSStar(T, key, vvalid, R);
        if (R.n_orphans) atomicAdd(&counters[kCntOrphans], (unsigned long long)R.n_orphans);
        for (int s = 0; s < 27; s++) {
            if (R.pair[s] == -2) continue;
            const int64_t id = G.Id(2 * x + T.off[s][0], 2 * y + T.off[s][1], 2 * z + T.off[s][2]);
            const int code = R.pair[s] < 0 ? dmt::kPairCritical : dmt::StarDir(T, s, R.pair[s]);
            grad[id] = dmt::MakeGradByte((uint8_t)code);
            if (step) step[id] = R.step[s];
            if (part) part[id] = R.part[s];
        }
    } else {
        const uint32_t bm = boundary_restricted ? dmt::StarBoundaryMask(T, G, x, y, z) : 0u;
        uint32_t L;
        dmt::StarCodeRecorder rec;
        if (boundary_restricted) dmt::PLSBRStarCodes(T, key, vvalid, bm, L, rec);
        else dmt::PLSStarCodes(T, key, vvalid, L, rec);
        if (rec.n_orphans) atomicAdd(&counters[kCntOrphans], (unsigned long long)rec.n_orphans);
        uint32_t m = rec.done & L;   // the cells PLSStar would write
        while (m) {
            const int s = __ffs((int)m) - 1;
            m &= m - 1;
            const int64_t id = G.Id(2 * x + (s % 3 - 1), 2 * y + ((s / 3) % 3 - 1), 2 * z + (s / 9 - 1));
            grad[id] = dmt::MakeGradByte((uint8_t)rec.Code(s));
        }
    }
}

// Per-cell census over the cell grid. finalize_unpaired = 0: count unassigned
// cells and critical cells per dimension (lower-star stats). = 1: every byte
// still 0 becomes critical and is counted per dimension (DMF final sweep).
__global__ void CellCensusKernel(dmt::Cubical3D G, uint8_t* __restrict__ grad, int finalize_unpaired,
                                 unsigned long long* counters) {
    __shared__ unsigned int sh[5];
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    if (lt < 5) sh[lt] = 0;
    __syncthreads();
    const int64_t cx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const int64_t cy = (int64_t)blockIdx.y * blockDim.y + threadIdx.y;
    const int64_t cz = (int64_t)blockIdx.z * blockDim.z + threadIdx.z;
    if (cx < G.c[0] && cy < G.c[1] && cz < G.c[2]) {
        const int64_t id = G.Id(cx, cy, cz);
        const int dim = (int)((cx & 1) + (cy & 1) + (cz & 1));
        const uint8_t b = grad[id];
        if (finalize_unpaired) {
            if (b == 0) {
                grad[id] = dmt::MakeGradByte(dmt::kPairCritical);
                atomicAdd(&sh[dim], 1u);
            }
        }
        else if (!dmt::GradAssigned(b)) {
            atomicAdd(&sh[4], 1u);
        }
        else if (dmt::GradCritical(b)) {
            atomicAdd(&sh[dim], 1u);
        }
    }
    __syncthreads();
    if (lt < 5 && sh[lt]) atomicAdd(&counters[lt == 4 ? kCntUnassigned : kCntCrit0 + lt], (unsigned long long)sh[lt]);
}

// ---------------------------------------------------------------------------
// DMF
// ---------------------------------------------------------------------------

// DMF kernels are templated on the parity pattern P of their cells (and on
// the selector of the pass), one 3D launch per parity class, so the key walks
// of dmt_dmf_keycache.h (DMFTopFacetsP, DMFLowestSelectingP) unroll into
// registers and a thread gets its cell's coordinates from its block indices.
// Top facets run d = 1..dmax in order (each reads the keys of dimension d-1
// back from top[]); vertices keep the 0x77 fill.
struct ClassGrid {
    int64_t ext[3];   // cells of the class per axis
    int par[3];
};

ClassGrid MakeClassGrid(const dmt::Cubical3D& G, int p) {
    ClassGrid c;
    for (int i = 0; i < 3; i++) {
        c.par[i] = (p >> i) & 1;
        c.ext[i] = c.par[i] ? G.n[i] - 1 : G.n[i];
    }
    return c;
}

__device__ __forceinline__ bool ClassCellAt3(const ClassGrid& K, int64_t& x, int64_t& y, int64_t& z) {
    const int64_t tx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    const int64_t ty = (int64_t)blockIdx.y * blockDim.y + threadIdx.y;
    const int64_t tz = (int64_t)blockIdx.z * blockDim.z + threadIdx.z;
    if (tx >= K.ext[0] || ty >= K.ext[1] || tz >= K.ext[2]) return false;
    x = 2 * tx + K.par[0];
    y = 2 * ty + K.par[1];
    z = 2 * tz + K.par[2];
    return true;
}

template <int P>
__global__ void DMFTopPKernel(dmt::Cubical3D G, ClassGrid K, const uint32_t* __restrict__ ranks, uint8_t* top) {
    int64_t x, y, z;
    if (!ClassCellAt3(K, x, y, z)) return;
    const dmt::DMFStrides<int64_t> S = dmt::MakeDMFStrides(G);
    const int64_t c = G.Id(x, y, z);
    const int64_t vb = (x >> 1) + G.n[0] * ((y >> 1) + G.n[1] * (z >> 1));
    uint8_t g0 = 7, h2 = 7;
    dmt::DMFTopFacetsP<P>(S, ranks, top, c, vb, g0, h2);
    top[c] = (uint8_t)(g0 | (h2 << 4));
}

// One algorithm (W = 0: algorithm 1, W = 1: algorithm 2) on one parity class.
// Bytes are 0 while a cell is unpaired, exactly as in the CPU driver.
template <int P, int W>
__global__ void DMFPassPKernel(dmt::Cubical3D G, ClassGrid K, const uint32_t* __restrict__ ranks,
                               const uint8_t* __restrict__ top, uint8_t* grad, uint8_t* stage,
                               unsigned long long* counters) {
    constexpr int alg = W + 1;
    __shared__ unsigned int s_guard, s_claim;
    const int lt = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
    if (lt == 0) { s_guard = 0; s_claim = 0; }
    __syncthreads();
    int64_t x, y, z;
    if (ClassCellAt3(K, x, y, z)) {
        const int64_t a = G.Id(x, y, z);
        if (!(alg == 2 && grad[a] != 0)) {
            const dmt::DMFStrides<int64_t> S = dmt::MakeDMFStrides(G);
            const int64_t vb = (x >> 1) + G.n[0] * ((y >> 1) + G.n[1] * (z >> 1));
            const int64_t q[3] = {x, y, z};
            const int dir = dmt::DMFLowestSelectingP<P, W>(G, S, ranks, top, a, vb, q);
            if (dir >= 0) {
                const int64_t st = dmt::DMFSelAxis(S.cs, dir >> 1);
                const int64_t b = (dir & 1) ? a - st : a + st;
                if (grad[a] != 0 || grad[b] != 0) {
                    atomicAdd(&s_guard, 1u);
                } else {
                    grad[a] = dmt::MakeGradByte((uint8_t)dir);
                    grad[b] = dmt::MakeGradByte((uint8_t)(dir ^ 1));
                    if (stage) { stage[a] = (uint8_t)alg; stage[b] = (uint8_t)alg; }
                    atomicAdd(&s_claim, 1u);
                }
            }
        }
    }
    __syncthreads();
    if (lt == 0) {
        if (s_guard) atomicAdd(&counters[alg == 1 ? kCntGuard1 : kCntGuard2], (unsigned long long)s_guard);
        if (s_claim) atomicAdd(&counters[alg == 1 ? kCntAlg1 : kCntAlg2], (unsigned long long)s_claim);
    }
}

bool LaunchDMFTopP(int p, dim3 g, const dmt::Cubical3D& G, const ClassGrid& K, const uint32_t* r, uint8_t* top) {
    switch (p) {
        case 1: DMFTopPKernel<1><<<g, kBlock3>>>(G, K, r, top); break;
        case 2: DMFTopPKernel<2><<<g, kBlock3>>>(G, K, r, top); break;
        case 3: DMFTopPKernel<3><<<g, kBlock3>>>(G, K, r, top); break;
        case 4: DMFTopPKernel<4><<<g, kBlock3>>>(G, K, r, top); break;
        case 5: DMFTopPKernel<5><<<g, kBlock3>>>(G, K, r, top); break;
        case 6: DMFTopPKernel<6><<<g, kBlock3>>>(G, K, r, top); break;
        case 7: DMFTopPKernel<7><<<g, kBlock3>>>(G, K, r, top); break;
        default: return false;
    }
    return cudaGetLastError() == cudaSuccess;
}

template <int W>
bool LaunchDMFPassP(int p, dim3 g, const dmt::Cubical3D& G, const ClassGrid& K, const uint32_t* r,
                    const uint8_t* top, uint8_t* grad, uint8_t* stage, unsigned long long* cnt) {
    switch (p) {
        case 0: DMFPassPKernel<0, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 1: DMFPassPKernel<1, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 2: DMFPassPKernel<2, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 3: DMFPassPKernel<3, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 4: DMFPassPKernel<4, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 5: DMFPassPKernel<5, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        case 6: DMFPassPKernel<6, W><<<g, kBlock3>>>(G, K, r, top, grad, stage, cnt); break;
        default: return false;
    }
    return cudaGetLastError() == cudaSuccess;
}

// (dimension, key) order: per dimension, LSD radix sort of the cell ids on
// key components 2^d - 1 down to 0 (stable passes, so the result is the
// lexicographic order of the flattened keys).
__global__ void ClassIdsKernel(dmt::Cubical3D G, DimClasses K, int64_t* ids) {
    const int64_t t = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (t < K.start[K.n]) ids[t] = ClassCell(G, K, t);
}

__global__ void KeyComponentKernel(dmt::Cubical3D G, const uint32_t* __restrict__ ranks,
                                   const uint8_t* __restrict__ top, const int64_t* __restrict__ ids,
                                   int64_t n, int comp, uint32_t* __restrict__ out) {
    const int64_t i = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = dmt::DMFKeyComponentFromTop(G, ranks, top, ids[i], comp);
}

__global__ void ScatterOrderKernel(const int64_t* __restrict__ ids, int64_t n, int64_t base,
                                   int64_t* __restrict__ order) {
    const int64_t i = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) order[ids[i]] = base + i;
}

unsigned Blocks1(int64_t n, int block) { return (unsigned)((n + block - 1) / block); }

// Significant bits of the largest rank: the radix sorts skip the rest. Ranks
// need not be a permutation of 0..nv-1, only pairwise distinct.
int RankBits(const uint32_t* ranks, int64_t nv) {
    uint32_t m = 0;
    for (int64_t i = 0; i < nv; i++) m = ranks[i] > m ? ranks[i] : m;
    int b = 1;
    while (b < 32 && (m >> b) != 0) b++;
    return b;
}

bool DeviceOk() {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

} // namespace

bool LowerStarGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                            bool boundary_restricted, uint8_t* out_grad,
                            uint8_t* out_step, uint8_t* out_part,
                            DmtLowerStarStats* stats, Dmt3DTimings* timings) {
    if (!ranks || !out_grad || X < 1 || Y < 1 || Z < 1 || !DeviceOk()) return false;
    const dmt::Cubical3D G = dmt::Cubical3D::Make(X, Y, Z);
    dmt::Star27Tables T;
    dmt::BuildStar27TablesCanonical(T);
    const int64_t nv = G.NumVertices(), nc = G.NumCells();
    dim3 gv, gc;
    if (!Grid3(X, Y, Z, gv) || !Grid3(G.c[0], G.c[1], G.c[2], gc)) return false;

    EventSet ev;
    if (!ev.ok) return false;
    DevBuf d_ranks, d_grad, d_step, d_part, d_cnt;
    DMT3D_CHECK(d_ranks.alloc(sizeof(uint32_t) * (size_t)nv));
    DMT3D_CHECK(d_grad.alloc((size_t)nc));
    if (out_step) DMT3D_CHECK(d_step.alloc((size_t)nc));
    if (out_part) DMT3D_CHECK(d_part.alloc((size_t)nc));
    DMT3D_CHECK(d_cnt.alloc(sizeof(unsigned long long) * kNumCounters));

    DMT3D_CHECK(cudaEventRecord(ev.e[0]));
    DMT3D_CHECK(cudaMemcpy(d_ranks.p, ranks, sizeof(uint32_t) * (size_t)nv, cudaMemcpyHostToDevice));
    DMT3D_CHECK(cudaEventRecord(ev.e[1]));
    DMT3D_CHECK(cudaMemset(d_grad.p, 0, (size_t)nc));
    if (out_step) DMT3D_CHECK(cudaMemset(d_step.p, 0, (size_t)nc));
    if (out_part) DMT3D_CHECK(cudaMemset(d_part.p, 0, (size_t)nc));
    DMT3D_CHECK(cudaMemset(d_cnt.p, 0, sizeof(unsigned long long) * kNumCounters));
    if (out_step || out_part)
        LowerStarKernel<true><<<gv, kBlock3>>>(T, G, d_ranks.as<uint32_t>(), boundary_restricted ? 1 : 0,
                                               d_grad.as<uint8_t>(), d_step.as<uint8_t>(), d_part.as<uint8_t>(),
                                               d_cnt.as<unsigned long long>());
    else
        LowerStarKernel<false><<<gv, kBlock3>>>(T, G, d_ranks.as<uint32_t>(), boundary_restricted ? 1 : 0,
                                                d_grad.as<uint8_t>(), nullptr, nullptr,
                                                d_cnt.as<unsigned long long>());
    DMT3D_CHECK(cudaGetLastError());
    CellCensusKernel<<<gc, kBlock3>>>(G, d_grad.as<uint8_t>(), 0, d_cnt.as<unsigned long long>());
    DMT3D_CHECK(cudaGetLastError());
    DMT3D_CHECK(cudaEventRecord(ev.e[2]));
    DMT3D_CHECK(cudaMemcpy(out_grad, d_grad.p, (size_t)nc, cudaMemcpyDeviceToHost));
    if (out_step) DMT3D_CHECK(cudaMemcpy(out_step, d_step.p, (size_t)nc, cudaMemcpyDeviceToHost));
    if (out_part) DMT3D_CHECK(cudaMemcpy(out_part, d_part.p, (size_t)nc, cudaMemcpyDeviceToHost));
    unsigned long long cnt[kNumCounters];
    DMT3D_CHECK(cudaMemcpy(cnt, d_cnt.p, sizeof(cnt), cudaMemcpyDeviceToHost));
    DMT3D_CHECK(cudaEventRecord(ev.e[3]));
    DMT3D_CHECK(cudaEventSynchronize(ev.e[3]));

    if (stats) {
        for (int d = 0; d < 4; d++) stats->critical[d] = (int64_t)cnt[kCntCrit0 + d];
        stats->orphans = (int64_t)cnt[kCntOrphans];
        stats->unassigned = (int64_t)cnt[kCntUnassigned];
    }
    if (timings) {
        timings->h2d_ms = ev.ms(0, 1);
        timings->kernel_ms = ev.ms(1, 2);
        timings->d2h_ms = ev.ms(2, 3);
    }
    return true;
}

bool ComputeDMFGradient3DGPU(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
                             bool second_pass, uint8_t* out_grad, uint8_t* out_stage,
                             int64_t* out_order, DmtDMFStats* stats, Dmt3DTimings* timings) {
    if (!ranks || !out_grad || X < 1 || Y < 1 || Z < 1 || !DeviceOk()) return false;
    const dmt::Cubical3D G = dmt::Cubical3D::Make(X, Y, Z);
    const int64_t nv = G.NumVertices(), nc = G.NumCells();
    const int dmax = MaxDim(G);
    dim3 gc;
    if (!Grid3(G.c[0], G.c[1], G.c[2], gc)) return false;
    DimClasses K[4];
    for (int d = 0; d < 4; d++) K[d] = MakeDimClasses(G, d);

    EventSet ev;
    if (!ev.ok) return false;
    DevBuf d_ranks, d_grad, d_stage, d_top, d_cnt, d_order;
    DMT3D_CHECK(d_ranks.alloc(sizeof(uint32_t) * (size_t)nv));
    DMT3D_CHECK(d_grad.alloc((size_t)nc));
    DMT3D_CHECK(d_top.alloc((size_t)nc));
    if (out_stage) DMT3D_CHECK(d_stage.alloc((size_t)nc));
    if (out_order) DMT3D_CHECK(d_order.alloc(sizeof(int64_t) * (size_t)nc));
    DMT3D_CHECK(d_cnt.alloc(sizeof(unsigned long long) * kNumCounters));

    DMT3D_CHECK(cudaEventRecord(ev.e[0]));
    DMT3D_CHECK(cudaMemcpy(d_ranks.p, ranks, sizeof(uint32_t) * (size_t)nv, cudaMemcpyHostToDevice));
    DMT3D_CHECK(cudaEventRecord(ev.e[1]));
    DMT3D_CHECK(cudaMemset(d_grad.p, 0, (size_t)nc));
    if (out_stage) DMT3D_CHECK(cudaMemset(d_stage.p, 0, (size_t)nc));
    DMT3D_CHECK(cudaMemset(d_cnt.p, 0, sizeof(unsigned long long) * kNumCounters));
    unsigned long long* cntp = d_cnt.as<unsigned long long>();

    const int kBlock1 = 256;
    DMT3D_CHECK(cudaMemset(d_top.p, 0x77, (size_t)nc));
    // one launch per parity class of the dimension
    auto forDim = [&](int d, auto launch) -> bool {
        for (int p = 0; p < 8; p++) {
            if (((p & 1) + ((p >> 1) & 1) + ((p >> 2) & 1)) != d) continue;
            const ClassGrid C = MakeClassGrid(G, p);
            if (C.ext[0] < 1 || C.ext[1] < 1 || C.ext[2] < 1) continue;
            dim3 g;
            if (!Grid3(C.ext[0], C.ext[1], C.ext[2], g)) return false;
            if (!launch(p, g, C)) return false;
        }
        return true;
    };
    for (int d = 1; d <= dmax; d++)
        if (!forDim(d, [&](int p, dim3 g, const ClassGrid& C) {
                return LaunchDMFTopP(p, g, G, C, d_ranks.as<uint32_t>(), d_top.as<uint8_t>());
            }))
            return false;
    for (int alg = 1; alg <= (second_pass ? 2 : 1); alg++)
        for (int d = 0; d < dmax; d++)
            if (!forDim(d, [&](int p, dim3 g, const ClassGrid& C) {
                    return alg == 1
                        ? LaunchDMFPassP<0>(p, g, G, C, d_ranks.as<uint32_t>(), d_top.as<uint8_t>(),
                                            d_grad.as<uint8_t>(), d_stage.as<uint8_t>(), cntp)
                        : LaunchDMFPassP<1>(p, g, G, C, d_ranks.as<uint32_t>(), d_top.as<uint8_t>(),
                                            d_grad.as<uint8_t>(), d_stage.as<uint8_t>(), cntp);
                }))
                return false;
    CellCensusKernel<<<gc, kBlock3>>>(G, d_grad.as<uint8_t>(), 1, cntp);
    DMT3D_CHECK(cudaGetLastError());

    if (out_order) {
        int64_t maxn = 0;
        for (int d = 0; d < 4; d++) maxn = K[d].start[K[d].n] > maxn ? K[d].start[K[d].n] : maxn;
        if (maxn > INT_MAX) return false;
        DevBuf ids0, ids1, key0, key1, temp;
        DMT3D_CHECK(ids0.alloc(sizeof(int64_t) * (size_t)maxn));
        DMT3D_CHECK(ids1.alloc(sizeof(int64_t) * (size_t)maxn));
        DMT3D_CHECK(key0.alloc(sizeof(uint32_t) * (size_t)maxn));
        DMT3D_CHECK(key1.alloc(sizeof(uint32_t) * (size_t)maxn));
        const int bits = RankBits(ranks, nv);
        size_t temp_bytes = 0;
        {
            cub::DoubleBuffer<uint32_t> kb(key0.as<uint32_t>(), key1.as<uint32_t>());
            cub::DoubleBuffer<int64_t> vb(ids0.as<int64_t>(), ids1.as<int64_t>());
            DMT3D_CHECK(cub::DeviceRadixSort::SortPairs(nullptr, temp_bytes, kb, vb, (int)maxn, 0, bits));
        }
        DMT3D_CHECK(temp.alloc(temp_bytes));
        int64_t base = 0;
        for (int d = 0; d < 4; d++) {
            const int64_t n = K[d].start[K[d].n];
            if (n == 0) continue;
            cub::DoubleBuffer<uint32_t> kb(key0.as<uint32_t>(), key1.as<uint32_t>());
            cub::DoubleBuffer<int64_t> vb(ids0.as<int64_t>(), ids1.as<int64_t>());
            ClassIdsKernel<<<Blocks1(n, kBlock1), kBlock1>>>(G, K[d], vb.Current());
            DMT3D_CHECK(cudaGetLastError());
            for (int comp = (1 << d) - 1; comp >= 0; comp--) {
                KeyComponentKernel<<<Blocks1(n, kBlock1), kBlock1>>>(G, d_ranks.as<uint32_t>(), d_top.as<uint8_t>(),
                                                                    vb.Current(), n, comp, kb.Current());
                DMT3D_CHECK(cudaGetLastError());
                size_t tb = temp_bytes;
                DMT3D_CHECK(cub::DeviceRadixSort::SortPairs(temp.p, tb, kb, vb, (int)n, 0, bits));
            }
            ScatterOrderKernel<<<Blocks1(n, kBlock1), kBlock1>>>(vb.Current(), n, base, d_order.as<int64_t>());
            DMT3D_CHECK(cudaGetLastError());
            base += n;
        }
    }
    DMT3D_CHECK(cudaEventRecord(ev.e[2]));
    DMT3D_CHECK(cudaMemcpy(out_grad, d_grad.p, (size_t)nc, cudaMemcpyDeviceToHost));
    if (out_stage) DMT3D_CHECK(cudaMemcpy(out_stage, d_stage.p, (size_t)nc, cudaMemcpyDeviceToHost));
    if (out_order) DMT3D_CHECK(cudaMemcpy(out_order, d_order.p, sizeof(int64_t) * (size_t)nc, cudaMemcpyDeviceToHost));
    unsigned long long cnt[kNumCounters];
    DMT3D_CHECK(cudaMemcpy(cnt, d_cnt.p, sizeof(cnt), cudaMemcpyDeviceToHost));
    DMT3D_CHECK(cudaEventRecord(ev.e[3]));
    DMT3D_CHECK(cudaEventSynchronize(ev.e[3]));

    if (stats) {
        stats->guard1 = (int64_t)cnt[kCntGuard1];
        stats->guard2 = (int64_t)cnt[kCntGuard2];
        stats->alg1 = (int64_t)cnt[kCntAlg1];
        stats->alg2 = (int64_t)cnt[kCntAlg2];
        for (int d = 0; d < 4; d++) stats->critical[d] = (int64_t)cnt[kCntCrit0 + d];
    }
    if (timings) {
        timings->h2d_ms = ev.ms(0, 1);
        timings->kernel_ms = ev.ms(1, 2);
        timings->d2h_ms = ev.ms(2, 3);
    }
    return true;
}

} // namespace gpu
} // namespace GInt
