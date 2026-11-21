// 3D lower-star (Robins) discrete gradient, fused and tiled.
//
// One thread per vertex of a view's computed box, from a 3D launch (no id to
// coordinate divisions). Each thread reproduces, for its vertex, the CPU's
// lower-star homotopy expansion
//   CLSSlidingGradient::ComputeLowerStar / HomotopyExpand (interior
//   vertices) and CLSGradient::ComputeLowerStar / HomotopyExpand
//   (boundary vertices; constrained pass),
// re-indexed onto the 27 cells around the vertex. No max/min label arrays:
// lower-star membership and lowest-vertex queries are evaluated on the fly
// from a shared-memory block of function values.
//
// Bit-exactness rests on replicating the CPU's orders and comparators:
//   - total order IsGreater(a,b): value first, then vertex index. Only
//     vertices of one 3x3x3 window are ever compared, and the window position
//     w = (ex+1)+3(ey+1)+9(ez+1) orders them exactly as the global index does,
//     so the tie-break needs no global (int64) index. The CPU max/min
//     labelings are argmax/argmin under this order (higher index wins a max
//     tie, lower index a min tie), which is what the fused tests compute.
//   - lower-star insertion order = AdjacentCellsIterator = ascending slot;
//     this fixes the d-cell scan order and the order of PickLowestCandidate
//     among unpaired d-cells (earliest wins a tie).
//   - candidate cofacets in CofacetsIterator order (table), earliest wins a
//     tie of lowest vertices.
//   - the vertex is paired with the edge whose other endpoint is lowest.
//   - each lower star is split into subsets keyed by boundaryValue(cell)
//     (+ class*maxDim() when constrained) and every subset is expanded on its
//     own; the subsets are disjoint and never interact, so the expansion order
//     of the subsets cannot change the result.
// The algorithm is comparisons-only, so matching orders and comparators makes
// the output memcmp-equal to the CPU, not just close.
//
// Expansion state lives in 27-bit masks: lower star, subset, paired, the
// "missing facets" counter as two bit planes (a lower-star cell has at most
// dim <= 3 facets inside the window - exactly the ones containing the
// vertex), and the pair code as three bit planes. No per-slot arrays, no
// local memory.
//
// This file must not include any gi_* headers; the tables arrive generated
// from the real meshes (cls3d_tables.h).

#include "dgrad_gpu_api.h"
#include "cls3d_internal.cuh"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>

namespace {

struct Tab3 {
    unsigned fac[27];    // facets inside the window (mask)
    unsigned cof[27];    // cofacets inside the window (mask)
    unsigned vert[27];   // vertices of the cell (mask of window positions)
    unsigned cofl[27];   // cofacets in CofacetsIterator order, 4 x 8 bits, 0xFF = none
    unsigned dimm[4];    // slots of dimension d
    unsigned axneg[3];   // slots with d_axis = -1
    unsigned axzero[3];  // slots with d_axis = 0
    unsigned axpos[3];   // slots with d_axis = +1
};

__constant__ Tab3 c_t3;

constexpr unsigned kAll27 = (1u << 27) - 1u;
constexpr int kCenter = 13;

// IsGreater for two window vertices: value first, then window position.
__device__ __forceinline__ bool greater_w(float va, int wa, float vb, int wb) {
    if (va > vb) return true;
    if (vb > va) return false;
    return wa > wb;
}

template <int BX, int BY, int BZ>
struct Win {
    const float* sh;      // shared value block
    const int* soff;      // window position -> shared offset
    int sb;               // this thread's vertex in sh

    __device__ __forceinline__ float val(int w) const { return sh[sb + soff[w]]; }

    // lowest vertex (window position) of the cell in slot s: argmin under
    // IsGreater. Positions are visited ascending, so a value tie keeps the
    // earlier (lower-index) vertex - the min labeling's tie rule.
    __device__ __forceinline__ int lowest(const unsigned* svert, int s, float& lv) const {
        unsigned m = svert[s];
        int best = __ffs(m) - 1;
        m &= m - 1;
        float bv = val(best);
        while (m) {
            const int w = __ffs(m) - 1;
            m &= m - 1;
            const float v = val(w);
            if (bv > v) { bv = v; best = w; }
        }
        lv = bv;
        return best;
    }
};

struct Planes {
    unsigned b0 = 0, b1 = 0;            // missing-facet counter bit planes
    unsigned c0 = 0, c1 = 0, c2 = 0;    // raw pair code bit planes
    unsigned paired = 0;
};

__device__ __forceinline__ void set_code(Planes& P, int s, unsigned code) {
    P.c0 |= (code & 1u) << s;
    P.c1 |= ((code >> 1) & 1u) << s;
    P.c2 |= ((code >> 2) & 1u) << s;
}

// counter-- on every cofacet of s inside the subset
__device__ __forceinline__ void dec_cofacets(Planes& P, const unsigned* scof, int s,
                                             unsigned sub) {
    const unsigned M = scof[s] & sub;
    const unsigned borrow = M & ~P.b0;
    P.b0 ^= M;
    P.b1 ^= borrow;
}

__device__ __forceinline__ void make_pair(Planes& P, const unsigned* scof, int s, int p,
                                          unsigned sub) {
    const int diff = p - s;
    const int ad = diff < 0 ? -diff : diff;
    const unsigned axis = ad == 1 ? 0u : (ad == 3 ? 1u : 2u);
    set_code(P, s, 2u * axis + (diff < 0 ? 1u : 0u));
    set_code(P, p, 2u * axis + (diff > 0 ? 1u : 0u));
    P.paired |= (1u << s) | (1u << p);
    dec_cofacets(P, scof, s, sub);
    dec_cofacets(P, scof, p, sub);
}

__device__ __forceinline__ void make_critical(Planes& P, const unsigned* scof, int s,
                                              unsigned sub) {
    set_code(P, s, 7u);
    P.paired |= 1u << s;
    dec_cofacets(P, scof, s, sub);
}

// HomotopyExpand on one subset 'sub' of the lower star.
template <class W>
__device__ void expand_subset(const W& win, const unsigned* sfac, const unsigned* scof,
                              const unsigned* svert, const unsigned* scofl, unsigned sub,
                              Planes& P) {
    // missing-facet counters: facets of s inside the subset (0..3)
    {
        unsigned m = sub;
        while (m) {
            const int s = __ffs(m) - 1;
            m &= m - 1;
            const unsigned n = (unsigned)__popc(sfac[s] & sub);
            P.b0 |= (n & 1u) << s;
            P.b1 |= ((n >> 1) & 1u) << s;
        }
    }

    // vertex phase: pair the vertex with the edge whose other endpoint is
    // lowest (that endpoint's window position is the edge's slot)
    if (sub >> kCenter & 1u) {
        unsigned E = sub & c_t3.dimm[1];
        if (E == 0) {
            make_critical(P, scof, kCenter, sub);
        } else {
            int best = __ffs(E) - 1;
            E &= E - 1;
            float bv = win.val(best);
            while (E) {
                const int s = __ffs(E) - 1;
                E &= E - 1;
                const float v = win.val(s);
                // CPU: switch iff Before(other, current) == IsGreater(current, other)
                if (greater_w(bv, best, v, s)) { bv = v; best = s; }
            }
            make_pair(P, scof, kCenter, best, sub);
        }
    }

    // dimensions 1 and 2: restart-scan pairing, else lowest unpaired critical
    for (int d = 1; d <= 2; d++) {
        unsigned U = sub & c_t3.dimm[d] & ~P.paired;
        while (U) {
            bool advanced = false;
            unsigned scan = U;
            while (scan) {
                const int s = __ffs(scan) - 1;
                scan &= scan - 1;
                const unsigned one = P.b0 & ~P.b1;   // counter == 1
                const unsigned cand = scof[s] & sub & one;
                if (cand == 0) continue;
                int p;
                if ((cand & (cand - 1)) == 0) {
                    p = __ffs(cand) - 1;
                } else {
                    // PickLowestCandidate over cofacets in iterator order
                    const unsigned L = scofl[s];
                    p = -1;
                    int pw = 0;
                    float pv = 0.f;
                    for (int k = 0; k < 4; k++) {
                        const int q = (int)((L >> (8 * k)) & 0xFFu);
                        if (q == 0xFF) break;
                        if (!(cand >> q & 1u)) continue;
                        float qv;
                        const int qw = win.lowest(svert, q, qv);
                        if (p < 0 || greater_w(pv, pw, qv, qw)) { p = q; pw = qw; pv = qv; }
                    }
                }
                make_pair(P, scof, s, p, sub);
                U &= ~(1u << s);
                advanced = true;
                break;  // the CPU restarts the d-cell scan after every pair
            }
            if (!advanced) {
                unsigned m = U;
                int best = __ffs(m) - 1;
                m &= m - 1;
                float bv;
                int bw = win.lowest(svert, best, bv);
                while (m) {
                    const int s = __ffs(m) - 1;
                    m &= m - 1;
                    float v;
                    const int w = win.lowest(svert, s, v);
                    if (greater_w(bv, bw, v, w)) { best = s; bw = w; bv = v; }
                }
                make_critical(P, scof, best, sub);
                U &= ~(1u << best);
            }
        }
    }
    // dimension 3: no cofacets, so every unpaired 3-cell becomes critical
    // (the CPU's one-at-a-time order only affects counters of cells that
    // do not exist)
    {
        const unsigned U = sub & c_t3.dimm[3] & ~P.paired;
        P.c0 |= U;
        P.c1 |= U;
        P.c2 |= U;
        P.paired |= U;
    }
}

// kInterior = true: the interior kernel. Its views hold only vertices strictly
// inside the domain (1 <= g <= N-2 on every axis) and no restriction classes,
// so every cell of every star exists and has boundaryValue 0: no validity
// mask, no boundary classes, no class bytes, one expansion per star.
// kInterior = false: the boundary / label kernel, everything else (the
// domain-face shell, and every vertex of a constrained run).
template <int BX, int BY, int BZ, bool kInterior>
__global__ void __launch_bounds__(BX * BY * BZ)
dgrad3d_kernel(const GInt::gpu::DeviceTileView3D* __restrict__ views, int nbz) {
    constexpr int SX = BX + 2, SY = BY + 2, SZ = BZ + 2;
    constexpr int NT = BX * BY * BZ;
    __shared__ float sh[SX * SY * SZ];
    __shared__ unsigned s_fac[27], s_cof[27], s_vert[27], s_cofl[27];
    __shared__ int s_off[27];

    const int vi = (int)blockIdx.z / nbz;
    const int bz = (int)blockIdx.z - vi * nbz;
    const GInt::gpu::DeviceTileView3D& v = views[vi];
    // blocks beyond this view's computed box (the grid covers the largest view
    // of the launch) leave before loading anything; uniform per block
    if ((int)blockIdx.x * BX >= v.cw || (int)blockIdx.y * BY >= v.ch || bz * BZ >= v.cd) return;

    const int tid = (int)threadIdx.x + BX * ((int)threadIdx.y + BY * (int)threadIdx.z);
    if (tid < 27) {
        s_fac[tid] = c_t3.fac[tid];
        s_cof[tid] = c_t3.cof[tid];
        s_vert[tid] = c_t3.vert[tid];
        s_cofl[tid] = c_t3.cofl[tid];
        s_off[tid] = (tid % 3 - 1) + ((tid / 3) % 3 - 1) * SX + (tid / 9 - 1) * SX * SY;
    }

    // block origin (global vertex coords of sh[0]); may be outside the domain
    const int bx0 = v.cx0 + (int)blockIdx.x * BX - 1;
    const int by0 = v.cy0 + (int)blockIdx.y * BY - 1;
    const int bz0 = v.cz0 + bz * BZ - 1;
    // cooperative load, clamped to the staged box (which holds every value a
    // computed vertex can need; clamped entries are only read for cells
    // outside the domain, which the validity mask removes)
    const int sx1 = v.sx0 + v.svx - 1, sy1 = v.sy0 + v.svy - 1, sz1 = v.sz0 + v.svz - 1;
    for (int i = tid; i < SX * SY * SZ; i += NT) {
        const int lx = i % SX, ly = (i / SX) % SY, lz = i / (SX * SY);
        const int gx = min(max(bx0 + lx, v.sx0), sx1);
        const int gy = min(max(by0 + ly, v.sy0), sy1);
        const int gz = min(max(bz0 + lz, v.sz0), sz1);
        sh[i] = v.vals[(gz - v.sz0) * v.vpz + (gy - v.sy0) * v.vpy + (gx - v.sx0)];
    }
    __syncthreads();

    const int lcx = (int)blockIdx.x * BX + (int)threadIdx.x;
    const int lcy = (int)blockIdx.y * BY + (int)threadIdx.y;
    const int lcz = bz * BZ + (int)threadIdx.z;
    if (lcx >= v.cw || lcy >= v.ch || lcz >= v.cd) return;
    const int gx = v.cx0 + lcx, gy = v.cy0 + lcy, gz = v.cz0 + lcz;

    Win<BX, BY, BZ> win;
    win.sh = sh;
    win.soff = s_off;
    win.sb = ((int)threadIdx.z + 1) * SX * SY + ((int)threadIdx.y + 1) * SX + ((int)threadIdx.x + 1);

    // G: window vertices strictly below the center under IsGreater
    const float cv = sh[win.sb];
    unsigned G = 0;
#pragma unroll
    for (int w = 0; w < 27; w++) {
        if (w == kCenter) continue;
        const float nv = sh[win.sb + (w % 3 - 1) + ((w / 3) % 3 - 1) * SX + (w / 9 - 1) * SX * SY];
        const bool gt = w < kCenter ? !(nv > cv) : (cv > nv);
        G |= (unsigned)gt << w;
    }
    // cells inside the domain (all of them for an interior vertex)
    unsigned invalid = 0;
    if (!kInterior) {
        if (gx == 0) invalid |= c_t3.axneg[0];
        if (gx == v.X - 1) invalid |= c_t3.axpos[0];
        if (gy == 0) invalid |= c_t3.axneg[1];
        if (gy == v.Y - 1) invalid |= c_t3.axpos[1];
        if (gz == 0) invalid |= c_t3.axneg[2];
        if (gz == v.Z - 1) invalid |= c_t3.axpos[2];
    }
    // lower star: every vertex of the cell other than the center is below it
    unsigned L = 1u << kCenter;
#pragma unroll
    for (int s = 0; s < 27; s++) {
        if (s == kCenter) continue;
        if ((c_t3.vert[s] & ~(1u << kCenter) & ~G) == 0) L |= 1u << s;
    }
    L &= ~invalid & kAll27;

    // cells this view owns: owner vertex floor(c/2) = g + min(d, 0) per axis
    auto own_axis = [&](int g, int o, int t, int a) -> unsigned {
        unsigned m = 0;
        if (g >= o && g < o + t) m |= c_t3.axzero[a] | c_t3.axpos[a];
        if (g - 1 >= o && g - 1 < o + t) m |= c_t3.axneg[a];
        return m;
    };
    const unsigned own = own_axis(gx, v.ox, v.tw, 0) & own_axis(gy, v.oy, v.th, 1) &
                         own_axis(gz, v.oz, v.td, 2);
    const unsigned W = L & own;
    if (W == 0) return;   // nothing of this lower star is written by this view

    Planes P;
    if (kInterior) {
        // one subset: no boundary cells, no classes
        expand_subset(win, s_fac, s_cof, s_vert, s_cofl, L, P);
    } else {
    // boundary classes: only even cell coordinates (d = 0) can be 0 or 2N-2
    const unsigned bxm = (gx == 0 || gx == v.X - 1) ? c_t3.axzero[0] : 0u;
    const unsigned bym = (gy == 0 || gy == v.Y - 1) ? c_t3.axzero[1] : 0u;
    const unsigned bzm = (gz == 0 || gz == v.Z - 1) ? c_t3.axzero[2] : 0u;

    if (v.cls == nullptr) {
        if ((bxm | bym | bzm) == 0) {
            expand_subset(win, s_fac, s_cof, s_vert, s_cofl, L, P);
        } else {
            const unsigned b3 = bxm & bym & bzm;
            const unsigned b23 = (bxm & bym) | (bxm & bzm) | (bym & bzm);
            const unsigned odd = bxm ^ bym ^ bzm;
            const unsigned subs[4] = {L & ~(bxm | bym | bzm), L & odd & ~b3, L & b23 & ~b3,
                                      L & b3};
            for (int k = 0; k < 4; k++)
                if (subs[k]) expand_subset(win, s_fac, s_cof, s_vert, s_cofl, subs[k], P);
        }
    } else {
        // key = class*maxDim + boundaryValue, clamped; stored as 8 bit planes
        unsigned kp[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        unsigned m = L;
        while (m) {
            const int s = __ffs(m) - 1;
            m &= m - 1;
            const int dx = s % 3 - 1, dy = (s / 3) % 3 - 1, dz = s / 9 - 1;
            const int cxl = 2 * gx + dx - v.kcx0, cyl = 2 * gy + dy - v.kcy0,
                      czl = 2 * gz + dz - v.kcz0;
            const int cls = (int)v.cls[czl * v.kpz + cyl * v.kpy + cxl];
            const int bvv = (int)(bxm >> s & 1u) + (int)(bym >> s & 1u) + (int)(bzm >> s & 1u);
            int key = cls * 3 + bvv;
            if (key > v.keyClamp) key = v.keyClamp;
#pragma unroll
            for (int b = 0; b < 8; b++) kp[b] |= (unsigned)((key >> b) & 1) << s;
        }
        unsigned rem = L;
        while (rem) {
            const int s0 = __ffs(rem) - 1;
            unsigned sub = rem;
#pragma unroll
            for (int b = 0; b < 8; b++) sub &= (kp[b] >> s0 & 1u) ? kp[b] : ~kp[b];
            expand_subset(win, s_fac, s_cof, s_vert, s_cofl, sub, P);
            rem &= ~sub;
        }
    }
    }   // !kInterior

    // scatter the owned bytes: assigned = 1, pair = code
    unsigned m = W;
    while (m) {
        const int s = __ffs(m) - 1;
        m &= m - 1;
        const unsigned raw = (P.c0 >> s & 1u) | ((P.c1 >> s & 1u) << 1) | ((P.c2 >> s & 1u) << 2);
        const unsigned code = v.dircode[raw];
        const int dx = s % 3 - 1, dy = (s / 3) % 3 - 1, dz = s / 9 - 1;
        const int lx = 2 * gx + dx - v.gcx0, ly = 2 * gy + dy - v.gcy0, lz = 2 * gz + dz - v.gcz0;
        v.grad[lz * v.gpz + ly * v.gpy + lx] = (unsigned char)(1u | (code << 2));
    }
}

#define K3_CUDA_CHECK(call)                                                     \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            return false;                                                       \
        }                                                                       \
    } while (0)

template <int BX, int BY, int BZ>
bool launch_shape(const GInt::gpu::DeviceTileView3D* d_views, int nviews, int maxCw,
                  int maxCh, int maxCd, bool interior, cudaStream_t stream) {
    const int nbx = (maxCw + BX - 1) / BX, nby = (maxCh + BY - 1) / BY,
              nbz = (maxCd + BZ - 1) / BZ;
    const long long gz = (long long)nbz * nviews;
    if (nby > 65535 || gz > 65535) {
        std::fprintf(stderr, "dgrad3d: launch grid too large (%d x %d x %lld)\n", nbx, nby, gz);
        return false;
    }
    dim3 block(BX, BY, BZ);
    dim3 grid((unsigned)nbx, (unsigned)nby, (unsigned)gz);
    if (interior) dgrad3d_kernel<BX, BY, BZ, true><<<grid, block, 0, stream>>>(d_views, nbz);
    else dgrad3d_kernel<BX, BY, BZ, false><<<grid, block, 0, stream>>>(d_views, nbz);
    K3_CUDA_CHECK(cudaGetLastError());
    return true;
}

__device__ __forceinline__ unsigned PackNib(unsigned b) { return ((b >> 2) & 7u) | ((b & 1u) << 3); }

// 8 gradient bytes -> 4 packed bytes per thread
__global__ void PackNibblesKernel(const uint2* __restrict__ src, unsigned* __restrict__ dst, size_t nWords) {
    const size_t k = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nWords) return;
    const uint2 v = src[k];
    // cells 0..3 are in v.x, cells 4..7 in v.y; packed byte i holds cells 2i, 2i+1
    const unsigned o = PackNib(v.x & 0xFFu) | (PackNib((v.x >> 8) & 0xFFu) << 4) |
        (PackNib((v.x >> 16) & 0xFFu) << 8) | (PackNib((v.x >> 24) & 0xFFu) << 12) |
        (PackNib(v.y & 0xFFu) << 16) | (PackNib((v.y >> 8) & 0xFFu) << 20) |
        (PackNib((v.y >> 16) & 0xFFu) << 24) | (PackNib((v.y >> 24) & 0xFFu) << 28);
    dst[k] = o;
}

} // namespace

namespace GInt {
namespace gpu {

bool LaunchPackNibbles(const unsigned char* src, unsigned char* dst, size_t nOut, cudaStream_t stream) {
    if (nOut == 0) return true;
    if (nOut % 8) return false;
    const size_t nWords = nOut / 4;   // packed 32-bit words; each reads 8 source bytes
    const unsigned blocks = (unsigned)((nWords + 255) / 256);
    PackNibblesKernel<<<blocks, 256, 0, stream>>>((const uint2*)src, (unsigned*)dst, nWords);
    K3_CUDA_CHECK(cudaGetLastError());
    return true;
}

bool UploadTables3D(const Dgrad3DTables& t) {
    Tab3 h;
    std::memset(&h, 0, sizeof(h));
    for (int s = 0; s < 27; s++) {
        // the ordered cofacet list is only scanned for edges and quads (<= 4
        // window cofacets); the vertex's 6 cofacets are used as a mask only
        if (t.dim[s] > 3 || t.cof_count[s] > 6 || (t.dim[s] > 0 && t.cof_count[s] > 4) ||
            t.fac_count[s] > 6 || t.vert_count[s] < 1 || t.vert_count[s] > 8)
            return false;
        for (int k = 0; k < t.fac_count[s]; k++) h.fac[s] |= 1u << t.fac[s][k];
        unsigned packed = 0xFFFFFFFFu;
        for (int k = 0; k < t.cof_count[s]; k++) {
            h.cof[s] |= 1u << t.cof[s][k];
            if (k < 4) {
                packed &= ~(0xFFu << (8 * k));
                packed |= (unsigned)t.cof[s][k] << (8 * k);
            }
        }
        h.cofl[s] = packed;
        for (int k = 0; k < t.vert_count[s]; k++) h.vert[s] |= 1u << t.vert[s][k];
        // the cell's own slot is its farthest vertex (bitmask lowest() relies
        // only on the mask, but the vertex phase uses edge slot == endpoint)
        if (!(h.vert[s] >> s & 1u)) return false;
        h.dimm[t.dim[s]] |= 1u << s;
        const int d[3] = {s % 3 - 1, (s / 3) % 3 - 1, s / 9 - 1};
        for (int a = 0; a < 3; a++) {
            if (d[a] < 0) h.axneg[a] |= 1u << s;
            else if (d[a] == 0) h.axzero[a] |= 1u << s;
            else h.axpos[a] |= 1u << s;
        }
        // a lower-star cell's window facets are the ones containing the
        // vertex: at most dim of them (the 2-bit counter relies on it)
        if (t.fac_count[s] > t.dim[s]) return false;
    }
    K3_CUDA_CHECK(cudaMemcpyToSymbol(c_t3, &h, sizeof(h)));
    return true;
}

bool LaunchDgrad3DViews(const DeviceTileView3D* d_views, int nviews, int maxCw, int maxCh,
                        int maxCd, int block_shape, cudaStream_t stream, bool interior) {
    if (nviews < 1 || maxCw < 1 || maxCh < 1 || maxCd < 1) return false;
    switch (block_shape) {
        case 1: return launch_shape<16, 8, 2>(d_views, nviews, maxCw, maxCh, maxCd, interior, stream);
        case 2: return launch_shape<8, 8, 4>(d_views, nviews, maxCw, maxCh, maxCd, interior, stream);
        // face shapes of the boundary kernel: one vertex thick across x, y, z
        case kShapeFaceX: return launch_shape<1, 32, 8>(d_views, nviews, maxCw, maxCh, maxCd, false, stream);
        case kShapeFaceY: return launch_shape<32, 1, 8>(d_views, nviews, maxCw, maxCh, maxCd, false, stream);
        case kShapeFaceZ: return launch_shape<32, 8, 1>(d_views, nviews, maxCw, maxCh, maxCd, false, stream);
        default: return launch_shape<32, 4, 2>(d_views, nviews, maxCw, maxCh, maxCd, interior, stream);
    }
}

} // namespace gpu
} // namespace GInt
