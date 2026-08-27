// 2D Robins lower-star discrete gradient, interior vertices only.
//
// Two kernels share one expansion core (a faithful re-implementation of
// MyRobinsNoalloc<TopologicalRegularGrid2D,...>::ComputeLowerStar for the
// interior single-subset case, re-indexed onto the 9 vertex-incident cells):
//
//   Stage 1 kernel: takes the CPU-precomputed max/min label byte arrays,
//     reads values/labels straight from global memory (L2 only). Kept for A/B.
//   Stage 2 kernel (fused): no label arrays at all - lower-star membership
//     and lowest-vertex queries are computed on the fly from an 18x18 tile of
//     function values staged in shared memory. This removes 2 bytes/cell of
//     labeling from storage and traffic (and its CPU precompute pass).
//
// Bit-exactness with the CPU requires replicating its iteration orders:
//   - lower-star insertion order = AdjacentCellsIterator order = ascending slot
//     (slot s = (dx+1)+3*(dy+1)); this fixes the d-cell scan order
//   - candidate cofacets are gathered in CofacetsIterator order (c_tab.cof)
//   - PickLowestCandidate keeps the EARLIEST candidate on a total-order tie
//   - the total order is IsGreater: value first, then grid vertex index; the
//     CPU max labeling breaks value ties toward the HIGHER index and the min
//     labeling toward the LOWER index (verified in gi_max_vertex_labeling.h
//     do_parallel_lines_max/min), i.e. exactly argmax/argmin under IsGreater,
//     so the fused on-the-fly computation reproduces the label semantics.
// The algorithm is comparisons-only, so matching orders and comparators makes
// the output memcmp-equal to the CPU, not just close.
//
// Per the cuda-gradient master plan (Appendix A.5) this file must not include
// any gi_* headers; the topology tables arrive generated from the real mesh.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>

namespace {

__constant__ GInt::gpu::Dgrad2DTables c_tab;

// A grid vertex under the (value, index) total order.
struct ValVn {
    float v;
    long long n;
};
// IsGreater(a,b) = a.v>b.v || (a.v==b.v && a.n>b.n); Before(a,b)=IsGreater(b,a).
__device__ __forceinline__ bool vv_greater(const ValVn& a, const ValVn& b) {
    if (a.v > b.v) return true;
    if (b.v > a.v) return false;
    return a.n > b.n;
}

// The homotopy-expansion core. CTX supplies:
//   ValVn lowvn(int s)      - lowest vertex of the lower-star cell in slot s
//   ValVn edge_other(int s) - other endpoint of the edge in slot s
template <class CTX>
__device__ void lower_star_core(const CTX& ctx, const unsigned lmask,
                                int (&pair)[9]) {
    int nmiss[9];
    for (int s = 0; s < 9; s++) {
        pair[s] = -1;
        int cnt = 0;
        if (lmask >> s & 1)
            for (int f = 0; f < c_tab.fac_count[s]; f++)
                if (lmask >> c_tab.fac[s][f] & 1) cnt++;
        nmiss[s] = cnt;
    }
    unsigned paired = 0;

    auto dec_cofacets = [&](int s) {
        for (int f = 0; f < c_tab.cof_count[s]; f++) {
            const int cf = c_tab.cof[s][f];
            if (lmask >> cf & 1) nmiss[cf]--;
        }
    };

    // ---- vertex phase: pair the vertex with its steepest lower-star edge ----
    {
        int best_s = -1;
        ValVn best{};
        for (int e = 0; e < 4; e++) {
            const int s = 2 * e + 1;               // edge slots 1,3,5,7 ascending
            if (!(lmask >> s & 1)) continue;
            const ValVn o = ctx.edge_other(s);
            // CPU: keep current iff Before(curr, other) == IsGreater(other, curr)
            if (best_s < 0 || !vv_greater(o, best)) {
                best_s = s;
                best = o;
            }
        }
        if (best_s < 0) {
            pair[4] = 4;                           // no lower edges: vertex critical
            paired |= 1u << 4;
            dec_cofacets(4);
        } else {
            pair[4] = best_s;
            pair[best_s] = 4;
            paired |= (1u << 4) | (1u << best_s);
            dec_cofacets(4);
            dec_cofacets(best_s);
        }
    }

    // ---- per-dimension expansion with restart-scan, exactly as the CPU ----
    for (int d = 0; d <= 2; d++) {
        int total = 0;
        for (int s = 0; s < 9; s++)
            if ((lmask >> s & 1) && c_tab.dim[s] == d && !(paired >> s & 1)) total++;
        int processed = 0;
        while (processed < total) {
            bool advanced = false;
            for (int s = 0; s < 9 && !advanced; s++) {
                if (!(lmask >> s & 1) || c_tab.dim[s] != d || (paired >> s & 1))
                    continue;
                int cand[4], ncand = 0;
                for (int f = 0; f < c_tab.cof_count[s]; f++) {
                    const int cf = c_tab.cof[s][f];
                    if ((lmask >> cf & 1) && nmiss[cf] == 1) cand[ncand++] = cf;
                }
                if (ncand > 0) {
                    int best = cand[0];
                    if (ncand > 1) {
                        ValVn blv = ctx.lowvn(best);
                        for (int i = 1; i < ncand; i++) {
                            const ValVn olv = ctx.lowvn(cand[i]);
                            // switch iff Before(olv, blv) == IsGreater(blv, olv)
                            if (vv_greater(blv, olv)) {
                                blv = olv;
                                best = cand[i];
                            }
                        }
                    }
                    pair[s] = best;
                    pair[best] = s;
                    paired |= (1u << s) | (1u << best);
                    dec_cofacets(s);
                    dec_cofacets(best);
                    processed++;
                    advanced = true;               // restart the d-cell scan
                }
            }
            if (!advanced) {
                int best = -1;
                ValVn blv{};
                for (int s = 0; s < 9; s++) {
                    if (!(lmask >> s & 1) || c_tab.dim[s] != d || (paired >> s & 1))
                        continue;
                    const ValVn olv = ctx.lowvn(s);
                    if (best < 0 || vv_greater(blv, olv)) {
                        best = s;
                        blv = olv;
                    }
                }
                pair[best] = best;                 // lowest unpaired becomes critical
                paired |= 1u << best;
                dec_cofacets(best);
                processed++;
            }
        }
    }
}

// Compose each GradBitfield byte once and store it (assigned=1, pair=code).
__device__ __forceinline__ void scatter_pairs(const unsigned lmask,
                                              const int (&pair)[9],
                                              const long long (&cid)[9],
                                              unsigned char* __restrict__ grad) {
    for (int s = 0; s < 9; s++) {
        if (!(lmask >> s & 1)) continue;
        const int p = pair[s];
        unsigned code;
        if (p == s) {
            code = 7;
        } else {
            const int ddx = (p % 3) - (s % 3);
            const int ddy = (p / 3) - (s / 3);
            code = (ddx == 1) ? c_tab.dir_px
                 : (ddx == -1) ? c_tab.dir_mx
                 : (ddy == 1) ? c_tab.dir_py
                 : c_tab.dir_my;
        }
        grad[cid[s]] = (unsigned char)(1u | (code << 2));
    }
}

// ---------------------------------------------------------------------------
// Stage 1 kernel: CPU-precomputed labels, global-memory reads.
// ---------------------------------------------------------------------------

struct LabelCtx {
    const float* vals;
    const unsigned char* minb;
    const long long* cid;
    int gx, gy, X;

    __device__ ValVn lowvn(int s) const {
        const int m = minb[cid[s]];
        const long long vx = (2LL * gx + (s % 3 - 1) + (m % 3 - 1)) >> 1;
        const long long vy = (2LL * gy + (s / 3 - 1) + (m / 3 - 1)) >> 1;
        const long long vn = vx + vy * (long long)X;
        return ValVn{vals[vn], vn};
    }
    __device__ ValVn edge_other(int s) const {
        const long long vn =
            (long long)(gx + (s % 3 - 1)) + (long long)(gy + (s / 3 - 1)) * X;
        return ValVn{vals[vn], vn};
    }
};

__global__ void dgrad2d_interior_kernel(const float* __restrict__ vals,
                                        const unsigned char* __restrict__ maxb,
                                        const unsigned char* __restrict__ minb,
                                        unsigned char* __restrict__ grad,
                                        int X, int Y) {
    const int gx = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int gy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (gx > X - 2 || gy > Y - 2) return;

    const long long W = 2LL * X - 1;
    const long long vmesh = 2LL * gx + 2LL * gy * W;

    // Lower-star membership: cell in slot s has highest vertex == center iff
    // its max byte points back at the center (offset symmetry: byte == 8-s).
    long long cid[9];
    unsigned lmask = 0;
    for (int s = 0; s < 9; s++) {
        const long long c = vmesh + (s % 3 - 1) + (long long)(s / 3 - 1) * W;
        cid[s] = c;
        if (maxb[c] == (unsigned char)(8 - s)) lmask |= 1u << s;
    }

    const LabelCtx ctx{vals, minb, cid, gx, gy, X};
    int pair[9];
    lower_star_core(ctx, lmask, pair);
    scatter_pairs(lmask, pair, cid, grad);
}

// ---------------------------------------------------------------------------
// Stage 2 kernel: label fusion. Block = 16x16 interior vertices staging an
// 18x18 tile of function values in shared memory; membership and lowest-vertex
// queries computed on the fly under the (value, index) total order.
// ---------------------------------------------------------------------------

constexpr int BX = 16, BY = 16;
constexpr int TX = BX + 2, TY = BY + 2;   // tile with 1-vertex halo
constexpr int TPITCH = TX + 1;            // pad a column to stagger banks

struct FusedCtx {
    const float (*sh)[TPITCH]; // shared tile
    int tx, ty;                // this thread's vertex at sh[ty+1][tx+1]
    int gx, gy, X;

    __device__ __forceinline__ ValVn vert(int ax, int ay) const {
        return ValVn{sh[ty + 1 + ay][tx + 1 + ax],
                     (long long)(gx + ax) + (long long)(gy + ay) * X};
    }
    // lowest vertex of the cell in slot s = argmin over its 1/2/4 vertices
    // under the total order (== the CPU min labeling's tie semantics)
    __device__ ValVn lowvn(int s) const {
        const int dx = s % 3 - 1, dy = s / 3 - 1;
        const int x0 = dx < 0 ? -1 : 0, x1 = dx > 0 ? 1 : 0;
        const int y0 = dy < 0 ? -1 : 0, y1 = dy > 0 ? 1 : 0;
        ValVn best = vert(x0, y0);
        for (int ay = y0; ay <= y1; ay++)
            for (int ax = x0; ax <= x1; ax++) {
                if (ax == x0 && ay == y0) continue;
                const ValVn o = vert(ax, ay);
                if (vv_greater(best, o)) best = o;
            }
        return best;
    }
    __device__ ValVn edge_other(int s) const {
        return vert(s % 3 - 1, s / 3 - 1);
    }
};

__global__ void dgrad2d_fused_kernel(const float* __restrict__ vals,
                                     unsigned char* __restrict__ grad,
                                     int X, int Y) {
    __shared__ float sh[TY][TPITCH];

    const int ox = blockIdx.x * BX;        // tile origin = vertex (ox, oy)
    const int oy = blockIdx.y * BY;

    // cooperative tile load (all threads participate, even ones with no vertex)
    for (int i = threadIdx.y * BX + threadIdx.x; i < TX * TY; i += BX * BY) {
        const int lx = i % TX, ly = i / TX;
        const int gxx = ox + lx, gyy = oy + ly;
        sh[ly][lx] = (gxx < X && gyy < Y) ? vals[(long long)gyy * X + gxx] : 0.f;
    }
    __syncthreads();

    const int gx = ox + threadIdx.x + 1;
    const int gy = oy + threadIdx.y + 1;
    if (gx > X - 2 || gy > Y - 2) return;

    const int tx = threadIdx.x, ty = threadIdx.y;
    const FusedCtx ctx{sh, tx, ty, gx, gy, X};

    // membership: center is the strict maximum of the cell's vertex set
    const ValVn center = ctx.vert(0, 0);
    const long long W = 2LL * X - 1;
    const long long vmesh = 2LL * gx + 2LL * gy * W;
    long long cid[9];
    unsigned lmask = 1u << 4;              // the vertex itself
    for (int s = 0; s < 9; s++) {
        const int dx = s % 3 - 1, dy = s / 3 - 1;
        cid[s] = vmesh + dx + (long long)dy * W;
        if (s == 4) continue;
        const int x0 = dx < 0 ? -1 : 0, x1 = dx > 0 ? 1 : 0;
        const int y0 = dy < 0 ? -1 : 0, y1 = dy > 0 ? 1 : 0;
        bool in = true;
        for (int ay = y0; ay <= y1 && in; ay++)
            for (int ax = x0; ax <= x1 && in; ax++) {
                if (ax == 0 && ay == 0) continue;
                if (!vv_greater(center, ctx.vert(ax, ay))) in = false;
            }
        if (in) lmask |= 1u << s;
    }

    int pair[9];
    lower_star_core(ctx, lmask, pair);
    scatter_pairs(lmask, pair, cid, grad);
}

#define API_CUDA_CHECK(call)                                                    \
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

namespace GInt {
namespace gpu {

bool ComputeDiscreteGradient2D(const float* values, int64_t X, int64_t Y,
                               const uint8_t* max_labels,
                               const uint8_t* min_labels,
                               const Dgrad2DTables& tables,
                               uint8_t* out_grad,
                               bool interior_only,
                               Dgrad2DTimings* timings) {
    if (!interior_only) {
        std::fprintf(stderr,
                     "ComputeDiscreteGradient2D: full-domain mode lands in Stage 3\n");
        return false;
    }
    if (X < 3 || Y < 3) return false;

    const int64_t n_verts = X * Y;
    const int64_t n_cells = (2 * X - 1) * (2 * Y - 1);
    const bool fused = (max_labels == nullptr);

    float* d_vals = nullptr;
    uint8_t *d_max = nullptr, *d_min = nullptr, *d_grad = nullptr;
    API_CUDA_CHECK(cudaMalloc(&d_vals, n_verts * sizeof(float)));
    API_CUDA_CHECK(cudaMalloc(&d_grad, n_cells));
    if (!fused) {
        API_CUDA_CHECK(cudaMalloc(&d_max, n_cells));
        API_CUDA_CHECK(cudaMalloc(&d_min, n_cells));
    }

    cudaEvent_t e0, e1, e2, e3;
    API_CUDA_CHECK(cudaEventCreate(&e0));
    API_CUDA_CHECK(cudaEventCreate(&e1));
    API_CUDA_CHECK(cudaEventCreate(&e2));
    API_CUDA_CHECK(cudaEventCreate(&e3));

    API_CUDA_CHECK(cudaEventRecord(e0));
    API_CUDA_CHECK(cudaMemcpy(d_vals, values, n_verts * sizeof(float),
                              cudaMemcpyHostToDevice));
    if (!fused) {
        API_CUDA_CHECK(cudaMemcpy(d_max, max_labels, n_cells, cudaMemcpyHostToDevice));
        API_CUDA_CHECK(cudaMemcpy(d_min, min_labels, n_cells, cudaMemcpyHostToDevice));
    }
    API_CUDA_CHECK(cudaMemset(d_grad, 0, n_cells));
    API_CUDA_CHECK(cudaMemcpyToSymbol(c_tab, &tables, sizeof(tables)));
    API_CUDA_CHECK(cudaEventRecord(e1));

    if (fused) {
        dim3 block(BX, BY);
        dim3 grid((unsigned)((X - 2 + BX - 1) / BX),
                  (unsigned)((Y - 2 + BY - 1) / BY));
        dgrad2d_fused_kernel<<<grid, block>>>(d_vals, d_grad, (int)X, (int)Y);
    } else {
        dim3 block(16, 16);
        dim3 grid((unsigned)((X - 2 + block.x - 1) / block.x),
                  (unsigned)((Y - 2 + block.y - 1) / block.y));
        dgrad2d_interior_kernel<<<grid, block>>>(d_vals, d_max, d_min, d_grad,
                                                 (int)X, (int)Y);
    }
    API_CUDA_CHECK(cudaGetLastError());
    API_CUDA_CHECK(cudaEventRecord(e2));

    API_CUDA_CHECK(cudaMemcpy(out_grad, d_grad, n_cells, cudaMemcpyDeviceToHost));
    API_CUDA_CHECK(cudaEventRecord(e3));
    API_CUDA_CHECK(cudaEventSynchronize(e3));

    if (timings) {
        API_CUDA_CHECK(cudaEventElapsedTime(&timings->h2d_ms, e0, e1));
        API_CUDA_CHECK(cudaEventElapsedTime(&timings->kernel_ms, e1, e2));
        API_CUDA_CHECK(cudaEventElapsedTime(&timings->d2h_ms, e2, e3));
    }

    cudaEventDestroy(e0);
    cudaEventDestroy(e1);
    cudaEventDestroy(e2);
    cudaEventDestroy(e3);
    cudaFree(d_vals);
    if (d_max) cudaFree(d_max);
    if (d_min) cudaFree(d_min);
    cudaFree(d_grad);
    return true;
}

bool ComputeDiscreteGradient2DFused(const float* values, int64_t X, int64_t Y,
                                    const Dgrad2DTables& tables,
                                    uint8_t* out_grad,
                                    bool interior_only,
                                    Dgrad2DTimings* timings) {
    return ComputeDiscreteGradient2D(values, X, Y, nullptr, nullptr, tables,
                                     out_grad, interior_only, timings);
}

} // namespace gpu
} // namespace GInt
