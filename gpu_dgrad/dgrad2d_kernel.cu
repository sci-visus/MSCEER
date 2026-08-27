// Stage 1: 2D Robins lower-star discrete gradient, interior vertices only.
//
// One thread per interior vertex. This is a faithful re-implementation of
// MyRobinsNoalloc<TopologicalRegularGrid2D,...>::ComputeLowerStar for the
// interior case (single boundary-class subset), re-indexed onto the 9 cells
// incident to the vertex. Bit-exactness with the CPU requires replicating its
// iteration orders exactly:
//   - lower-star insertion order = AdjacentCellsIterator order = ascending slot
//     (slot s = (dx+1)+3*(dy+1)); this fixes the d-cell scan order below
//   - candidate cofacets are gathered in CofacetsIterator order (c_tab.cof)
//   - PickLowestCandidate keeps the EARLIEST candidate on a total-order tie
//   - the total order is IsGreater: value first, then grid vertex index
// The algorithm is comparisons-only (no FP arithmetic), so matching orders
// and comparators makes the output memcmp-equal to the CPU, not just close.
//
// Per the cuda-gradient master plan (Appendix A.5) this file must not include
// any gi_* headers; the topology tables arrive generated from the real mesh.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>

namespace {

__constant__ GInt::gpu::Dgrad2DTables c_tab;

// The (value, index) total order on grid vertices:
// RegularGridBilinearFunction::IsGreater(a,b) = v[a]>v[b] || (v[a]==v[b] && a>b).
// Before(a,b) in the labeling is IsGreater(b,a).
__device__ __forceinline__ bool is_greater(const float* __restrict__ v,
                                           long long a, long long b) {
    const float va = v[a], vb = v[b];
    if (va > vb) return true;
    if (vb > va) return false;
    return a > b;
}

__global__ void dgrad2d_interior_kernel(const float* __restrict__ vals,
                                        const unsigned char* __restrict__ maxb,
                                        const unsigned char* __restrict__ minb,
                                        unsigned char* __restrict__ grad,
                                        int X, int Y) {
    const int gx = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int gy = blockIdx.y * blockDim.y + threadIdx.y + 1;
    if (gx > X - 2 || gy > Y - 2) return;

    const long long W = 2LL * X - 1;               // doubled-grid row stride
    const long long vmesh = 2LL * gx + 2LL * gy * W;

    // Lower-star membership: cell in slot s has highest vertex == center iff
    // its max byte points back at the center. Offsets are symmetric across the
    // table (off9[8-s] == -off9[s]), so that byte is exactly 8-s.
    long long cid[9];
    unsigned lmask = 0;
    for (int s = 0; s < 9; s++) {
        const long long c = vmesh + (s % 3 - 1) + (long long)(s / 3 - 1) * W;
        cid[s] = c;
        if (maxb[c] == (unsigned char)(8 - s)) lmask |= 1u << s;
    }

    int nmiss[9];
    int pair[9];
    for (int s = 0; s < 9; s++) {
        pair[s] = -1;
        int cnt = 0;
        if (lmask >> s & 1)
            for (int f = 0; f < c_tab.fac_count[s]; f++)
                if (lmask >> c_tab.fac[s][f] & 1) cnt++;
        nmiss[s] = cnt;
    }
    unsigned paired = 0;

    // Lowest vertex (grid number) of the cell in slot s, via the min labels.
    auto lowvn = [&](int s) -> long long {
        const int m = minb[cid[s]];
        const long long vx = (2LL * gx + (s % 3 - 1) + (m % 3 - 1)) >> 1;
        const long long vy = (2LL * gy + (s / 3 - 1) + (m / 3 - 1)) >> 1;
        return vx + vy * (long long)X;
    };
    auto dec_cofacets = [&](int s) {
        for (int f = 0; f < c_tab.cof_count[s]; f++) {
            const int cf = c_tab.cof[s][f];
            if (lmask >> cf & 1) nmiss[cf]--;
        }
    };

    // ---- vertex phase: pair the vertex with its steepest lower-star edge ----
    {
        int best_s = -1;
        long long best_vn = -1;
        for (int e = 0; e < 4; e++) {
            const int s = 2 * e + 1;               // edge slots 1,3,5,7 ascending
            if (!(lmask >> s & 1)) continue;
            const long long ovn =
                (long long)(gx + (s % 3 - 1)) + (long long)(gy + (s / 3 - 1)) * X;
            // CPU: keep current iff Before(curr, other) == IsGreater(other, curr)
            if (best_s < 0 || !is_greater(vals, ovn, best_vn)) {
                best_s = s;
                best_vn = ovn;
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
                        long long blv = lowvn(best);
                        for (int i = 1; i < ncand; i++) {
                            const long long olv = lowvn(cand[i]);
                            // switch iff Before(olv, blv) == IsGreater(blv, olv)
                            if (is_greater(vals, blv, olv)) {
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
                long long blv = 0;
                for (int s = 0; s < 9; s++) {
                    if (!(lmask >> s & 1) || c_tab.dim[s] != d || (paired >> s & 1))
                        continue;
                    const long long olv = lowvn(s);
                    if (best < 0 || is_greater(vals, blv, olv)) {
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

    // ---- scatter: compose each GradBitfield byte once and store it ----
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
        grad[cid[s]] = (unsigned char)(1u | (code << 2)); // assigned=1, pair=code
    }
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

    float* d_vals = nullptr;
    uint8_t *d_max = nullptr, *d_min = nullptr, *d_grad = nullptr;
    API_CUDA_CHECK(cudaMalloc(&d_vals, n_verts * sizeof(float)));
    API_CUDA_CHECK(cudaMalloc(&d_max, n_cells));
    API_CUDA_CHECK(cudaMalloc(&d_min, n_cells));
    API_CUDA_CHECK(cudaMalloc(&d_grad, n_cells));

    cudaEvent_t e0, e1, e2, e3;
    API_CUDA_CHECK(cudaEventCreate(&e0));
    API_CUDA_CHECK(cudaEventCreate(&e1));
    API_CUDA_CHECK(cudaEventCreate(&e2));
    API_CUDA_CHECK(cudaEventCreate(&e3));

    API_CUDA_CHECK(cudaEventRecord(e0));
    API_CUDA_CHECK(cudaMemcpy(d_vals, values, n_verts * sizeof(float),
                              cudaMemcpyHostToDevice));
    API_CUDA_CHECK(cudaMemcpy(d_max, max_labels, n_cells, cudaMemcpyHostToDevice));
    API_CUDA_CHECK(cudaMemcpy(d_min, min_labels, n_cells, cudaMemcpyHostToDevice));
    API_CUDA_CHECK(cudaMemset(d_grad, 0, n_cells));
    API_CUDA_CHECK(cudaMemcpyToSymbol(c_tab, &tables, sizeof(tables)));
    API_CUDA_CHECK(cudaEventRecord(e1));

    dim3 block(16, 16);
    dim3 grid((unsigned)((X - 2 + block.x - 1) / block.x),
              (unsigned)((Y - 2 + block.y - 1) / block.y));
    dgrad2d_interior_kernel<<<grid, block>>>(d_vals, d_max, d_min, d_grad,
                                             (int)X, (int)Y);
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
    cudaFree(d_max);
    cudaFree(d_min);
    cudaFree(d_grad);
    return true;
}

} // namespace gpu
} // namespace GInt
