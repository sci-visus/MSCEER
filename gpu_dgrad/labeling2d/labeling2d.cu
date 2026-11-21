// Stage 4: flow-terminal labeling of the 2D discrete gradient by pointer
// doubling over the offset-doubled successor.
//
// The discrete flow operator is cell -> pair -> other facet/cofacet. Because
// like-dimensional continuation in the doubled grid is always a straight
// two-step, the intermediate edge never needs to be visited:
//   descending (vertex) flow:  v paired with e = v + d  =>  succ(v) = v + 2d
//   ascending  (quad)  flow:   q paired with e = q + d  =>  succ(q) = q + 2d
// One gradient byte read and an add per step - this is the ALU saving the
// offset-doubling observation buys. On the vertex lattice 2d becomes a +-1
// step in x or y (and likewise on the quad lattice), so successors are plain
// int32 lattice indices and the terminal labeling is in-place pointer
// doubling: L[i] = L[L[i]] until fixpoint. In-place halving is race-tolerant
// (any interleaving still reads an ancestor) and the forest roots make the
// result deterministic.
//
// Boundary safety needs no masks here: a vertex's paired edge lies in the
// vertex's own boundary stratum, so its other endpoint is in-domain; and all
// quads have boundaryValue 0, so a quad's paired edge is interior and the
// opposite cofacet always exists.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>

namespace {

struct StepTable {
    int step[8]; // lattice step per 3-bit pair code; code 7 = critical/terminal
};

__global__ void init_vertex_succ(const unsigned char* __restrict__ grad,
                                 int32_t* __restrict__ L, int X, int Y,
                                 StepTable st) {
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i >= (long long)X * Y) return;
    const long long gx = i % X, gy = i / X;
    const long long W = 2LL * X - 1;
    const unsigned code = (grad[2 * gx + 2 * gy * W] >> 2) & 7u;
    L[i] = (code == 7u) ? (int32_t)i : (int32_t)(i + st.step[code]);
}

__global__ void init_quad_succ(const unsigned char* __restrict__ grad,
                               int32_t* __restrict__ L, int X, int Y,
                               StepTable st) {
    const int QX = X - 1, QY = Y - 1;
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i >= (long long)QX * QY) return;
    const long long qx = i % QX, qy = i / QX;
    const long long W = 2LL * X - 1;
    const unsigned code = (grad[(2 * qx + 1) + (2 * qy + 1) * W] >> 2) & 7u;
    L[i] = (code == 7u) ? (int32_t)i : (int32_t)(i + st.step[code]);
}

__global__ void pointer_jump(int32_t* __restrict__ L, long long n,
                             int* __restrict__ changed) {
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i >= n) return;
    const int32_t a = L[i];
    const int32_t b = L[a];
    if (b != a) {
        L[i] = b;
        *changed = 1; // benign race: any lane setting it is enough
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

bool run_doubling(int32_t* d_L, long long n, int* d_changed) {
    const int TPB = 256;
    const unsigned nb = (unsigned)((n + TPB - 1) / TPB);
    int h_changed = 1;
    while (h_changed) {
        API_CUDA_CHECK(cudaMemset(d_changed, 0, sizeof(int)));
        pointer_jump<<<nb, TPB>>>(d_L, n, d_changed);
        API_CUDA_CHECK(cudaGetLastError());
        API_CUDA_CHECK(cudaMemcpy(&h_changed, d_changed, sizeof(int),
                                  cudaMemcpyDeviceToHost));
    }
    return true;
}

} // namespace

namespace GInt {
namespace gpu {

bool Label2DFlowTerminals(const uint8_t* grad, int64_t X, int64_t Y,
                          const Dgrad2DTables& tables,
                          int32_t* out_vertex_min, int32_t* out_quad_max,
                          Dgrad2DTimings* timings) {
    if (X < 2 || Y < 2) return false;
    const long long n_cells = (2 * X - 1) * (2 * Y - 1);
    const long long n_verts = (long long)X * Y;
    const long long n_quads = (long long)(X - 1) * (Y - 1);

    StepTable vstep{}, qstep{};
    vstep.step[tables.dir_px] = 1;
    vstep.step[tables.dir_mx] = -1;
    vstep.step[tables.dir_py] = (int)X;
    vstep.step[tables.dir_my] = -(int)X;
    qstep.step[tables.dir_px] = 1;
    qstep.step[tables.dir_mx] = -1;
    qstep.step[tables.dir_py] = (int)(X - 1);
    qstep.step[tables.dir_my] = -(int)(X - 1);

    uint8_t* d_grad = nullptr;
    int32_t *d_LV = nullptr, *d_LQ = nullptr;
    int* d_changed = nullptr;
    API_CUDA_CHECK(cudaMalloc(&d_grad, n_cells));
    API_CUDA_CHECK(cudaMalloc(&d_changed, sizeof(int)));
    if (out_vertex_min)
        API_CUDA_CHECK(cudaMalloc(&d_LV, n_verts * sizeof(int32_t)));
    if (out_quad_max)
        API_CUDA_CHECK(cudaMalloc(&d_LQ, n_quads * sizeof(int32_t)));

    cudaEvent_t e0, e1, e2, e3;
    API_CUDA_CHECK(cudaEventCreate(&e0));
    API_CUDA_CHECK(cudaEventCreate(&e1));
    API_CUDA_CHECK(cudaEventCreate(&e2));
    API_CUDA_CHECK(cudaEventCreate(&e3));

    API_CUDA_CHECK(cudaEventRecord(e0));
    API_CUDA_CHECK(cudaMemcpy(d_grad, grad, n_cells, cudaMemcpyHostToDevice));
    API_CUDA_CHECK(cudaEventRecord(e1));

    const int TPB = 256;
    bool ok = true;
    if (d_LV) {
        init_vertex_succ<<<(unsigned)((n_verts + TPB - 1) / TPB), TPB>>>(
            d_grad, d_LV, (int)X, (int)Y, vstep);
        API_CUDA_CHECK(cudaGetLastError());
        ok = ok && run_doubling(d_LV, n_verts, d_changed);
    }
    if (ok && d_LQ) {
        init_quad_succ<<<(unsigned)((n_quads + TPB - 1) / TPB), TPB>>>(
            d_grad, d_LQ, (int)X, (int)Y, qstep);
        API_CUDA_CHECK(cudaGetLastError());
        ok = ok && run_doubling(d_LQ, n_quads, d_changed);
    }
    API_CUDA_CHECK(cudaEventRecord(e2));

    if (ok && d_LV)
        API_CUDA_CHECK(cudaMemcpy(out_vertex_min, d_LV,
                                  n_verts * sizeof(int32_t),
                                  cudaMemcpyDeviceToHost));
    if (ok && d_LQ)
        API_CUDA_CHECK(cudaMemcpy(out_quad_max, d_LQ, n_quads * sizeof(int32_t),
                                  cudaMemcpyDeviceToHost));
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
    cudaFree(d_grad);
    cudaFree(d_changed);
    if (d_LV) cudaFree(d_LV);
    if (d_LQ) cudaFree(d_LQ);
    return ok;
}

} // namespace gpu
} // namespace GInt
