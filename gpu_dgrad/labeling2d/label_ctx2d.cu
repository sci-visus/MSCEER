// Persistent device label context for 2D per-pixel relabeling.
//
// The base labeling assigns every pixel a compact base-region id, and it is a
// function of the discrete gradient alone: it does not change when the user
// moves the persistence threshold. What changes is the (tiny) map from base
// region to whatever the caller wants to display - a living node id, a class
// id, a color. Every other entry point in this library is
// alloc -> upload -> compute -> download -> free, which for this workload
// would re-upload the largest buffer in the problem on every interactive
// slider tick, at ~20 GB/s pinned (see CPU_TO_GPU_GUIDE.md) - far more than
// the gather itself costs.
//
// So this is the one persistent handle: the base labeling stays resident, a
// relabel uploads m int32 and runs one fully coalesced gather. The context
// holds independent buffers for the ascending and descending directions
// because callers routinely paint both at one threshold.
//
// The gather is deliberately bounds-checked (a compact id outside [0, m)
// yields -1). It cannot happen when the remap is sized by the region count,
// but the check is one predicated compare against a value already in a
// register, it makes the contract total, and it keeps a caller bug from
// reading out of bounds on the device where the failure would be silent.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>

namespace {

// out[i] = base[i] < 0 ? -1 : remap[base[i]], -1 also for out-of-range ids.
__global__ void gather_relabel(const int32_t* __restrict__ base,
                               const int32_t* __restrict__ remap,
                               int32_t* __restrict__ out, long long n,
                               int32_t m) {
    const long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i >= n) return;
    const int32_t b = base[i];
    out[i] = (b < 0 || b >= m) ? -1 : remap[b];
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

// Direction index: 0 = ascending, 1 = descending. Buffers are lazily
// allocated per direction, so a caller that only ever paints one direction
// pays for one.
struct LabelCtx2D {
    int64_t X;
    int64_t Y;
    long long n;
    int32_t* d_base[2];
    int32_t* d_out[2];
    int32_t* d_remap[2];
    int32_t remap_cap[2];
};

LabelCtx2D* CreateLabelCtx2D(int64_t X, int64_t Y) {
    if (X <= 0 || Y <= 0) return NULL;
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) return NULL;
    // Forces context creation now, so a machine with a device that cannot be
    // initialized fails here rather than on the first relabel.
    if (cudaFree(0) != cudaSuccess) return NULL;

    LabelCtx2D* ctx = new LabelCtx2D();
    ctx->X = X;
    ctx->Y = Y;
    ctx->n = X * Y;
    for (int d = 0; d < 2; d++) {
        ctx->d_base[d] = NULL;
        ctx->d_out[d] = NULL;
        ctx->d_remap[d] = NULL;
        ctx->remap_cap[d] = 0;
    }
    return ctx;
}

void DestroyLabelCtx2D(LabelCtx2D* ctx) {
    if (!ctx) return;
    for (int d = 0; d < 2; d++) {
        if (ctx->d_base[d]) cudaFree(ctx->d_base[d]);
        if (ctx->d_out[d]) cudaFree(ctx->d_out[d]);
        if (ctx->d_remap[d]) cudaFree(ctx->d_remap[d]);
    }
    delete ctx;
}

bool CtxSetBaseLabels(LabelCtx2D* ctx, bool ascending, const int32_t* base_compact) {
    if (!ctx || !base_compact) return false;
    const int d = ascending ? 0 : 1;
    const size_t bytes = (size_t)ctx->n * sizeof(int32_t);
    if (!ctx->d_base[d]) API_CUDA_CHECK(cudaMalloc(&ctx->d_base[d], bytes));
    API_CUDA_CHECK(cudaMemcpy(ctx->d_base[d], base_compact, bytes,
                              cudaMemcpyHostToDevice));
    return true;
}

const void* CtxBaseLabelsDevice(const LabelCtx2D* ctx, bool ascending) {
    if (!ctx) return NULL;
    return ctx->d_base[ascending ? 0 : 1];
}

namespace {

// Shared body of CtxRelabel / CtxRelabelDevice: upload remap, launch, and
// optionally download. Splitting it here keeps the two entry points from
// drifting apart; the event bracketing matches the rest of the library
// (e0 | H2D | e1 | kernel | e2 | D2H | e3).
bool ctx_relabel_impl(LabelCtx2D* ctx, bool ascending, const int32_t* remap,
                      int32_t m, int32_t* out_labels_host,
                      Dgrad2DTimings* timings) {
    if (!ctx || m < 0) return false;
    if (m > 0 && !remap) return false;
    const int d = ascending ? 0 : 1;
    if (!ctx->d_base[d]) return false; // no base labeling uploaded yet

    const size_t out_bytes = (size_t)ctx->n * sizeof(int32_t);
    if (!ctx->d_out[d]) API_CUDA_CHECK(cudaMalloc(&ctx->d_out[d], out_bytes));
    if (m > ctx->remap_cap[d]) {
        if (ctx->d_remap[d]) cudaFree(ctx->d_remap[d]);
        ctx->d_remap[d] = NULL;
        ctx->remap_cap[d] = 0;
        API_CUDA_CHECK(cudaMalloc(&ctx->d_remap[d], (size_t)m * sizeof(int32_t)));
        ctx->remap_cap[d] = m;
    }

    cudaEvent_t e0, e1, e2, e3;
    API_CUDA_CHECK(cudaEventCreate(&e0));
    API_CUDA_CHECK(cudaEventCreate(&e1));
    API_CUDA_CHECK(cudaEventCreate(&e2));
    API_CUDA_CHECK(cudaEventCreate(&e3));

    API_CUDA_CHECK(cudaEventRecord(e0));
    if (m > 0)
        API_CUDA_CHECK(cudaMemcpy(ctx->d_remap[d], remap,
                                  (size_t)m * sizeof(int32_t),
                                  cudaMemcpyHostToDevice));
    API_CUDA_CHECK(cudaEventRecord(e1));

    const int TPB = 256;
    gather_relabel<<<(unsigned)((ctx->n + TPB - 1) / TPB), TPB>>>(
        ctx->d_base[d], ctx->d_remap[d], ctx->d_out[d], ctx->n, m);
    API_CUDA_CHECK(cudaGetLastError());
    API_CUDA_CHECK(cudaEventRecord(e2));

    if (out_labels_host)
        API_CUDA_CHECK(cudaMemcpy(out_labels_host, ctx->d_out[d], out_bytes,
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
    return true;
}

} // namespace

bool CtxRelabel(LabelCtx2D* ctx, bool ascending, const int32_t* remap, int32_t m,
                int32_t* out_labels_host, Dgrad2DTimings* timings) {
    if (!out_labels_host) return false;
    return ctx_relabel_impl(ctx, ascending, remap, m, out_labels_host, timings);
}

bool CtxRelabelDevice(LabelCtx2D* ctx, bool ascending, const int32_t* remap,
                      int32_t m, const void** out_dev_labels,
                      Dgrad2DTimings* timings) {
    if (!out_dev_labels) return false;
    if (!ctx_relabel_impl(ctx, ascending, remap, m, NULL, timings)) return false;
    *out_dev_labels = ctx->d_out[ascending ? 0 : 1];
    return true;
}

bool QueryDeviceMemory(int64_t* free_bytes, int64_t* total_bytes) {
    size_t f = 0, t = 0;
    if (cudaMemGetInfo(&f, &t) != cudaSuccess) return false;
    if (free_bytes) *free_bytes = (int64_t)f;
    if (total_bytes) *total_bytes = (int64_t)t;
    return true;
}

} // namespace gpu
} // namespace GInt
