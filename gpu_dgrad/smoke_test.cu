// Stage 0 smoke test for the GPU discrete gradient path.
//
// 1. Queries and prints device 0 properties through the gpu_dgrad API.
// 2. Runs a trivial byte-fill kernel and verifies the result on the host
//    (exercises kernel launch + D2H on this toolchain/driver combination).
// 3. Measures pinned-memory H2D and D2H bandwidth and device-to-device
//    effective bandwidth. These calibrate the performance model in the
//    cuda-gradient master plan (Appendix D).
//
// Exit code 0 on success, nonzero on any failure.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#define CUDA_CHECK(call)                                                        \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            std::exit(1);                                                       \
        }                                                                       \
    } while (0)

__global__ void fill_pattern(unsigned char* buf, size_t n) {
    size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
    if (i < n) buf[i] = (unsigned char)((i * 31u + 7u) & 0xFF);
}

static float timed_memcpy_gbps(void* dst, const void* src, size_t bytes,
                               cudaMemcpyKind kind, int iters) {
    cudaEvent_t t0, t1;
    CUDA_CHECK(cudaEventCreate(&t0));
    CUDA_CHECK(cudaEventCreate(&t1));
    CUDA_CHECK(cudaMemcpy(dst, src, bytes, kind)); // warm-up
    CUDA_CHECK(cudaEventRecord(t0));
    for (int i = 0; i < iters; i++) CUDA_CHECK(cudaMemcpy(dst, src, bytes, kind));
    CUDA_CHECK(cudaEventRecord(t1));
    CUDA_CHECK(cudaEventSynchronize(t1));
    float ms = 0.f;
    CUDA_CHECK(cudaEventElapsedTime(&ms, t0, t1));
    CUDA_CHECK(cudaEventDestroy(t0));
    CUDA_CHECK(cudaEventDestroy(t1));
    return (float)((double)bytes * iters / (ms * 1e-3) / 1e9);
}

int main() {
    // --- 1. Device query through the API ---
    GInt::gpu::DeviceInfo info;
    if (!GInt::gpu::QueryDevice(info)) {
        std::fprintf(stderr, "FAIL: no usable CUDA device\n");
        return 1;
    }
    int runtime_ver = 0, driver_ver = 0;
    CUDA_CHECK(cudaRuntimeGetVersion(&runtime_ver));
    CUDA_CHECK(cudaDriverGetVersion(&driver_ver));
    std::printf("device        : %s\n", info.name);
    std::printf("compute cap   : %d.%d\n", info.cc_major, info.cc_minor);
    std::printf("SMs           : %d\n", info.sm_count);
    std::printf("global memory : %.2f GB\n", info.total_mem_bytes / 1e9);
    std::printf("CUDA runtime  : %d.%d\n", runtime_ver / 1000, (runtime_ver % 1000) / 10);
    std::printf("CUDA driver   : %d.%d\n", driver_ver / 1000, (driver_ver % 1000) / 10);

    // --- 2. Trivial kernel round-trip ---
    const size_t n = 64ull << 20; // 64 MB byte array
    unsigned char* d_buf = nullptr;
    CUDA_CHECK(cudaMalloc(&d_buf, n));
    fill_pattern<<<(unsigned)((n + 255) / 256), 256>>>(d_buf, n);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    unsigned char* h_buf = nullptr;
    CUDA_CHECK(cudaMallocHost(&h_buf, n));
    CUDA_CHECK(cudaMemcpy(h_buf, d_buf, n, cudaMemcpyDeviceToHost));
    for (size_t i = 0; i < n; i += 1013) { // stride-sample the verification
        unsigned char expect = (unsigned char)((i * 31u + 7u) & 0xFF);
        if (h_buf[i] != expect) {
            std::fprintf(stderr, "FAIL: kernel output mismatch at %zu: %u != %u\n",
                         i, h_buf[i], expect);
            return 1;
        }
    }
    std::printf("kernel check  : PASS (64 MB pattern round-trip)\n");

    // --- 3. Bandwidth measurements ---
    const size_t bw_bytes = 256ull << 20; // 256 MB
    const int iters = 10;
    unsigned char *h_pin = nullptr, *d_a = nullptr, *d_b = nullptr;
    CUDA_CHECK(cudaMallocHost(&h_pin, bw_bytes));
    CUDA_CHECK(cudaMalloc(&d_a, bw_bytes));
    CUDA_CHECK(cudaMalloc(&d_b, bw_bytes));

    float h2d = timed_memcpy_gbps(d_a, h_pin, bw_bytes, cudaMemcpyHostToDevice, iters);
    float d2h = timed_memcpy_gbps(h_pin, d_a, bw_bytes, cudaMemcpyDeviceToHost, iters);
    float d2d = timed_memcpy_gbps(d_b, d_a, bw_bytes, cudaMemcpyDeviceToDevice, iters);

    std::printf("pinned H2D    : %.1f GB/s\n", h2d);
    std::printf("pinned D2H    : %.1f GB/s\n", d2h);
    // A D2D copy reads and writes every byte: effective DRAM traffic is 2x the size.
    std::printf("device DRAM   : %.1f GB/s effective (D2D copy, read+write counted)\n", d2d * 2.f);

    CUDA_CHECK(cudaFreeHost(h_pin));
    CUDA_CHECK(cudaFreeHost(h_buf));
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
    CUDA_CHECK(cudaFree(d_buf));

    std::printf("SMOKE TEST    : PASS\n");
    return 0;
}
