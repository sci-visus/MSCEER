// Pinned-buffer + multi-stream transfer ring. See gpu_transfer.h.

#include "gpu_transfer.h"

#include <cstdio>

namespace {
#define RING_CUDA_CHECK(call)                                                   \
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

bool TransferRing::Create(int nslots, size_t cap_in, size_t cap_out) {
    Destroy();
    if (nslots < 1) return false;
    m_slots = new TransferSlot[nslots];
    m_n = nslots;
    bool ok = true;
    auto fail = [&](cudaError_t e) {
        if (e != cudaSuccess) {
            std::fprintf(stderr, "TransferRing::Create: %s\n", cudaGetErrorString(e));
            ok = false;
        }
        return ok;
    };
    fail(cudaStreamCreateWithFlags(&m_s_h2d, cudaStreamNonBlocking));
    fail(cudaStreamCreateWithFlags(&m_s_kern, cudaStreamNonBlocking));
    fail(cudaStreamCreateWithFlags(&m_s_d2h, cudaStreamNonBlocking));
    for (int i = 0; i < nslots && ok; i++) {
        TransferSlot& s = m_slots[i];
        s.cap_in = cap_in;
        s.cap_out = cap_out;
        if (cap_in) {
            fail(cudaMallocHost(&s.h_in, cap_in));
            fail(cudaMalloc(&s.d_in, cap_in));
        }
        if (cap_out) {
            fail(cudaMalloc(&s.d_out, cap_out));
            fail(cudaMallocHost(&s.h_out, cap_out));
        }
        fail(cudaEventCreateWithFlags(&s.ev_start, cudaEventDefault));
        fail(cudaEventCreateWithFlags(&s.ev_h2d, cudaEventDefault));
        fail(cudaEventCreateWithFlags(&s.ev_kernel, cudaEventDefault));
        fail(cudaEventCreateWithFlags(&s.ev_d2h, cudaEventDefault));
        m_pinned += cap_in + cap_out;
        m_device += cap_in + cap_out;
    }
    if (!ok) Destroy();
    return ok;
}

void TransferRing::Destroy() {
    for (int i = 0; i < m_n; i++) {
        TransferSlot& s = m_slots[i];
        if (s.in_flight) cudaEventSynchronize(s.ev_d2h);
        if (s.h_in) cudaFreeHost(s.h_in);
        if (s.d_in) cudaFree(s.d_in);
        if (s.d_out) cudaFree(s.d_out);
        if (s.h_out) cudaFreeHost(s.h_out);
        if (s.ev_start) cudaEventDestroy(s.ev_start);
        if (s.ev_h2d) cudaEventDestroy(s.ev_h2d);
        if (s.ev_kernel) cudaEventDestroy(s.ev_kernel);
        if (s.ev_d2h) cudaEventDestroy(s.ev_d2h);
    }
    delete[] m_slots;
    m_slots = nullptr;
    m_n = 0;
    if (m_s_h2d) cudaStreamDestroy(m_s_h2d);
    if (m_s_kern) cudaStreamDestroy(m_s_kern);
    if (m_s_d2h) cudaStreamDestroy(m_s_d2h);
    m_s_h2d = m_s_kern = m_s_d2h = nullptr;
    m_pinned = m_device = 0;
}

bool TransferRing::BeginUpload(int slot, size_t bytes) {
    TransferSlot& s = m_slots[slot];
    if (bytes > s.cap_in) return false;
    RING_CUDA_CHECK(cudaEventRecord(s.ev_start, m_s_h2d));
    if (bytes)
        RING_CUDA_CHECK(cudaMemcpyAsync(s.d_in, s.h_in, bytes,
                                        cudaMemcpyHostToDevice, m_s_h2d));
    RING_CUDA_CHECK(cudaEventRecord(s.ev_h2d, m_s_h2d));
    RING_CUDA_CHECK(cudaStreamWaitEvent(m_s_kern, s.ev_h2d, 0));
    s.in_flight = true;
    return true;
}

bool TransferRing::EndComputeAndDownload(int slot, size_t bytes) {
    TransferSlot& s = m_slots[slot];
    if (bytes > s.cap_out) return false;
    RING_CUDA_CHECK(cudaEventRecord(s.ev_kernel, m_s_kern));
    RING_CUDA_CHECK(cudaStreamWaitEvent(m_s_d2h, s.ev_kernel, 0));
    if (bytes)
        RING_CUDA_CHECK(cudaMemcpyAsync(s.h_out, s.d_out, bytes,
                                        cudaMemcpyDeviceToHost, m_s_d2h));
    RING_CUDA_CHECK(cudaEventRecord(s.ev_d2h, m_s_d2h));
    return true;
}

bool TransferRing::Retire(int slot, float* h2d_ms, float* kernel_ms, float* d2h_ms) {
    TransferSlot& s = m_slots[slot];
    if (!s.in_flight) return true;
    RING_CUDA_CHECK(cudaEventSynchronize(s.ev_d2h));
    if (h2d_ms) RING_CUDA_CHECK(cudaEventElapsedTime(h2d_ms, s.ev_start, s.ev_h2d));
    if (kernel_ms) RING_CUDA_CHECK(cudaEventElapsedTime(kernel_ms, s.ev_h2d, s.ev_kernel));
    if (d2h_ms) RING_CUDA_CHECK(cudaEventElapsedTime(d2h_ms, s.ev_kernel, s.ev_d2h));
    s.in_flight = false;
    return true;
}

} // namespace gpu
} // namespace GInt
