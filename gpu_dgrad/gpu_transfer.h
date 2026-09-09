#ifndef GPU_TRANSFER_H
#define GPU_TRANSFER_H

// Pinned-buffer + multi-stream transfer ring (Track C1 / Stage A0).
//
// Every entry point before this stage moved data with synchronous pageable
// cudaMemcpy on the default stream (~11/9 GB/s measured) although the link
// does ~21/20 GB/s pinned. This ring gives a driver N slots, each with a
// pinned input buffer, a device input buffer, a device output buffer and a
// pinned output buffer, plus three streams (H2D, compute, D2H) and the events
// that order them. A slot is reused only after its previous D2H has landed and
// the host has consumed it (the driver calls Retire, which waits).
//
// CUDA types are used here, so this header is library-private (never included
// by a host-only TU).

#include <cuda_runtime.h>

#include <cstddef>
#include <cstdint>

namespace GInt {
namespace gpu {

struct TransferSlot {
    void* h_in = nullptr;      // pinned
    void* d_in = nullptr;
    void* d_out = nullptr;
    void* h_out = nullptr;     // pinned
    size_t cap_in = 0, cap_out = 0;
    cudaEvent_t ev_start = nullptr;   // recorded on H2D stream before the copy
    cudaEvent_t ev_h2d = nullptr;     // H2D complete
    cudaEvent_t ev_kernel = nullptr;  // kernel complete
    cudaEvent_t ev_d2h = nullptr;     // D2H complete
    bool in_flight = false;
    int64_t tag = -1;                 // driver-defined (e.g. first tile id)
};

class TransferRing {
public:
    TransferRing() {}
    ~TransferRing() { Destroy(); }
    TransferRing(const TransferRing&) = delete;
    TransferRing& operator=(const TransferRing&) = delete;

    // Allocates nslots slots with the given capacities. Returns false on any
    // CUDA failure (the object is then fully released).
    bool Create(int nslots, size_t cap_in, size_t cap_out);
    void Destroy();

    int NumSlots() const { return m_n; }
    TransferSlot& Slot(int i) { return m_slots[i]; }
    cudaStream_t StreamH2D() const { return m_s_h2d; }
    cudaStream_t StreamKernel() const { return m_s_kern; }
    cudaStream_t StreamD2H() const { return m_s_d2h; }

    // Enqueues H2D of 'bytes' from the slot's pinned input and records
    // ev_start/ev_h2d; the compute stream is made to wait on ev_h2d.
    bool BeginUpload(int slot, size_t bytes);
    // Records ev_kernel on the compute stream; the D2H stream waits on it and
    // copies 'bytes' from d_out to h_out, recording ev_d2h.
    bool EndComputeAndDownload(int slot, size_t bytes);
    // Blocks until the slot's D2H has completed. Returns the three elapsed
    // legs (ms) via the events if requested.
    bool Retire(int slot, float* h2d_ms, float* kernel_ms, float* d2h_ms);

    size_t PinnedBytes() const { return m_pinned; }
    size_t DeviceBytes() const { return m_device; }

private:
    TransferSlot* m_slots = nullptr;
    int m_n = 0;
    cudaStream_t m_s_h2d = nullptr, m_s_kern = nullptr, m_s_d2h = nullptr;
    size_t m_pinned = 0, m_device = 0;
};

} // namespace gpu
} // namespace GInt

#endif
