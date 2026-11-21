#ifndef GPU_DEVICE_H
#define GPU_DEVICE_H

// Device query and host page-locking helpers of the gpu_dgrad library. Plain C++.

#include <cstddef>
#include <cstdint>

namespace GInt {
namespace gpu {

struct DeviceInfo {
    char name[256];
    int sm_count;
    long long total_mem_bytes;
    int cc_major;
    int cc_minor;
};

// Fills 'out' with the properties of CUDA device 0.
// Returns false if no usable CUDA device is present.
bool QueryDevice(DeviceInfo& out);

// Free / total device memory in bytes for the current device. Either output
// may be null. Returns false if there is no usable device. Exists so leak
// gates can be written in a plain host TU without including <cuda_runtime.h>.
bool QueryDeviceMemory(int64_t* free_bytes, int64_t* total_bytes);

// Page-lock / release a host range (cudaHostRegister); IsHostMemoryPinned
// tells whether [ptr, ptr+bytes) is page-locked memory known to CUDA.
bool PinHostMemory(void* ptr, size_t bytes);
bool UnpinHostMemory(void* ptr);
bool IsHostMemoryPinned(const void* ptr, size_t bytes);

} // namespace gpu
} // namespace GInt

#endif
