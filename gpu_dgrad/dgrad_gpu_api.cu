#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>

namespace GInt {
namespace gpu {

bool QueryDevice(DeviceInfo& out) {
    std::memset(&out, 0, sizeof(out));
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) return false;
    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) return false;
    std::snprintf(out.name, sizeof(out.name), "%s", prop.name);
    out.sm_count = prop.multiProcessorCount;
    out.total_mem_bytes = static_cast<long long>(prop.totalGlobalMem);
    out.cc_major = prop.major;
    out.cc_minor = prop.minor;
    return true;
}

bool ComputeDiscreteGradient3D(const float*, int64_t, int64_t, int64_t,
                               const uint8_t*, const uint8_t*, uint8_t*, bool) {
    // Implemented in Stage 1 of the cuda-gradient plan.
    std::fprintf(stderr, "ComputeDiscreteGradient3D: not implemented yet (Stage 1)\n");
    return false;
}

} // namespace gpu
} // namespace GInt
