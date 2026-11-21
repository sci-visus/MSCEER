#ifndef EXECUTION_STRATEGY_H
#define EXECUTION_STRATEGY_H

#include <vector>
#include <string>
#include <memory>

#include "gi_basic_types.h"
#include "block_decomposition.h"
#include "convolution_kernel.h"

namespace GInt {

/**
 * @brief Enumeration of execution targets
 * 
 * Note on stdpar: C++17 parallel algorithms run on CPU with MSVC/GCC.
 * Only NVIDIA HPC SDK (nvc++ -stdpar=gpu) offloads to GPU.
 */
enum class ExecutionTarget {
    CPU_OPENMP,     // OpenMP: CPU threads
    CPU_STDPAR,     // C++17 parallel algorithms: CPU threads (MSVC/GCC) or GPU (NVIDIA HPC SDK)
    GPU_SYCL        // SYCL: GPU offloading (Intel oneAPI, etc.)
};

/**
 * @brief Parse execution target from string
 */
inline ExecutionTarget parseExecutionTarget(const std::string& name) {
    if (name == "cpu" || name == "openmp" || name == "CPU_OPENMP") return ExecutionTarget::CPU_OPENMP;
    if (name == "stdpar" || name == "STDPAR" || name == "CPU_STDPAR") return ExecutionTarget::CPU_STDPAR;
    if (name == "sycl" || name == "SYCL" || name == "GPU_SYCL") return ExecutionTarget::GPU_SYCL;
    throw std::invalid_argument("Unknown execution target: " + name);
}

/**
 * @brief Get string name of execution target
 */
inline std::string executionTargetName(ExecutionTarget target) {
    switch (target) {
        case ExecutionTarget::CPU_OPENMP: return "CPU_OPENMP";
        case ExecutionTarget::CPU_STDPAR: return "CPU_STDPAR";
        case ExecutionTarget::GPU_SYCL: return "GPU_SYCL";
        default: return "UNKNOWN";
    }
}

/**
 * @brief Abstract base class for execution strategies
 */
class ExecutionStrategy {
public:
    virtual ~ExecutionStrategy() = default;
    
    /**
     * @brief Execute 3D convolution on all blocks in parallel
     * @param input_blocks Input blocks (with ghost cells populated)
     * @param output_blocks Output blocks (will be written to interior only)
     * @param kernel Flattened kernel array
     * @param kernel_size Kernel dimension
     * @param ghost_zone_size Ghost zone size used
     */
    virtual void executeConvolution3D(
        const std::vector<BlockData<float>>& input_blocks,
        std::vector<BlockData<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) = 0;
    
    /**
     * @brief Execute 2D convolution on all blocks in parallel
     */
    virtual void executeConvolution2D(
        const std::vector<BlockData2D<float>>& input_blocks,
        std::vector<BlockData2D<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) = 0;
    
    /**
     * @brief Get the name of this execution strategy
     */
    virtual std::string name() const = 0;
    
    /**
     * @brief Get the execution target type
     */
    virtual ExecutionTarget target() const = 0;
    
    /**
     * @brief Check if this strategy is available on the current system
     */
    virtual bool isAvailable() const = 0;
};

/**
 * @brief Factory function to create execution strategies
 */
std::unique_ptr<ExecutionStrategy> createExecutionStrategy(ExecutionTarget target);

} // namespace GInt

#endif // EXECUTION_STRATEGY_H
