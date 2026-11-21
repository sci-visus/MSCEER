#ifndef EXECUTION_STRATEGY_CPU_H
#define EXECUTION_STRATEGY_CPU_H

#include "execution_strategy.h"
#include "convolution_kernel.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GInt {

/**
 * @brief OpenMP-based CPU execution strategy
 */
class ExecutionStrategyCPU : public ExecutionStrategy {
public:
    ExecutionStrategyCPU() = default;
    ~ExecutionStrategyCPU() override = default;
    
    void executeConvolution3D(
        const std::vector<BlockData<float>>& input_blocks,
        std::vector<BlockData<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) override {
        (void)ghost_zone_size; // Used for validation only
        
        // Ensure output blocks have same structure as input blocks
        if (output_blocks.size() != input_blocks.size()) {
            output_blocks.clear();
            output_blocks.reserve(input_blocks.size());
            for (const auto& input : input_blocks) {
                output_blocks.emplace_back(input.interior_size, input.global_start, 
                                          input.block_index, input.ghost_zone_size);
            }
        }
        
        const INDEX_TYPE num_blocks = static_cast<INDEX_TYPE>(input_blocks.size());
        
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (INDEX_TYPE i = 0; i < num_blocks; i++) {
            convolveBlock3D(input_blocks[i], output_blocks[i], kernel, kernel_size);
        }
    }
    
    void executeConvolution2D(
        const std::vector<BlockData2D<float>>& input_blocks,
        std::vector<BlockData2D<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) override {
        (void)ghost_zone_size;
        
        // Ensure output blocks have same structure as input blocks
        if (output_blocks.size() != input_blocks.size()) {
            output_blocks.clear();
            output_blocks.reserve(input_blocks.size());
            for (const auto& input : input_blocks) {
                output_blocks.emplace_back(input.interior_size, input.global_start,
                                          input.block_index, input.ghost_zone_size);
            }
        }
        
        const INDEX_TYPE num_blocks = static_cast<INDEX_TYPE>(input_blocks.size());
        
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
#endif
        for (INDEX_TYPE i = 0; i < num_blocks; i++) {
            convolveBlock2D(input_blocks[i], output_blocks[i], kernel, kernel_size);
        }
    }
    
    std::string name() const override {
        return "CPU_OpenMP";
    }
    
    ExecutionTarget target() const override {
        return ExecutionTarget::CPU_OPENMP;
    }
    
    bool isAvailable() const override {
        // CPU execution is always available
        return true;
    }
};

} // namespace GInt

#endif // EXECUTION_STRATEGY_CPU_H
