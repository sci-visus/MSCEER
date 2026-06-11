#ifndef EXECUTION_STRATEGY_STDPAR_H
#define EXECUTION_STRATEGY_STDPAR_H

#include "execution_strategy.h"

#ifdef GINT_USE_STDPAR
#include <execution>
#include <algorithm>
#include <numeric>
#endif

namespace GInt {

/**
 * @brief C++17 parallel algorithms (stdpar) execution strategy
 * 
 * Uses std::execution::par_unseq for parallel execution.
 * Requires GINT_USE_STDPAR to be defined and a compiler with stdpar support.
 * Works with NVIDIA HPC SDK, Intel oneAPI, or GCC with TBB.
 */
class ExecutionStrategyStdpar : public ExecutionStrategy {
#ifdef GINT_USE_STDPAR

public:
    ExecutionStrategyStdpar() = default;
    ~ExecutionStrategyStdpar() override = default;
    
    void executeConvolution3D(
        const std::vector<BlockData<float>>& input_blocks,
        std::vector<BlockData<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) override {
        // Ensure output blocks have same structure
        if (output_blocks.size() != input_blocks.size()) {
            output_blocks.clear();
            output_blocks.reserve(input_blocks.size());
            for (const auto& input : input_blocks) {
                output_blocks.emplace_back(input.interior_size, input.global_start,
                                          input.block_index, input.ghost_zone_size);
            }
        }
        
        const size_t num_blocks = input_blocks.size();
        std::vector<size_t> block_indices(num_blocks);
        std::iota(block_indices.begin(), block_indices.end(), 0);
        
        // Process blocks in parallel
        std::for_each(std::execution::par_unseq, block_indices.begin(), block_indices.end(),
            [&](size_t b) {
                const auto& input = input_blocks[b];
                auto& output = output_blocks[b];
                
                const int half = kernel_size / 2;
                INDEX_TYPE total_x = input.total_size[0];
                INDEX_TYPE total_y = input.total_size[1];
                int gz = input.ghost_zone_size;
                
                // Create index vector for interior elements
                INDEX_TYPE num_interior = input.interior_size[0] * input.interior_size[1] * input.interior_size[2];
                std::vector<INDEX_TYPE> indices(num_interior);
                std::iota(indices.begin(), indices.end(), 0);
                
                // Process each interior element in parallel
                std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                    [&](INDEX_TYPE idx) {
                        INDEX_TYPE ix = idx % input.interior_size[0];
                        INDEX_TYPE iy = (idx / input.interior_size[0]) % input.interior_size[1];
                        INDEX_TYPE iz = idx / (input.interior_size[0] * input.interior_size[1]);
                        
                        INDEX_TYPE lx = ix + gz;
                        INDEX_TYPE ly = iy + gz;
                        INDEX_TYPE lz = iz + gz;
                        
                        float sum = 0.0f;
                        for (int kz = -half; kz <= half; kz++) {
                            for (int ky = -half; ky <= half; ky++) {
                                for (int kx = -half; kx <= half; kx++) {
                                    INDEX_TYPE in_idx = (lx + kx) + (ly + ky) * total_x + 
                                                        (lz + kz) * total_x * total_y;
                                    INDEX_TYPE k_idx = (kz + half) * kernel_size * kernel_size +
                                                       (ky + half) * kernel_size + (kx + half);
                                    sum += kernel[k_idx] * input.data[in_idx];
                                }
                            }
                        }
                        
                        INDEX_TYPE out_idx = (ix + gz) + (iy + gz) * total_x + 
                                            (iz + gz) * total_x * total_y;
                        output.data[out_idx] = sum;
                    });
            });
    }
    
    void executeConvolution2D(
        const std::vector<BlockData2D<float>>& input_blocks,
        std::vector<BlockData2D<float>>& output_blocks,
        const float* kernel,
        int kernel_size,
        int ghost_zone_size
    ) override {
        if (output_blocks.size() != input_blocks.size()) {
            output_blocks.clear();
            output_blocks.reserve(input_blocks.size());
            for (const auto& input : input_blocks) {
                output_blocks.emplace_back(input.interior_size, input.global_start,
                                          input.block_index, input.ghost_zone_size);
            }
        }
        
        const size_t num_blocks = input_blocks.size();
        std::vector<size_t> block_indices(num_blocks);
        std::iota(block_indices.begin(), block_indices.end(), 0);
        
        std::for_each(std::execution::par_unseq, block_indices.begin(), block_indices.end(),
            [&](size_t b) {
                const auto& input = input_blocks[b];
                auto& output = output_blocks[b];
                
                const int half = kernel_size / 2;
                INDEX_TYPE total_x = input.total_size[0];
                int gz = input.ghost_zone_size;
                
                INDEX_TYPE num_interior = input.interior_size[0] * input.interior_size[1];
                std::vector<INDEX_TYPE> indices(num_interior);
                std::iota(indices.begin(), indices.end(), 0);
                
                std::for_each(std::execution::par_unseq, indices.begin(), indices.end(),
                    [&](INDEX_TYPE idx) {
                        INDEX_TYPE ix = idx % input.interior_size[0];
                        INDEX_TYPE iy = idx / input.interior_size[0];
                        
                        INDEX_TYPE lx = ix + gz;
                        INDEX_TYPE ly = iy + gz;
                        
                        float sum = 0.0f;
                        for (int ky = -half; ky <= half; ky++) {
                            for (int kx = -half; kx <= half; kx++) {
                                INDEX_TYPE in_idx = (lx + kx) + (ly + ky) * total_x;
                                INDEX_TYPE k_idx = (ky + half) * kernel_size + (kx + half);
                                sum += kernel[k_idx] * input.data[in_idx];
                            }
                        }
                        
                        INDEX_TYPE out_idx = (ix + gz) + (iy + gz) * total_x;
                        output.data[out_idx] = sum;
                    });
            });
    }
    
    std::string name() const override { return "CPU_Stdpar"; }
    ExecutionTarget target() const override { return ExecutionTarget::CPU_STDPAR; }
    bool isAvailable() const override { return true; }
    
#else // !GINT_USE_STDPAR

public:
    ExecutionStrategyStdpar() {
        throw std::runtime_error("C++17 stdpar support not compiled. Define GINT_USE_STDPAR and use a compatible compiler.");
    }
    ~ExecutionStrategyStdpar() override = default;
    
    void executeConvolution3D(
        const std::vector<BlockData<float>>&,
        std::vector<BlockData<float>>&,
        const float*, int, int
    ) override {
        throw std::runtime_error("Stdpar not available");
    }
    
    void executeConvolution2D(
        const std::vector<BlockData2D<float>>&,
        std::vector<BlockData2D<float>>&,
        const float*, int, int
    ) override {
        throw std::runtime_error("Stdpar not available");
    }
    
    std::string name() const override { return "CPU_Stdpar (not available)"; }
    ExecutionTarget target() const override { return ExecutionTarget::CPU_STDPAR; }
    bool isAvailable() const override { return false; }
    
#endif // GINT_USE_STDPAR
};

} // namespace GInt

#endif // EXECUTION_STRATEGY_STDPAR_H
