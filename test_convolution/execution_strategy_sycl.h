#ifndef EXECUTION_STRATEGY_SYCL_H
#define EXECUTION_STRATEGY_SYCL_H

#include "execution_strategy.h"

#ifdef GINT_USE_SYCL
// Support both old SYCL 1.2.1 format and modern DPC++ format
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#elif __has_include(<CL/sycl.hpp>)
#include <CL/sycl.hpp>
#else
#error "SYCL headers not found. Please install Intel oneAPI DPC++ or another SYCL implementation."
#endif
#endif

namespace GInt {

/**
 * @brief SYCL-based GPU execution strategy
 * 
 * This implementation uses SYCL for GPU offloading.
 * Requires GINT_USE_SYCL to be defined at compile time.
 */
class ExecutionStrategySYCL : public ExecutionStrategy {
#ifdef GINT_USE_SYCL
private:
    sycl::queue m_queue;
    
public:
    ExecutionStrategySYCL() {
        try {
            // Try to use GPU selector (modern SYCL 2020)
            #if defined(SYCL_LANGUAGE_VERSION) && SYCL_LANGUAGE_VERSION >= 202001
                m_queue = sycl::queue(sycl::gpu_selector_v);
            #else
                // Fallback for older SYCL implementations
                sycl::gpu_selector selector;
                m_queue = sycl::queue(selector);
            #endif
            
            // Print device info
            auto device = m_queue.get_device();
            std::cout << "SYCL device: " << device.get_info<sycl::info::device::name>() << std::endl;
            std::cout << "  Vendor: " << device.get_info<sycl::info::device::vendor>() << std::endl;
            std::cout << "  Type: " << (device.is_gpu() ? "GPU" : "Other") << std::endl;
        } catch (const sycl::exception& e) {
            std::cerr << "SYCL GPU device selection failed: " << e.what() << std::endl;
            std::cerr << "Falling back to default device selector..." << std::endl;
            m_queue = sycl::queue();  // Use default selector
            auto device = m_queue.get_device();
            std::cout << "SYCL device (default): " << device.get_info<sycl::info::device::name>() << std::endl;
        }
    }
    
    ~ExecutionStrategySYCL() override = default;
    
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
        
        const int half = kernel_size / 2;
        
        // Create kernel buffer
        sycl::buffer<float, 1> kernel_buf(kernel, sycl::range<1>(kernel_size * kernel_size * kernel_size));
        
        // Process each block
        for (size_t b = 0; b < input_blocks.size(); b++) {
            const auto& input = input_blocks[b];
            auto& output = output_blocks[b];
            
            // Create buffers for this block
            sycl::buffer<float, 1> input_buf(input.data.data(), 
                                            sycl::range<1>(input.data.size()));
            sycl::buffer<float, 1> output_buf(output.data.data(), 
                                             sycl::range<1>(output.data.size()));
            
            INDEX_TYPE interior_x = input.interior_size[0];
            INDEX_TYPE interior_y = input.interior_size[1];
            INDEX_TYPE interior_z = input.interior_size[2];
            INDEX_TYPE total_x = input.total_size[0];
            INDEX_TYPE total_y = input.total_size[1];
            int gz = input.ghost_zone_size;
            
            m_queue.submit([&](sycl::handler& h) {
                auto in = input_buf.get_access<sycl::access::mode::read>(h);
                auto out = output_buf.get_access<sycl::access::mode::write>(h);
                auto kern = kernel_buf.get_access<sycl::access::mode::read>(h);
                
                h.parallel_for(sycl::range<3>(interior_z, interior_y, interior_x),
                    [=](sycl::id<3> idx) {
                        INDEX_TYPE iz = idx[0];
                        INDEX_TYPE iy = idx[1];
                        INDEX_TYPE ix = idx[2];
                        
                        // Local coordinates with ghost offset
                        INDEX_TYPE lx = ix + gz;
                        INDEX_TYPE ly = iy + gz;
                        INDEX_TYPE lz = iz + gz;
                        
                        float sum = 0.0f;
                        for (int kz = -half; kz <= half; kz++) {
                            for (int ky = -half; ky <= half; ky++) {
                                for (int kx = -half; kx <= half; kx++) {
                                    INDEX_TYPE idx_in = (lx + kx) + (ly + ky) * total_x + 
                                                        (lz + kz) * total_x * total_y;
                                    INDEX_TYPE idx_k = (kz + half) * kernel_size * kernel_size +
                                                       (ky + half) * kernel_size + (kx + half);
                                    sum += kern[idx_k] * in[idx_in];
                                }
                            }
                        }
                        
                        // Write to output at interior position (with ghost offset)
                        INDEX_TYPE out_idx = (ix + gz) + (iy + gz) * total_x + 
                                            (iz + gz) * total_x * total_y;
                        out[out_idx] = sum;
                    });
            }).wait();
        }
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
        
        const int half = kernel_size / 2;
        sycl::buffer<float, 1> kernel_buf(kernel, sycl::range<1>(kernel_size * kernel_size));
        
        for (size_t b = 0; b < input_blocks.size(); b++) {
            const auto& input = input_blocks[b];
            auto& output = output_blocks[b];
            
            sycl::buffer<float, 1> input_buf(input.data.data(), sycl::range<1>(input.data.size()));
            sycl::buffer<float, 1> output_buf(output.data.data(), sycl::range<1>(output.data.size()));
            
            INDEX_TYPE interior_x = input.interior_size[0];
            INDEX_TYPE interior_y = input.interior_size[1];
            INDEX_TYPE total_x = input.total_size[0];
            int gz = input.ghost_zone_size;
            
            m_queue.submit([&](sycl::handler& h) {
                auto in = input_buf.get_access<sycl::access::mode::read>(h);
                auto out = output_buf.get_access<sycl::access::mode::write>(h);
                auto kern = kernel_buf.get_access<sycl::access::mode::read>(h);
                
                h.parallel_for(sycl::range<2>(interior_y, interior_x),
                    [=](sycl::id<2> idx) {
                        INDEX_TYPE iy = idx[0];
                        INDEX_TYPE ix = idx[1];
                        INDEX_TYPE lx = ix + gz;
                        INDEX_TYPE ly = iy + gz;
                        
                        float sum = 0.0f;
                        for (int ky = -half; ky <= half; ky++) {
                            for (int kx = -half; kx <= half; kx++) {
                                INDEX_TYPE idx_in = (lx + kx) + (ly + ky) * total_x;
                                INDEX_TYPE idx_k = (ky + half) * kernel_size + (kx + half);
                                sum += kern[idx_k] * in[idx_in];
                            }
                        }
                        
                        INDEX_TYPE out_idx = (ix + gz) + (iy + gz) * total_x;
                        out[out_idx] = sum;
                    });
            }).wait();
        }
    }
    
    std::string name() const override { return "GPU_SYCL"; }
    ExecutionTarget target() const override { return ExecutionTarget::GPU_SYCL; }
    bool isAvailable() const override { return true; }
    
#else // !GINT_USE_SYCL

public:
    ExecutionStrategySYCL() {
        throw std::runtime_error("SYCL support not compiled. Define GINT_USE_SYCL and link SYCL libraries.");
    }
    ~ExecutionStrategySYCL() override = default;
    
    void executeConvolution3D(
        const std::vector<BlockData<float>>&,
        std::vector<BlockData<float>>&,
        const float*, int, int
    ) override {
        throw std::runtime_error("SYCL not available");
    }
    
    void executeConvolution2D(
        const std::vector<BlockData2D<float>>&,
        std::vector<BlockData2D<float>>&,
        const float*, int, int
    ) override {
        throw std::runtime_error("SYCL not available");
    }
    
    std::string name() const override { return "GPU_SYCL (not available)"; }
    ExecutionTarget target() const override { return ExecutionTarget::GPU_SYCL; }
    bool isAvailable() const override { return false; }
    
#endif // GINT_USE_SYCL
};

} // namespace GInt

#endif // EXECUTION_STRATEGY_SYCL_H
