/**
 * @file test_convolution.cxx
 * @brief Test application for parallel convolution with multiple execution strategies
 * 
 * This application demonstrates:
 * - Block decomposition of 2D/3D data with variable ghost zones
 * - Parallel convolution using OpenMP, SYCL, or C++17 stdpar
 * - Reading and writing raw float32 files
 * 
 * Usage:
 *   test_convolution <input_file> <output_file> <X> <Y> <Z> <execution_target> 
 *                    [kernel_type] [kernel_size] [block_size_x] [block_size_y] [block_size_z] [ghost_zone_size]
 * 
 * Arguments:
 *   input_file       - Raw float32 input file
 *   output_file      - Raw float32 output file
 *   X, Y, Z          - Grid dimensions
 *   execution_target - cpu, sycl, or stdpar
 *   kernel_type      - identity, gaussian, sobel_x, sobel_y, sobel_z, laplacian, box_blur, sharpen (default: gaussian)
 *   kernel_size      - Kernel dimension (3, 5, 7, etc., default depends on kernel type)
 *   block_size_*     - Block dimensions (default: 64)
 *   ghost_zone_size  - Ghost cell layers (default: auto from kernel_size)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <cstring>
#include <stdexcept>

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"
#include "block_decomposition.h"
#include "convolution_kernel.h"
#include "execution_strategy.h"
#include "execution_strategy_cpu.h"
#include "execution_strategy_sycl.h"
#include "execution_strategy_stdpar.h"

using namespace GInt;

// Factory function implementation (declared in execution_strategy.h)
namespace GInt {
std::unique_ptr<ExecutionStrategy> createExecutionStrategy(ExecutionTarget target) {
    switch (target) {
        case ExecutionTarget::CPU_OPENMP:
            return std::make_unique<ExecutionStrategyCPU>();
        case ExecutionTarget::CPU_STDPAR:
            return std::make_unique<ExecutionStrategyStdpar>();
        case ExecutionTarget::GPU_SYCL:
            return std::make_unique<ExecutionStrategySYCL>();
        default:
            throw std::invalid_argument("Unknown execution target");
    }
}
}

void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <input_file> <output_file> <X> <Y> <Z> <execution_target> "
              << "[kernel_type] [kernel_size] [block_size_x] [block_size_y] [block_size_z] [ghost_zone_size]\n\n"
              << "Arguments:\n"
              << "  input_file       - Raw float32 input file\n"
              << "  output_file      - Raw float32 output file\n"
              << "  X, Y, Z          - Grid dimensions\n"
              << "  execution_target - cpu, sycl, or stdpar\n"
              << "  kernel_type      - identity, gaussian, sobel_x, sobel_y, sobel_z, laplacian, box_blur, sharpen\n"
              << "  kernel_size      - Kernel dimension (3, 5, 7, etc.)\n"
              << "  block_size_*     - Block dimensions (default: 64)\n"
              << "  ghost_zone_size  - Ghost cell layers (default: auto from kernel_size)\n";
}

/**
 * @brief Read raw float32 data from file
 */
bool readRawFile(const std::string& filename, std::vector<float>& data, INDEX_TYPE expected_size) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open input file: " << filename << std::endl;
        return false;
    }
    
    // Get file size
    file.seekg(0, std::ios::end);
    size_t file_size = file.tellg();
    file.seekg(0, std::ios::beg);
    
    size_t expected_bytes = expected_size * sizeof(float);
    if (file_size != expected_bytes) {
        std::cerr << "Error: File size mismatch. Expected " << expected_bytes 
                  << " bytes, got " << file_size << " bytes.\n";
        return false;
    }
    
    data.resize(expected_size);
    file.read(reinterpret_cast<char*>(data.data()), expected_bytes);
    
    if (!file) {
        std::cerr << "Error: Failed to read data from file.\n";
        return false;
    }
    
    return true;
}

/**
 * @brief Write raw float32 data to file
 */
bool writeRawFile(const std::string& filename, const std::vector<float>& data) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error: Could not open output file: " << filename << std::endl;
        return false;
    }
    
    file.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
    
    if (!file) {
        std::cerr << "Error: Failed to write data to file.\n";
        return false;
    }
    
    return true;
}

int main(int argc, char* argv[]) {
    // Parse command line arguments
    if (argc < 7) {
        printUsage(argv[0]);
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    INDEX_TYPE dim_x = std::stoll(argv[3]);
    INDEX_TYPE dim_y = std::stoll(argv[4]);
    INDEX_TYPE dim_z = std::stoll(argv[5]);
    std::string exec_target_str = argv[6];
    
    // Optional arguments
    std::string kernel_type_str = (argc > 7) ? argv[7] : "gaussian";
    int kernel_size = (argc > 8) ? std::stoi(argv[8]) : 0;  // 0 = use default for kernel type
    INDEX_TYPE block_size_x = (argc > 9) ? std::stoll(argv[9]) : 64;
    INDEX_TYPE block_size_y = (argc > 10) ? std::stoll(argv[10]) : 64;
    INDEX_TYPE block_size_z = (argc > 11) ? std::stoll(argv[11]) : 64;
    int ghost_zone_size = (argc > 12) ? std::stoi(argv[12]) : -1;  // -1 = auto
    
    // Determine if 2D or 3D
    bool is_2d = (dim_z == 1);
    
    std::cout << "=== Parallel Convolution Test Application ===\n\n";
    std::cout << "Configuration:\n";
    std::cout << "  Input file:  " << input_file << "\n";
    std::cout << "  Output file: " << output_file << "\n";
    std::cout << "  Dimensions:  " << dim_x << " x " << dim_y << " x " << dim_z;
    if (is_2d) std::cout << " (2D mode)";
    std::cout << "\n";
    
    try {
        // Parse execution target
        ExecutionTarget exec_target = parseExecutionTarget(exec_target_str);
        
        // Parse kernel type
        KernelType kernel_type = parseKernelType(kernel_type_str);
        
        // Get default kernel size if not specified
        if (kernel_size == 0) {
            kernel_size = getDefaultKernelSize(kernel_type);
        }
        
        // Calculate ghost zone size if not specified
        if (ghost_zone_size < 0) {
            ghost_zone_size = getRequiredGhostZoneSize(kernel_size);
        }
        
        // Validate ghost zone size
        int required_ghost = getRequiredGhostZoneSize(kernel_size);
        if (ghost_zone_size < required_ghost) {
            std::cerr << "Warning: Specified ghost zone size (" << ghost_zone_size 
                      << ") is smaller than required (" << required_ghost 
                      << ") for kernel size " << kernel_size << ". Using required size.\n";
            ghost_zone_size = required_ghost;
        }
        
        std::cout << "  Execution:   " << executionTargetName(exec_target) << "\n";
        std::cout << "  Kernel:      " << kernel_type_str << " (" << kernel_size << "x" << kernel_size;
        if (!is_2d) std::cout << "x" << kernel_size;
        std::cout << ")\n";
        std::cout << "  Block size:  " << block_size_x << " x " << block_size_y;
        if (!is_2d) std::cout << " x " << block_size_z;
        std::cout << "\n";
        std::cout << "  Ghost zone:  " << ghost_zone_size << " layers\n\n";
        
        // Create execution strategy
        auto strategy = createExecutionStrategy(exec_target);
        if (!strategy->isAvailable()) {
            std::cerr << "Error: Execution strategy " << strategy->name() << " is not available.\n";
            return 1;
        }
        
        // Calculate total elements
        INDEX_TYPE total_elements = dim_x * dim_y * dim_z;
        
        // Read input data
        std::cout << "Reading input file...\n";
        std::vector<float> input_data;
        if (!readRawFile(input_file, input_data, total_elements)) {
            return 1;
        }
        std::cout << "  Read " << input_data.size() << " elements.\n\n";
        
        // Generate kernel
        std::cout << "Generating convolution kernel...\n";
        std::vector<float> kernel;
        if (is_2d) {
            kernel = generateKernel2D(kernel_type, kernel_size);
        } else {
            kernel = generateKernel3D(kernel_type, kernel_size);
        }
        std::cout << "  Kernel size: " << kernel.size() << " elements.\n\n";
        
        // Prepare output data
        std::vector<float> output_data(total_elements, 0.0f);
        
        // Timing
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (is_2d) {
            // 2D processing
            Vec2l global_dims(dim_x, dim_y);
            Vec2l block_size(block_size_x, block_size_y);
            
            std::cout << "Creating 2D block decomposition...\n";
            BlockDecomposition2D<float> decomp(global_dims, block_size, ghost_zone_size);
            std::cout << "  Blocks: " << decomp.blocksPerAxis()[0] << " x " << decomp.blocksPerAxis()[1] 
                      << " = " << decomp.numBlocks() << " total\n\n";
            
            // Decompose input into blocks
            std::cout << "Decomposing input data into blocks...\n";
            std::vector<BlockData2D<float>> input_blocks;
            decomp.decompose(input_data.data(), input_blocks);
            std::cout << "  Created " << input_blocks.size() << " blocks.\n\n";
            
            // Create output blocks
            std::vector<BlockData2D<float>> output_blocks;
            
            // Execute convolution
            std::cout << "Executing 2D convolution using " << strategy->name() << "...\n";
            auto conv_start = std::chrono::high_resolution_clock::now();
            
            strategy->executeConvolution2D(input_blocks, output_blocks, kernel.data(), kernel_size, ghost_zone_size);
            
            auto conv_end = std::chrono::high_resolution_clock::now();
            auto conv_duration = std::chrono::duration_cast<std::chrono::milliseconds>(conv_end - conv_start);
            std::cout << "  Convolution completed in " << conv_duration.count() << " ms.\n\n";
            
            // Gather results
            std::cout << "Gathering results...\n";
            decomp.gather(output_blocks, output_data.data());
            
        } else {
            // 3D processing
            Vec3l global_dims(dim_x, dim_y, dim_z);
            Vec3l block_size(block_size_x, block_size_y, block_size_z);
            
            std::cout << "Creating 3D block decomposition...\n";
            BlockDecomposition<float> decomp(global_dims, block_size, ghost_zone_size);
            std::cout << "  Blocks: " << decomp.blocksPerAxis()[0] << " x " 
                      << decomp.blocksPerAxis()[1] << " x " << decomp.blocksPerAxis()[2]
                      << " = " << decomp.numBlocks() << " total\n\n";
            
            // Decompose input into blocks
            std::cout << "Decomposing input data into blocks...\n";
            std::vector<BlockData<float>> input_blocks;
            decomp.decompose(input_data.data(), input_blocks);
            std::cout << "  Created " << input_blocks.size() << " blocks.\n\n";
            
            // Create output blocks
            std::vector<BlockData<float>> output_blocks;
            
            // Execute convolution
            std::cout << "Executing 3D convolution using " << strategy->name() << "...\n";
            auto conv_start = std::chrono::high_resolution_clock::now();
            
            strategy->executeConvolution3D(input_blocks, output_blocks, kernel.data(), kernel_size, ghost_zone_size);
            
            auto conv_end = std::chrono::high_resolution_clock::now();
            auto conv_duration = std::chrono::duration_cast<std::chrono::milliseconds>(conv_end - conv_start);
            std::cout << "  Convolution completed in " << conv_duration.count() << " ms.\n\n";
            
            // Gather results
            std::cout << "Gathering results...\n";
            decomp.gather(output_blocks, output_data.data());
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        
        // Write output
        std::cout << "Writing output file...\n";
        if (!writeRawFile(output_file, output_data)) {
            return 1;
        }
        std::cout << "  Wrote " << output_data.size() << " elements.\n\n";
        
        std::cout << "=== Processing Complete ===\n";
        std::cout << "Total time: " << total_duration.count() << " ms.\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
