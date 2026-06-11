#ifndef CONVOLUTION_KERNEL_H
#define CONVOLUTION_KERNEL_H

#include <vector>
#include <cmath>
#include <stdexcept>

#include "gi_basic_types.h"
#include "block_decomposition.h"

namespace GInt {

/**
 * @brief Calculate required ghost zone size for a given kernel size
 * @param kernel_size The size of the kernel (must be odd)
 * @return The number of ghost cell layers needed
 */
inline int getRequiredGhostZoneSize(int kernel_size) {
    if (kernel_size < 1 || kernel_size % 2 == 0) {
        throw std::invalid_argument("Kernel size must be a positive odd number");
    }
    return (kernel_size - 1) / 2;
}

/**
 * @brief Validate that block's ghost zone is sufficient for kernel
 */
inline void validateGhostZone(int block_ghost_size, int kernel_size) {
    int required = getRequiredGhostZoneSize(kernel_size);
    if (block_ghost_size < required) {
        throw std::invalid_argument("Block ghost zone size insufficient for kernel");
    }
}

/**
 * @brief Enumeration of built-in kernel types
 */
enum class KernelType {
    IDENTITY,
    GAUSSIAN,
    SOBEL_X,
    SOBEL_Y,
    SOBEL_Z,
    LAPLACIAN,
    BOX_BLUR,
    SHARPEN
};

/**
 * @brief Get default kernel size for a kernel type
 */
inline int getDefaultKernelSize(KernelType type) {
    switch (type) {
        case KernelType::IDENTITY:
        case KernelType::SOBEL_X:
        case KernelType::SOBEL_Y:
        case KernelType::SOBEL_Z:
        case KernelType::LAPLACIAN:
        case KernelType::BOX_BLUR:
        case KernelType::SHARPEN:
            return 3;
        case KernelType::GAUSSIAN:
            return 5;
        default:
            return 3;
    }
}

/**
 * @brief Generate a 2D convolution kernel
 * @param type Kernel type
 * @param size Kernel size (must be odd)
 * @return Flattened kernel array (row-major)
 */
inline std::vector<float> generateKernel2D(KernelType type, int size = 0) {
    if (size == 0) size = getDefaultKernelSize(type);
    
    std::vector<float> kernel(size * size, 0.0f);
    int half = size / 2;
    
    switch (type) {
        case KernelType::IDENTITY:
            kernel[half * size + half] = 1.0f;
            break;
            
        case KernelType::BOX_BLUR:
            {
                float val = 1.0f / (size * size);
                for (int i = 0; i < size * size; i++) {
                    kernel[i] = val;
                }
            }
            break;
            
        case KernelType::GAUSSIAN:
            {
                float sigma = size / 6.0f;
                float sum = 0.0f;
                for (int y = -half; y <= half; y++) {
                    for (int x = -half; x <= half; x++) {
                        float val = std::exp(-(x*x + y*y) / (2.0f * sigma * sigma));
                        kernel[(y + half) * size + (x + half)] = val;
                        sum += val;
                    }
                }
                // Normalize
                for (int i = 0; i < size * size; i++) {
                    kernel[i] /= sum;
                }
            }
            break;
            
        case KernelType::SOBEL_X:
            if (size != 3) throw std::invalid_argument("Sobel kernel must be 3x3");
            kernel = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
            break;
            
        case KernelType::SOBEL_Y:
            if (size != 3) throw std::invalid_argument("Sobel kernel must be 3x3");
            kernel = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
            break;
            
        case KernelType::LAPLACIAN:
            if (size != 3) throw std::invalid_argument("Laplacian kernel must be 3x3");
            kernel = {0, 1, 0, 1, -4, 1, 0, 1, 0};
            break;
            
        case KernelType::SHARPEN:
            if (size != 3) throw std::invalid_argument("Sharpen kernel must be 3x3");
            kernel = {0, -1, 0, -1, 5, -1, 0, -1, 0};
            break;
            
        default:
            kernel[half * size + half] = 1.0f;
    }
    
    return kernel;
}

/**
 * @brief Generate a 3D convolution kernel
 */
inline std::vector<float> generateKernel3D(KernelType type, int size = 0) {
    if (size == 0) size = getDefaultKernelSize(type);
    
    std::vector<float> kernel(size * size * size, 0.0f);
    int half = size / 2;
    
    switch (type) {
        case KernelType::IDENTITY:
            kernel[half * size * size + half * size + half] = 1.0f;
            break;
            
        case KernelType::BOX_BLUR:
            {
                float val = 1.0f / (size * size * size);
                for (int i = 0; i < size * size * size; i++) {
                    kernel[i] = val;
                }
            }
            break;
            
        case KernelType::GAUSSIAN:
            {
                float sigma = size / 6.0f;
                float sum = 0.0f;
                for (int z = -half; z <= half; z++) {
                    for (int y = -half; y <= half; y++) {
                        for (int x = -half; x <= half; x++) {
                            float val = std::exp(-(x*x + y*y + z*z) / (2.0f * sigma * sigma));
                            kernel[(z + half) * size * size + (y + half) * size + (x + half)] = val;
                            sum += val;
                        }
                    }
                }
                // Normalize
                for (int i = 0; i < size * size * size; i++) {
                    kernel[i] /= sum;
                }
            }
            break;
            
        case KernelType::LAPLACIAN:
            if (size != 3) throw std::invalid_argument("3D Laplacian kernel must be 3x3x3");
            // 6-connected Laplacian
            for (int i = 0; i < 27; i++) kernel[i] = 0.0f;
            kernel[4] = 1.0f;   // (0,1,1) -> -z face
            kernel[10] = 1.0f;  // (1,0,1) -> -y face  
            kernel[12] = 1.0f;  // (1,1,0) -> -x face
            kernel[13] = -6.0f; // center
            kernel[14] = 1.0f;  // (1,1,2) -> +x face
            kernel[16] = 1.0f;  // (1,2,1) -> +y face
            kernel[22] = 1.0f;  // (2,1,1) -> +z face
            break;
            
        case KernelType::SOBEL_X:
        case KernelType::SOBEL_Y:
        case KernelType::SOBEL_Z:
            // 3D Sobel is more complex, use simple gradient for now
            if (size != 3) throw std::invalid_argument("3D Sobel kernel must be 3x3x3");
            for (int i = 0; i < 27; i++) kernel[i] = 0.0f;
            if (type == KernelType::SOBEL_X) {
                kernel[12] = -1.0f; kernel[14] = 1.0f; // x gradient at center
            } else if (type == KernelType::SOBEL_Y) {
                kernel[10] = -1.0f; kernel[16] = 1.0f; // y gradient at center
            } else {
                kernel[4] = -1.0f; kernel[22] = 1.0f;  // z gradient at center
            }
            break;
            
        default:
            kernel[half * size * size + half * size + half] = 1.0f;
    }
    
    return kernel;
}

/**
 * @brief Apply 2D convolution to a single block
 * @param input Input block (with ghost cells)
 * @param output Output block (will write to interior only)
 * @param kernel Flattened kernel array (row-major, size kernel_size * kernel_size)
 * @param kernel_size Kernel dimension (must be odd)
 */
template<typename T>
void convolveBlock2D(const BlockData2D<T>& input, BlockData2D<T>& output,
                     const float* kernel, int kernel_size) {
    validateGhostZone(input.ghost_zone_size, kernel_size);
    
    int half = kernel_size / 2;
    
    // Only write to interior region
    for (INDEX_TYPE iy = 0; iy < input.interior_size[1]; iy++) {
        for (INDEX_TYPE ix = 0; ix < input.interior_size[0]; ix++) {
            // Local coordinates in the block (including ghost offset)
            INDEX_TYPE lx = ix + input.ghost_zone_size;
            INDEX_TYPE ly = iy + input.ghost_zone_size;
            
            float sum = 0.0f;
            for (int ky = -half; ky <= half; ky++) {
                for (int kx = -half; kx <= half; kx++) {
                    float kval = kernel[(ky + half) * kernel_size + (kx + half)];
                    sum += kval * static_cast<float>(input.at(lx + kx, ly + ky));
                }
            }
            output.interiorAt(ix, iy) = static_cast<T>(sum);
        }
    }
}

/**
 * @brief Apply 3D convolution to a single block
 */
template<typename T>
void convolveBlock3D(const BlockData<T>& input, BlockData<T>& output,
                     const float* kernel, int kernel_size) {
    validateGhostZone(input.ghost_zone_size, kernel_size);
    
    int half = kernel_size / 2;
    
    // Only write to interior region
    for (INDEX_TYPE iz = 0; iz < input.interior_size[2]; iz++) {
        for (INDEX_TYPE iy = 0; iy < input.interior_size[1]; iy++) {
            for (INDEX_TYPE ix = 0; ix < input.interior_size[0]; ix++) {
                // Local coordinates in the block (including ghost offset)
                INDEX_TYPE lx = ix + input.ghost_zone_size;
                INDEX_TYPE ly = iy + input.ghost_zone_size;
                INDEX_TYPE lz = iz + input.ghost_zone_size;
                
                float sum = 0.0f;
                for (int kz = -half; kz <= half; kz++) {
                    for (int ky = -half; ky <= half; ky++) {
                        for (int kx = -half; kx <= half; kx++) {
                            float kval = kernel[(kz + half) * kernel_size * kernel_size + 
                                                (ky + half) * kernel_size + (kx + half)];
                            sum += kval * static_cast<float>(input.at(lx + kx, ly + ky, lz + kz));
                        }
                    }
                }
                output.interiorAt(ix, iy, iz) = static_cast<T>(sum);
            }
        }
    }
}

/**
 * @brief Parse kernel type from string
 */
inline KernelType parseKernelType(const std::string& name) {
    if (name == "identity") return KernelType::IDENTITY;
    if (name == "gaussian") return KernelType::GAUSSIAN;
    if (name == "sobel_x") return KernelType::SOBEL_X;
    if (name == "sobel_y") return KernelType::SOBEL_Y;
    if (name == "sobel_z") return KernelType::SOBEL_Z;
    if (name == "laplacian") return KernelType::LAPLACIAN;
    if (name == "box_blur" || name == "box") return KernelType::BOX_BLUR;
    if (name == "sharpen") return KernelType::SHARPEN;
    throw std::invalid_argument("Unknown kernel type: " + name);
}

} // namespace GInt

#endif // CONVOLUTION_KERNEL_H
