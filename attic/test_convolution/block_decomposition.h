#ifndef BLOCK_DECOMPOSITION_H
#define BLOCK_DECOMPOSITION_H

#include <vector>
#include <cstring>
#include <algorithm>
#include <stdexcept>

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"

namespace GInt {

/**
 * @brief Data for a single block in the decomposition
 * 
 * Contains the block-local data array with ghost cells, and metadata
 * about the block's position in the global grid.
 */
template<typename T>
struct BlockData {
    // Grid for the local block (includes ghost cells)
    RegularGrid3D local_grid;
    
    // Block-local data array (includes ghost cells)
    std::vector<T> data;
    
    // Starting global coordinates of the interior region (excludes ghost cells)
    Vec3l global_start;
    
    // Block index in the decomposition (i, j, k)
    Vec3l block_index;
    
    // Ghost zone size (number of layers on each side)
    int ghost_zone_size;
    
    // Interior size (without ghost cells)
    Vec3l interior_size;
    
    // Total size including ghost cells
    Vec3l total_size;
    
    BlockData() : local_grid(Vec3l(1,1,1), Vec3b(false,false,false)), ghost_zone_size(0) {}
    
    BlockData(Vec3l interior, Vec3l global_start_, Vec3l block_idx, int ghost_size)
        : interior_size(interior),
          global_start(global_start_),
          block_index(block_idx),
          ghost_zone_size(ghost_size),
          total_size(interior[0] + 2*ghost_size, interior[1] + 2*ghost_size, interior[2] + 2*ghost_size),
          local_grid(total_size, Vec3b(false, false, false))
    {
        data.resize(total_size[0] * total_size[1] * total_size[2], T(0));
    }
    
    // Get local index (with ghost cells) from local coordinates
    INDEX_TYPE localIndex(INDEX_TYPE x, INDEX_TYPE y, INDEX_TYPE z) const {
        return x + y * total_size[0] + z * total_size[0] * total_size[1];
    }
    
    // Get local index from interior coordinates (adds ghost offset)
    INDEX_TYPE interiorToLocalIndex(INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz) const {
        return localIndex(ix + ghost_zone_size, iy + ghost_zone_size, iz + ghost_zone_size);
    }
    
    // Access data at local coordinates (including ghost zone)
    T& at(INDEX_TYPE x, INDEX_TYPE y, INDEX_TYPE z) {
        return data[localIndex(x, y, z)];
    }
    
    const T& at(INDEX_TYPE x, INDEX_TYPE y, INDEX_TYPE z) const {
        return data[localIndex(x, y, z)];
    }
    
    // Access data at interior coordinates (automatically offsets by ghost zone)
    T& interiorAt(INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz) {
        return data[interiorToLocalIndex(ix, iy, iz)];
    }
    
    const T& interiorAt(INDEX_TYPE ix, INDEX_TYPE iy, INDEX_TYPE iz) const {
        return data[interiorToLocalIndex(ix, iy, iz)];
    }
    
    // Get pointer to raw data
    T* dataPtr() { return data.data(); }
    const T* dataPtr() const { return data.data(); }
};

/**
 * @brief Decomposes a 3D grid into blocks with configurable ghost zones
 * 
 * Supports variable ghost zone sizes to accommodate different kernel sizes.
 * Ghost zone size should be at least (kernel_size - 1) / 2 for convolution.
 */
template<typename T>
class BlockDecomposition {
private:
    RegularGrid3D m_global_grid;
    Vec3l m_block_size;          // Interior block size (without ghost cells)
    Vec3l m_blocks_per_axis;
    int m_ghost_zone_size;       // Number of ghost cell layers (same for all dimensions)
    INDEX_TYPE m_num_blocks;
    
public:
    /**
     * @brief Construct block decomposition
     * @param global_dims Global grid dimensions (X, Y, Z)
     * @param block_size Interior block size (without ghost cells)
     * @param ghost_zone_size Number of ghost cell layers on each side
     */
    BlockDecomposition(Vec3l global_dims, Vec3l block_size, int ghost_zone_size)
        : m_global_grid(global_dims, Vec3b(false, false, false)),
          m_block_size(block_size),
          m_ghost_zone_size(ghost_zone_size)
    {
        // Calculate number of blocks per axis
        for (int i = 0; i < 3; i++) {
            m_blocks_per_axis[i] = (global_dims[i] + block_size[i] - 1) / block_size[i];
        }
        m_num_blocks = m_blocks_per_axis[0] * m_blocks_per_axis[1] * m_blocks_per_axis[2];
        
        // Validate ghost zone size
        if (ghost_zone_size < 0) {
            throw std::invalid_argument("Ghost zone size must be non-negative");
        }
    }
    
    /**
     * @brief Construct block decomposition with auto-calculated ghost zone from kernel size
     * @param global_dims Global grid dimensions
     * @param block_size Interior block size
     * @param kernel_size Kernel dimension (ghost_zone = (kernel_size - 1) / 2)
     * @param from_kernel_size Tag to distinguish from other constructor
     */
    struct FromKernelSize {};
    BlockDecomposition(Vec3l global_dims, Vec3l block_size, int kernel_size, FromKernelSize)
        : BlockDecomposition(global_dims, block_size, (kernel_size - 1) / 2)
    {}
    
    // Accessors
    const RegularGrid3D& globalGrid() const { return m_global_grid; }
    Vec3l blocksPerAxis() const { return m_blocks_per_axis; }
    INDEX_TYPE numBlocks() const { return m_num_blocks; }
    int ghostZoneSize() const { return m_ghost_zone_size; }
    Vec3l blockSize() const { return m_block_size; }
    
    /**
     * @brief Get the interior size for a specific block
     * The last block in each dimension may be smaller
     */
    Vec3l getBlockInteriorSize(Vec3l block_idx) const {
        Vec3l interior;
        Vec3l global_dims = m_global_grid.XYZ();
        
        for (int i = 0; i < 3; i++) {
            INDEX_TYPE start = block_idx[i] * m_block_size[i];
            INDEX_TYPE end = std::min(start + m_block_size[i], global_dims[i]);
            interior[i] = end - start;
        }
        return interior;
    }
    
    /**
     * @brief Get the global starting coordinate for a block's interior
     */
    Vec3l getBlockGlobalStart(Vec3l block_idx) const {
        return Vec3l(
            block_idx[0] * m_block_size[0],
            block_idx[1] * m_block_size[1],
            block_idx[2] * m_block_size[2]
        );
    }
    
    /**
     * @brief Create an empty BlockData for a given block index
     */
    BlockData<T> createBlock(Vec3l block_idx) const {
        Vec3l interior = getBlockInteriorSize(block_idx);
        Vec3l global_start = getBlockGlobalStart(block_idx);
        return BlockData<T>(interior, global_start, block_idx, m_ghost_zone_size);
    }
    
    /**
     * @brief Decompose global data into blocks with ghost cells
     * @param global_data Pointer to contiguous global data array (X fastest, then Y, then Z)
     * @param blocks Output vector of block data
     * @param boundary_value Value to use for ghost cells outside the global domain
     */
    void decompose(const T* global_data, std::vector<BlockData<T>>& blocks, T boundary_value = T(0)) const {
        blocks.clear();
        blocks.reserve(m_num_blocks);
        
        Vec3l global_dims = m_global_grid.XYZ();
        
        // Create all blocks
        for (INDEX_TYPE bk = 0; bk < m_blocks_per_axis[2]; bk++) {
            for (INDEX_TYPE bj = 0; bj < m_blocks_per_axis[1]; bj++) {
                for (INDEX_TYPE bi = 0; bi < m_blocks_per_axis[0]; bi++) {
                    Vec3l block_idx(bi, bj, bk);
                    blocks.push_back(createBlock(block_idx));
                    BlockData<T>& block = blocks.back();
                    
                    // Copy data from global array to block (including ghost cells)
                    for (INDEX_TYPE lz = 0; lz < block.total_size[2]; lz++) {
                        for (INDEX_TYPE ly = 0; ly < block.total_size[1]; ly++) {
                            for (INDEX_TYPE lx = 0; lx < block.total_size[0]; lx++) {
                                // Convert local coords to global coords
                                INDEX_TYPE gx = block.global_start[0] + lx - m_ghost_zone_size;
                                INDEX_TYPE gy = block.global_start[1] + ly - m_ghost_zone_size;
                                INDEX_TYPE gz = block.global_start[2] + lz - m_ghost_zone_size;
                                
                                // Check if within global bounds
                                if (gx >= 0 && gx < global_dims[0] &&
                                    gy >= 0 && gy < global_dims[1] &&
                                    gz >= 0 && gz < global_dims[2]) {
                                    INDEX_TYPE global_idx = gx + gy * global_dims[0] + gz * global_dims[0] * global_dims[1];
                                    block.at(lx, ly, lz) = global_data[global_idx];
                                } else {
                                    // Outside global bounds - use boundary value
                                    block.at(lx, ly, lz) = boundary_value;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    /**
     * @brief Gather interior results from blocks into contiguous global array
     * @param blocks Vector of block data
     * @param global_data Output pointer to contiguous global data array
     */
    void gather(const std::vector<BlockData<T>>& blocks, T* global_data) const {
        Vec3l global_dims = m_global_grid.XYZ();
        
        for (const auto& block : blocks) {
            // Copy interior data from block to global array
            for (INDEX_TYPE iz = 0; iz < block.interior_size[2]; iz++) {
                for (INDEX_TYPE iy = 0; iy < block.interior_size[1]; iy++) {
                    for (INDEX_TYPE ix = 0; ix < block.interior_size[0]; ix++) {
                        INDEX_TYPE gx = block.global_start[0] + ix;
                        INDEX_TYPE gy = block.global_start[1] + iy;
                        INDEX_TYPE gz = block.global_start[2] + iz;
                        
                        INDEX_TYPE global_idx = gx + gy * global_dims[0] + gz * global_dims[0] * global_dims[1];
                        global_data[global_idx] = block.interiorAt(ix, iy, iz);
                    }
                }
            }
        }
    }
    
    /**
     * @brief Exchange ghost cell data between adjacent blocks
     * Call this after decompose() if blocks need updated ghost cells from neighbors
     * @param blocks Vector of block data (modified in place)
     */
    void exchangeGhostCells(std::vector<BlockData<T>>& blocks) const {
        // For now, the decompose() function already fills ghost cells from global data
        // This method would be used in iterative algorithms where blocks are updated
        // and need fresh ghost data from neighbors
        
        // Implementation would copy ghost regions from neighboring blocks' interior
        // This is more complex and depends on the specific use case
        // For single-pass convolution, decompose() is sufficient
    }
    
    /**
     * @brief Get block index from flat index
     */
    Vec3l blockIdxFromFlat(INDEX_TYPE flat_idx) const {
        INDEX_TYPE bx = flat_idx % m_blocks_per_axis[0];
        INDEX_TYPE by = (flat_idx / m_blocks_per_axis[0]) % m_blocks_per_axis[1];
        INDEX_TYPE bz = flat_idx / (m_blocks_per_axis[0] * m_blocks_per_axis[1]);
        return Vec3l(bx, by, bz);
    }
    
    /**
     * @brief Get flat index from block index
     */
    INDEX_TYPE flatIdxFromBlock(Vec3l block_idx) const {
        return block_idx[0] + block_idx[1] * m_blocks_per_axis[0] + 
               block_idx[2] * m_blocks_per_axis[0] * m_blocks_per_axis[1];
    }
};

/**
 * @brief 2D version of BlockData
 */
template<typename T>
struct BlockData2D {
    RegularGrid2D local_grid;
    std::vector<T> data;
    Vec2l global_start;
    Vec2l block_index;
    int ghost_zone_size;
    Vec2l interior_size;
    Vec2l total_size;
    
    BlockData2D() : local_grid(Vec2l(1,1), Vec2b(false,false)), ghost_zone_size(0) {}
    
    BlockData2D(Vec2l interior, Vec2l global_start_, Vec2l block_idx, int ghost_size)
        : interior_size(interior),
          global_start(global_start_),
          block_index(block_idx),
          ghost_zone_size(ghost_size),
          total_size(interior[0] + 2*ghost_size, interior[1] + 2*ghost_size),
          local_grid(total_size, Vec2b(false, false))
    {
        data.resize(total_size[0] * total_size[1], T(0));
    }
    
    INDEX_TYPE localIndex(INDEX_TYPE x, INDEX_TYPE y) const {
        return x + y * total_size[0];
    }
    
    INDEX_TYPE interiorToLocalIndex(INDEX_TYPE ix, INDEX_TYPE iy) const {
        return localIndex(ix + ghost_zone_size, iy + ghost_zone_size);
    }
    
    T& at(INDEX_TYPE x, INDEX_TYPE y) { return data[localIndex(x, y)]; }
    const T& at(INDEX_TYPE x, INDEX_TYPE y) const { return data[localIndex(x, y)]; }
    
    T& interiorAt(INDEX_TYPE ix, INDEX_TYPE iy) { return data[interiorToLocalIndex(ix, iy)]; }
    const T& interiorAt(INDEX_TYPE ix, INDEX_TYPE iy) const { return data[interiorToLocalIndex(ix, iy)]; }
    
    T* dataPtr() { return data.data(); }
    const T* dataPtr() const { return data.data(); }
};

/**
 * @brief 2D version of BlockDecomposition
 */
template<typename T>
class BlockDecomposition2D {
private:
    RegularGrid2D m_global_grid;
    Vec2l m_block_size;
    Vec2l m_blocks_per_axis;
    int m_ghost_zone_size;
    INDEX_TYPE m_num_blocks;
    
public:
    BlockDecomposition2D(Vec2l global_dims, Vec2l block_size, int ghost_zone_size)
        : m_global_grid(global_dims, Vec2b(false, false)),
          m_block_size(block_size),
          m_ghost_zone_size(ghost_zone_size)
    {
        for (int i = 0; i < 2; i++) {
            m_blocks_per_axis[i] = (global_dims[i] + block_size[i] - 1) / block_size[i];
        }
        m_num_blocks = m_blocks_per_axis[0] * m_blocks_per_axis[1];
    }
    
    struct FromKernelSize {};
    BlockDecomposition2D(Vec2l global_dims, Vec2l block_size, int kernel_size, FromKernelSize)
        : BlockDecomposition2D(global_dims, block_size, (kernel_size - 1) / 2)
    {}
    
    const RegularGrid2D& globalGrid() const { return m_global_grid; }
    Vec2l blocksPerAxis() const { return m_blocks_per_axis; }
    INDEX_TYPE numBlocks() const { return m_num_blocks; }
    int ghostZoneSize() const { return m_ghost_zone_size; }
    Vec2l blockSize() const { return m_block_size; }
    
    Vec2l getBlockInteriorSize(Vec2l block_idx) const {
        Vec2l interior;
        Vec2l global_dims = m_global_grid.XY();
        for (int i = 0; i < 2; i++) {
            INDEX_TYPE start = block_idx[i] * m_block_size[i];
            INDEX_TYPE end = std::min(start + m_block_size[i], global_dims[i]);
            interior[i] = end - start;
        }
        return interior;
    }
    
    Vec2l getBlockGlobalStart(Vec2l block_idx) const {
        return Vec2l(block_idx[0] * m_block_size[0], block_idx[1] * m_block_size[1]);
    }
    
    BlockData2D<T> createBlock(Vec2l block_idx) const {
        return BlockData2D<T>(getBlockInteriorSize(block_idx), getBlockGlobalStart(block_idx), 
                              block_idx, m_ghost_zone_size);
    }
    
    void decompose(const T* global_data, std::vector<BlockData2D<T>>& blocks, T boundary_value = T(0)) const {
        blocks.clear();
        blocks.reserve(m_num_blocks);
        Vec2l global_dims = m_global_grid.XY();
        
        for (INDEX_TYPE bj = 0; bj < m_blocks_per_axis[1]; bj++) {
            for (INDEX_TYPE bi = 0; bi < m_blocks_per_axis[0]; bi++) {
                Vec2l block_idx(bi, bj);
                blocks.push_back(createBlock(block_idx));
                BlockData2D<T>& block = blocks.back();
                
                for (INDEX_TYPE ly = 0; ly < block.total_size[1]; ly++) {
                    for (INDEX_TYPE lx = 0; lx < block.total_size[0]; lx++) {
                        INDEX_TYPE gx = block.global_start[0] + lx - m_ghost_zone_size;
                        INDEX_TYPE gy = block.global_start[1] + ly - m_ghost_zone_size;
                        
                        if (gx >= 0 && gx < global_dims[0] && gy >= 0 && gy < global_dims[1]) {
                            block.at(lx, ly) = global_data[gx + gy * global_dims[0]];
                        } else {
                            block.at(lx, ly) = boundary_value;
                        }
                    }
                }
            }
        }
    }
    
    void gather(const std::vector<BlockData2D<T>>& blocks, T* global_data) const {
        Vec2l global_dims = m_global_grid.XY();
        for (const auto& block : blocks) {
            for (INDEX_TYPE iy = 0; iy < block.interior_size[1]; iy++) {
                for (INDEX_TYPE ix = 0; ix < block.interior_size[0]; ix++) {
                    INDEX_TYPE gx = block.global_start[0] + ix;
                    INDEX_TYPE gy = block.global_start[1] + iy;
                    global_data[gx + gy * global_dims[0]] = block.interiorAt(ix, iy);
                }
            }
        }
    }
};

} // namespace GInt

#endif // BLOCK_DECOMPOSITION_H
