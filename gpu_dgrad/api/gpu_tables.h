#ifndef GPU_TABLES_H
#define GPU_TABLES_H

// Local topology tables and timing structs shared by the 2D and 3D GPU entry points.
// They describe the grid, not a gradient method. Plain C++.

#include <cstdint>

namespace GInt {
namespace gpu {

// Local topology of the 3x3 cell neighborhood of a vertex in the doubled grid.
// Slot s = (dx+1) + 3*(dy+1) for dx,dy in {-1,0,1} - the same enumeration as
// TopologicalRegularGrid2D::m_adjacent_cell_offsets, so slot 4 is the vertex.
// The tables MUST be generated from the real mesh's iterators (see the
// validation harness) so the device kernel provably replicates the CPU
// iteration orders; they are never hand-typed.
struct Dgrad2DTables {
    uint8_t dim[9];                            // cell dimension per slot
    uint8_t fac_count[9]; uint8_t fac[9][2];   // facets inside the 3x3, FacetsIterator order
    uint8_t cof_count[9]; uint8_t cof[9][4];   // cofacets inside the 3x3, CofacetsIterator order
    uint8_t dir_px, dir_mx, dir_py, dir_my;    // GradBitfield pair codes for +x,-x,+y,-y
};

struct Dgrad2DTimings {
    float h2d_ms;
    float kernel_ms;
    float d2h_ms;
};

// Local topology of the 27 cells around a vertex of the doubled grid. Slot
// s = (dx+1) + 3*(dy+1) + 9*(dz+1), the AdjacentCellsIterator order of both
// TopologicalRegularGrid3D and Explicit3x3x3SmallRegularGrid, so slot 13 is
// the vertex. Vertex-window positions w use the same encoding for the grid
// offset of a vertex (the small RegularGrid3D index of the CPU sliding
// window), so w orders the 27 window vertices exactly as the global vertex
// index does. Generated from the real iterators (cls3d_tables.h), never
// hand-typed.
struct Dgrad3DTables {
    uint8_t dim[27];
    uint8_t fac_count[27]; uint8_t fac[27][6];   // facets inside the window, FacetsIterator order
    uint8_t cof_count[27]; uint8_t cof[27][6];   // cofacets inside the window, CofacetsIterator order
    uint8_t vert_count[27]; uint8_t vert[27][8]; // vertices of the cell as window positions
};

struct Dgrad3DTimings {
    float h2d_ms;      // summed per-batch event windows
    float kernel_ms;
    float d2h_ms;
    float wall_ms;     // end to end, host copies included
};

} // namespace gpu
} // namespace GInt

#endif
