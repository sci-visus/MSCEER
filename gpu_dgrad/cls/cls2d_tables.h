#ifndef DGRAD2D_TABLES_H
#define DGRAD2D_TABLES_H

// Generates the 9-slot topology tables for the 2D GPU kernels from a real
// TopologicalRegularGrid2D, by running the SAME iterators the CPU algorithm
// uses - the single source of truth binding the two implementations (the
// other half of the contract is the byte-exact output compare).
//
// Header-only and templated on the mesh type so gpu_dgrad itself stays free
// of gi_* includes; the caller (validation harness, msc_2d_lib) instantiates
// it with the concrete mesh.
//
// Returns false if the mesh's offset enumeration does not match the slot
// conventions the kernels assume (slot s = (dx+1)+3*(dy+1), symmetric table).
// Requires a mesh of at least 3x3 vertices (an interior vertex must exist).

#include "dgrad_gpu_api.h"

#include <cstring>

namespace GInt {
namespace gpu {

template <class MeshType>
bool BuildDgrad2DTablesFromMesh(MeshType& mesh, Dgrad2DTables& t) {
    typedef decltype(mesh.numCells()) IndexT;
    const IndexT W = mesh.numCellsAxis(0);
    if (mesh.numCellsAxis(0) < 5 || mesh.numCellsAxis(1) < 5) return false;

    // interior vertex at mesh coords (2,2) = grid vertex (1,1)
    typename MeshType::ICoordType vc;
    vc[0] = 2;
    vc[1] = 2;
    const IndexT base = mesh.coords2Cellid(vc);

    IndexT off9[9];
    for (int s = 0; s < 9; s++) {
        off9[s] = mesh.get8NeighborOffset(s);
        const IndexT expect = (s % 3 - 1) + (IndexT)(s / 3 - 1) * W;
        if (off9[s] != expect) return false;
    }
    for (int s = 0; s < 9; s++)
        if (off9[s] != -off9[8 - s]) return false; // maxbyte==8-s membership test

    std::memset(&t, 0, sizeof(t));
    for (int s = 0; s < 9; s++) {
        const IndexT c = base + off9[s];
        t.dim[s] = (uint8_t)mesh.dimension(c);

        typename MeshType::FacetsIterator fit(&mesh);
        int nf = 0;
        for (fit.begin(c); fit.valid(); fit.advance()) {
            const IndexT d = fit.value() - base;
            for (int q = 0; q < 9; q++)
                if (off9[q] == d) t.fac[s][nf++] = (uint8_t)q;
            // facets outside the 3x3 can never be in a lower star; skipped
        }
        t.fac_count[s] = (uint8_t)nf;

        typename MeshType::CofacetsIterator cfit(&mesh);
        int nc = 0;
        for (cfit.begin(c); cfit.valid(); cfit.advance()) {
            const IndexT d = cfit.value() - base;
            for (int q = 0; q < 9; q++)
                if (off9[q] == d) t.cof[s][nc++] = (uint8_t)q;
        }
        t.cof_count[s] = (uint8_t)nc;
    }
    t.dir_px = mesh.Compress6NeighborOffsetToByte(base, base + 1);
    t.dir_mx = mesh.Compress6NeighborOffsetToByte(base, base - 1);
    t.dir_py = mesh.Compress6NeighborOffsetToByte(base, base + W);
    t.dir_my = mesh.Compress6NeighborOffsetToByte(base, base - W);
    return true;
}

} // namespace gpu
} // namespace GInt

#endif
