#ifndef DGRAD3D_TABLES_H
#define DGRAD3D_TABLES_H

// Generates the 27-slot topology tables of the 3D GPU lower-star kernel from
// the iterators the CPU algorithm runs:
//   - Explicit3x3x3SmallRegularGrid (the sliding window of
//     SlidingWindowRobinsNoalloc::ComputePairing_sliding, interior vertices),
//   - TopologicalRegularGrid3D (MyRobinsNoalloc::ComputeLowerStar, used by the
//     sliding driver for boundary vertices and by the constrained pass).
// The tables are built from the small grid and then required to be identical
// to what the big mesh's iterators produce around an interior vertex, so one
// table serves both CPU code paths. The other half of the contract is the
// byte-exact output compare in cls3d_validate.
//
// Header-only and templated on the mesh types so gpu_dgrad stays free of gi_*
// includes; the caller instantiates it with the concrete meshes.
//
// Returns false if any order or convention the kernel relies on differs.

#include "dgrad_gpu_api.h"

#include <cstring>

namespace GInt {
namespace gpu {

namespace dgrad3d_tables_detail {
inline int SlotOf(int dx, int dy, int dz) { return (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1); }
} // namespace dgrad3d_tables_detail

// small: an Explicit3x3x3SmallRegularGrid (5x5x5 cells; the window vertex is
// at cell coords (2,2,2)). VecType: the 3-vector of long the mesh API takes.
template <class SmallMeshType, class VecType>
bool BuildDgrad3DTablesFromSmallMesh(SmallMeshType& sm, Dgrad3DTables& t) {
    using namespace dgrad3d_tables_detail;
    typedef long long IndexT;
    VecType vc;
    vc[0] = 2; vc[1] = 2; vc[2] = 2;
    const IndexT base = (IndexT)sm.coords2Cellid(vc);
    IndexT off[27];
    for (int s = 0; s < 27; s++) {
        off[s] = (IndexT)sm.get27NeighborOffset(s);
        const int dx = s % 3 - 1, dy = (s / 3) % 3 - 1, dz = s / 9 - 1;
        if (off[s] != dx + 5 * dy + 25 * dz) return false;  // slot convention
    }
    // AdjacentCellsIterator order == ascending slot
    {
        typename SmallMeshType::AdjacentCellsIterator it(&sm);
        int s = 0;
        for (it.begin(base); it.valid(); it.advance(), s++)
            if (s >= 27 || (IndexT)it.value() != base + off[s]) return false;
        if (s != 27) return false;
    }
    std::memset(&t, 0, sizeof(t));
    for (int s = 0; s < 27; s++) {
        const IndexT c = base + off[s];
        t.dim[s] = (uint8_t)sm.dimension(c);
        const int dx = s % 3 - 1, dy = (s / 3) % 3 - 1, dz = s / 9 - 1;
        if (t.dim[s] != (uint8_t)((dx != 0) + (dy != 0) + (dz != 0))) return false;

        typename SmallMeshType::FacetsIterator fit(&sm);
        int nf = 0;
        for (fit.begin(c); fit.valid(); fit.advance()) {
            const IndexT d = (IndexT)fit.value() - base;
            for (int q = 0; q < 27; q++)
                if (off[q] == d) {
                    if (nf >= 6) return false;
                    t.fac[s][nf++] = (uint8_t)q;
                }
            // facets outside the window can never be in this lower star
        }
        t.fac_count[s] = (uint8_t)nf;

        typename SmallMeshType::CofacetsIterator cfit(&sm);
        int nc = 0;
        for (cfit.begin(c); cfit.valid(); cfit.advance()) {
            const IndexT d = (IndexT)cfit.value() - base;
            for (int q = 0; q < 27; q++)
                if (off[q] == d) {
                    if (nc >= 6) return false;
                    t.cof[s][nc++] = (uint8_t)q;
                }
            // Out-of-window values are skipped: the small grid's iterators
            // visit one entry past each id list (m_end = 1 + count), i.e. the
            // zero padding or the next row's count, an id below 7 that is
            // never in a lower star. Every true cofacet of a window cell is
            // in the window; CheckDgrad3DTablesAgainstMesh confirms the list.
        }
        t.cof_count[s] = (uint8_t)nc;

        // vertices: grid offset e with e_i in {0, d_i}; window position w(e).
        // Cross-checked against the small mesh's vertex numbering (the index
        // the CPU's Before() tie-break uses).
        int nv = 0;
        for (int ez = (dz < 0 ? dz : 0); ez <= (dz > 0 ? dz : 0); ez++)
            for (int ey = (dy < 0 ? dy : 0); ey <= (dy > 0 ? dy : 0); ey++)
                for (int ex = (dx < 0 ? dx : 0); ex <= (dx > 0 ? dx : 0); ex++) {
                    const int w = SlotOf(ex, ey, ez);
                    VecType vv;
                    vv[0] = 2 + 2 * ex; vv[1] = 2 + 2 * ey; vv[2] = 2 + 2 * ez;
                    if ((int)sm.VertexNumberFromCellID(sm.coords2Cellid(vv)) != w)
                        return false;
                    t.vert[s][nv++] = (uint8_t)w;
                }
        t.vert_count[s] = (uint8_t)nv;
    }
    return true;
}

// Builds the same tables from a TopologicalRegularGrid3D (at least 3 vertices
// per axis, non-periodic) around the vertex at grid (1,1,1) and compares them
// with 't'. Returns true iff identical.
template <class MeshType, class VecType>
bool CheckDgrad3DTablesAgainstMesh(MeshType& mesh, const Dgrad3DTables& t) {
    typedef long long IndexT;
    VecType vc;
    vc[0] = 2; vc[1] = 2; vc[2] = 2;
    const IndexT base = (IndexT)mesh.coords2Cellid(vc);
    IndexT off[27];
    for (int s = 0; s < 27; s++) off[s] = (IndexT)mesh.get27NeighborOffset(s);
    {
        typename MeshType::AdjacentCellsIterator it(&mesh);
        int s = 0;
        for (it.begin(base); it.valid(); it.advance(), s++)
            if (s >= 27 || (IndexT)it.value() != base + off[s]) return false;
        if (s != 27) return false;
    }
    for (int s = 0; s < 27; s++) {
        const IndexT c = base + off[s];
        if ((int)mesh.dimension(c) != t.dim[s]) return false;
        typename MeshType::FacetsIterator fit(&mesh);
        int nf = 0;
        for (fit.begin(c); fit.valid(); fit.advance()) {
            const IndexT d = (IndexT)fit.value() - base;
            for (int q = 0; q < 27; q++)
                if (off[q] == d) {
                    if (nf >= t.fac_count[s] || t.fac[s][nf] != q) return false;
                    nf++;
                }
        }
        if (nf != t.fac_count[s]) return false;
        typename MeshType::CofacetsIterator cfit(&mesh);
        int nc = 0;
        for (cfit.begin(c); cfit.valid(); cfit.advance()) {
            const IndexT d = (IndexT)cfit.value() - base;
            for (int q = 0; q < 27; q++)
                if (off[q] == d) {
                    if (nc >= t.cof_count[s] || t.cof[s][nc] != q) return false;
                    nc++;
                }
        }
        if (nc != t.cof_count[s]) return false;
    }
    return true;
}

} // namespace gpu
} // namespace GInt

#endif
