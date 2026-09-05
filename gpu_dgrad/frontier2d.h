#ifndef DGRAD_FRONTIER2D_H
#define DGRAD_FRONTIER2D_H

// Stage A1/A2: host-side frontier graph of the tiled successor resolver.
// Plain C++ (no CUDA types) so later stages and validators can reuse it.
//
// After the in-tile doubling, every cell of a tile carries a root word:
// TERMINAL(v) (a terminal inside the tile) or FRONTIER(r) where r is the
// position, on the tile's OUTER ring (the one-cell layer around its lattice
// box), of the exit target the chain leaves through. A successor is never
// more than one step outside the box (a +-1 flow step, or the corner of a
// voxel the numeric hop was in when it crossed the box), so the outer ring
// enumerates every possible target; the target's global index is implied by
// (tile, r), and an exit target is always a cell on the INNER ring of the
// neighbouring tile, whose root word phase 1 downloads. The frontier slots
// (tile, r) therefore form a successor function on a perimeter-sized set,
// rooted at terminals; ResolveFrontier2D follows it with memoization. The
// successor field is a forest by construction, so no cycles exist; a step
// cap guards against a corrupt input anyway.

#include "tile2d.h"

#include <cstdint>
#include <vector>

namespace GInt {
namespace gpu {

// Word encoding shared with the device (successor2d.cu keeps a copy).
enum : uint32_t {
    kWordTagShift = 30,
    kWordValueMask = 0x3FFFFFFFu,
    kWordTagSucc = 0u,
    kWordTagTerminal = 2u,
    kWordTagFrontier = 3u,
};
inline uint32_t WordTag(int32_t w) { return ((uint32_t)w) >> kWordTagShift; }
inline uint32_t WordValue(int32_t w) { return ((uint32_t)w) & kWordValueMask; }

// How a TERMINAL(v) word maps to its final payload.
enum PayloadMode {
    kPayloadLatticeIndex = 0,  // v = local lattice index of the terminal cell
    kPayloadCompactId = 1,     // v = scan-order rank; payload = tile base + v
    kPayloadDirect = 2,        // v IS the payload (e.g. a living region id)
};

struct TileRingRecord2D {
    int32_t count = 0;               // terminals owned by the tile
    int64_t base = 0;                // compact-id base (prefix sum in tile order)
    int w = 0, h = 0;                // lattice box of the tile
    int64_t ox = 0, oy = 0;
    int R = 0;                       // inner ring cells
    int Rout = 0;                    // outer ring slots
    std::vector<int32_t> word;       // [R] root word of each inner ring cell
    std::vector<int64_t> resolved;   // [Rout] final payload of each outer slot, -1 = no target
};

struct FrontierGraph2D {
    int64_t NX = 0, NY = 0;
    TileGrid2D grid;
    std::vector<TileRingRecord2D> tiles;   // [numTiles]

    int64_t TileOfLattice(int64_t g, int& lx, int& ly) const {
        const int64_t gx = g % NX, gy = g / NX;
        const int tx = (int)(gx / grid.T), ty = (int)(gy / grid.T);
        const int64_t t = (int64_t)tx + (int64_t)ty * grid.ntx;
        lx = (int)(gx - tiles[(size_t)t].ox);
        ly = (int)(gy - tiles[(size_t)t].oy);
        return t;
    }
    // global lattice index of outer slot r of tile t, or -1 outside the domain
    int64_t ExitTarget(int64_t t, int r) const {
        const TileRingRecord2D& tr = tiles[(size_t)t];
        int x, y;
        RingToLocal2D(r, tr.w + 2, tr.h + 2, x, y);
        const int64_t gx = tr.ox + x - 1, gy = tr.oy + y - 1;
        if (gx < 0 || gy < 0 || gx >= NX || gy >= NY) return -1;
        return gx + gy * NX;
    }
};

// Payload of a TERMINAL word of tile t in the requested id space.
inline int64_t TerminalPayload2D(const FrontierGraph2D& g, int64_t t, uint32_t v,
                                 PayloadMode mode) {
    const TileRingRecord2D& tr = g.tiles[(size_t)t];
    if (mode == kPayloadDirect) return (int64_t)v;
    if (mode == kPayloadCompactId) return tr.base + (int64_t)v;
    const int64_t lx = (int64_t)v % tr.w, ly = (int64_t)v / tr.w;
    return (tr.oy + ly) * g.NX + (tr.ox + lx);
}

// Resolves every outer slot with a target to a final payload. Returns false
// on a malformed graph (a target not on the neighbour's inner ring, an
// unresolved word there, or a chain longer than the slot count). steps
// receives the number of slot hops taken.
inline bool ResolveFrontier2D(FrontierGraph2D& g, PayloadMode mode, int64_t* steps) {
    int64_t totalSlots = 0;
    for (size_t t = 0; t < g.tiles.size(); t++) {
        TileRingRecord2D& tr = g.tiles[t];
        tr.resolved.assign((size_t)tr.Rout, -2);   // -2 = not yet visited
        for (int r = 0; r < tr.Rout; r++)
            if (g.ExitTarget((int64_t)t, r) < 0) tr.resolved[(size_t)r] = -1;
        totalSlots += tr.Rout;
    }
    int64_t hops = 0;
    std::vector<std::pair<int64_t, int>> path;
    for (size_t t0 = 0; t0 < g.tiles.size(); t0++) {
        TileRingRecord2D& tr0 = g.tiles[t0];
        for (int r0 = 0; r0 < tr0.Rout; r0++) {
            if (tr0.resolved[(size_t)r0] != -2) continue;
            path.clear();
            int64_t t = (int64_t)t0;
            int r = r0;
            int64_t payload = -1;
            while (true) {
                TileRingRecord2D& tr = g.tiles[(size_t)t];
                const int64_t memo = tr.resolved[(size_t)r];
                if (memo >= 0) { payload = memo; break; }
                if (memo == -1) return false;   // chain into a slot with no target
                path.push_back(std::make_pair(t, r));
                if ((int64_t)path.size() > totalSlots + 1) return false;
                const int64_t gidx = g.ExitTarget(t, r);
                int lx = 0, ly = 0;
                const int64_t tn = g.TileOfLattice(gidx, lx, ly);
                const TileRingRecord2D& trn = g.tiles[(size_t)tn];
                const int rn = RingIndex2D(lx, ly, trn.w, trn.h);
                if (rn < 0) return false;
                hops++;
                const int32_t wn = trn.word[(size_t)rn];
                if (WordTag(wn) == kWordTagTerminal) {
                    payload = TerminalPayload2D(g, tn, WordValue(wn), mode);
                    break;
                }
                if (WordTag(wn) != kWordTagFrontier) return false;
                t = tn;
                r = (int)WordValue(wn);
            }
            for (size_t i = 0; i < path.size(); i++)
                g.tiles[(size_t)path[i].first].resolved[(size_t)path[i].second] = payload;
        }
    }
    if (steps) *steps = hops;
    return true;
}

} // namespace gpu
} // namespace GInt

#endif
