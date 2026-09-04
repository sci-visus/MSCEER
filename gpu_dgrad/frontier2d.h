#ifndef DGRAD_FRONTIER2D_H
#define DGRAD_FRONTIER2D_H

// Stage A1: host-side frontier graph of the tiled successor resolver.
// Plain C++ (no CUDA types) so later stages and validators can reuse it.
//
// After the in-tile doubling, every ring cell of a tile carries a root word:
// TERMINAL(v) (a terminal inside the tile) or FRONTIER(r) (its chain leaves
// the tile through ring cell r, whose exit target is the neighbour lattice
// point exitGlobal[r]). The exit target is at most one step outside, so it is
// a ring cell of the neighbouring tile, whose root word is again known. The
// frontier slots (tile, r) therefore form a successor function on a
// perimeter-sized set, rooted at terminals; ResolveFrontier2D follows it
// with memoization. The discrete flow is a forest, so no cycles exist; a
// step cap guards against a corrupt input anyway.

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

struct TileRingRecord2D {
    int32_t count = 0;               // terminals owned by the tile
    int64_t base = 0;                // compact-id base (prefix sum in tile order)
    int w = 0, h = 0;                // lattice box of the tile
    int64_t ox = 0, oy = 0;
    int R = 0;                       // ring cells
    std::vector<int32_t> word;       // [R] root word of each ring cell
    std::vector<int64_t> exitGlobal; // [R] exit target (global lattice index) or -1
    std::vector<int64_t> resolved;   // [R] final payload of exiting slots, -1 otherwise
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
};

// Payload of a TERMINAL word of tile t in the requested id space.
inline int64_t TerminalPayload2D(const FrontierGraph2D& g, int64_t t, uint32_t v,
                                 bool compactIds) {
    const TileRingRecord2D& tr = g.tiles[(size_t)t];
    if (compactIds) return tr.base + (int64_t)v;   // v = rank inside the tile
    // v = local lattice index inside the tile box
    const int64_t lx = (int64_t)v % tr.w, ly = (int64_t)v / tr.w;
    return (tr.oy + ly) * g.NX + (tr.ox + lx);
}

// Resolves every exiting frontier slot to a final payload. Returns false on
// a malformed graph (exit target not on a ring, or a chain longer than the
// slot count). steps receives the number of slot hops taken.
inline bool ResolveFrontier2D(FrontierGraph2D& g, bool compactIds, int64_t* steps) {
    int64_t totalSlots = 0;
    for (size_t t = 0; t < g.tiles.size(); t++) {
        TileRingRecord2D& tr = g.tiles[t];
        tr.resolved.assign((size_t)tr.R, -2);   // -2 = not yet visited
        for (int r = 0; r < tr.R; r++)
            if (tr.exitGlobal[(size_t)r] < 0) tr.resolved[(size_t)r] = -1;
        totalSlots += tr.R;
    }
    int64_t hops = 0;
    std::vector<std::pair<int64_t, int>> path;
    for (size_t t0 = 0; t0 < g.tiles.size(); t0++) {
        TileRingRecord2D& tr0 = g.tiles[t0];
        for (int r0 = 0; r0 < tr0.R; r0++) {
            if (tr0.resolved[(size_t)r0] != -2) continue;
            path.clear();
            int64_t t = (int64_t)t0;
            int r = r0;
            int64_t payload = -1;
            while (true) {
                TileRingRecord2D& tr = g.tiles[(size_t)t];
                const int64_t memo = tr.resolved[(size_t)r];
                if (memo >= 0) { payload = memo; break; }
                path.push_back(std::make_pair(t, r));
                if ((int64_t)path.size() > totalSlots + 1) return false;
                // hop to the exit target's tile and ring cell
                const int64_t gidx = tr.exitGlobal[(size_t)r];
                if (gidx < 0) return false;
                int lx = 0, ly = 0;
                const int64_t tn = g.TileOfLattice(gidx, lx, ly);
                const TileRingRecord2D& trn = g.tiles[(size_t)tn];
                const int rn = RingIndex2D(lx, ly, trn.w, trn.h);
                if (rn < 0) return false;
                hops++;
                const int32_t wn = trn.word[(size_t)rn];
                if (WordTag(wn) == kWordTagTerminal) {
                    payload = TerminalPayload2D(g, tn, WordValue(wn), compactIds);
                    break;
                }
                if (WordTag(wn) != kWordTagFrontier) return false;
                t = tn;
                r = (int)WordValue(wn);
                // the exiting ring cell's own chain: continue from (tn, r)
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
