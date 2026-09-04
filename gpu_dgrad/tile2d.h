#ifndef DGRAD_TILE2D_H
#define DGRAD_TILE2D_H

// 2D tile geometry for the tiled GPU pipeline (Stage A0 of the accurate-
// gradient plan). Plain C++: no CUDA types, includable by any host TU.
//
// The ALGORITHMIC tile is the unit of ownership: tile (tx,ty) owns the grid
// vertices [ox, ox+tw) x [oy, oy+th), and a doubled-lattice cell c is owned by
// the tile owning vertex floor(c/2). That is exactly the 1-high "bridge"
// convention of steepest_lstar_tiled/localized_morse_complex.h: the owned cell
// box of a tile is 2*tw wide (its last cell is the edge bridging to the next
// tile), 2*tw-1 for the last tile of an axis. Every cell is owned exactly once.
//
// How many tiles are resident on the device, and how many go into one launch,
// are scheduling choices (TiledOptions2D); they never change the output. The
// standing gate of every tiled stage is: tiled output == in-core output,
// byte for byte, for every tile size and residency.

#include <cstdint>
#include <functional>

namespace GInt {
namespace gpu {

struct TileGrid2D {
    int64_t X = 0, Y = 0;    // image size in vertices
    int T = 0;               // nominal tile size (vertices per axis)
    int ntx = 0, nty = 0;    // tiles per axis
    int64_t numTiles = 0;
};

struct TileDesc2D {
    int tx = 0, ty = 0;      // tile coordinates
    int64_t ox = 0, oy = 0;  // first owned vertex
    int tw = 0, th = 0;      // owned vertices per axis
};

// An axis-aligned box in some lattice: [x0, x0+w) x [y0, y0+h).
struct Box2D {
    int64_t x0 = 0, y0 = 0;
    int w = 0, h = 0;
};

inline TileGrid2D MakeTileGrid2D(int64_t X, int64_t Y, int T) {
    TileGrid2D g;
    g.X = X;
    g.Y = Y;
    g.T = T < 1 ? 1 : T;
    g.ntx = (int)((X + g.T - 1) / g.T);
    g.nty = (int)((Y + g.T - 1) / g.T);
    g.numTiles = (int64_t)g.ntx * g.nty;
    return g;
}

// tileId = tx + ty*ntx (row-major over tiles).
inline TileDesc2D TileAt(const TileGrid2D& g, int64_t tileId) {
    TileDesc2D d;
    d.tx = (int)(tileId % g.ntx);
    d.ty = (int)(tileId / g.ntx);
    d.ox = (int64_t)d.tx * g.T;
    d.oy = (int64_t)d.ty * g.T;
    d.tw = (int)((g.X - d.ox) < g.T ? (g.X - d.ox) : g.T);
    d.th = (int)((g.Y - d.oy) < g.T ? (g.Y - d.oy) : g.T);
    return d;
}

// Owned doubled-lattice cell box: c owned iff floor(c/2) in [ox, ox+tw), i.e.
// x in [2ox, 2(ox+tw)) clipped to the lattice extent 2X-1 (exclusive): 2tw
// cells wide, 2tw-1 for the last tile of the axis.
inline Box2D OwnedCellBox2D(const TileGrid2D& g, const TileDesc2D& d) {
    Box2D b;
    b.x0 = 2 * d.ox;
    b.y0 = 2 * d.oy;
    const int64_t cx1 = 2 * (d.ox + d.tw) < 2 * g.X - 1 ? 2 * (d.ox + d.tw) : 2 * g.X - 1;
    const int64_t cy1 = 2 * (d.oy + d.th) < 2 * g.Y - 1 ? 2 * (d.oy + d.th) : 2 * g.Y - 1;
    b.w = (int)(cx1 - b.x0);
    b.h = (int)(cy1 - b.y0);
    return b;
}

// Vertices whose lower stars the kernel evaluates for this tile: owned + a
// 1-vertex ring (the bridge cells lie in the lower stars of the neighbouring
// tile's first vertices), clipped to the domain.
inline Box2D ComputedVertexBox2D(const TileGrid2D& g, const TileDesc2D& d) {
    Box2D b;
    const int64_t x0 = d.ox - 1 < 0 ? 0 : d.ox - 1;
    const int64_t y0 = d.oy - 1 < 0 ? 0 : d.oy - 1;
    const int64_t x1 = d.ox + d.tw + 1 > g.X ? g.X : d.ox + d.tw + 1;
    const int64_t y1 = d.oy + d.th + 1 > g.Y ? g.Y : d.oy + d.th + 1;
    b.x0 = x0;
    b.y0 = y0;
    b.w = (int)(x1 - x0);
    b.h = (int)(y1 - y0);
    return b;
}

// Function values the kernel needs: computed vertices + their 1-ring (the
// fused lower-star membership test reads the 3x3 neighbourhood), clipped.
inline Box2D StagedVertexBox2D(const TileGrid2D& g, const TileDesc2D& d) {
    const Box2D c = ComputedVertexBox2D(g, d);
    Box2D b;
    const int64_t x0 = c.x0 - 1 < 0 ? 0 : c.x0 - 1;
    const int64_t y0 = c.y0 - 1 < 0 ? 0 : c.y0 - 1;
    const int64_t x1 = c.x0 + c.w + 1 > g.X ? g.X : c.x0 + c.w + 1;
    const int64_t y1 = c.y0 + c.h + 1 > g.Y ? g.Y : c.y0 + c.h + 1;
    b.x0 = x0;
    b.y0 = y0;
    b.w = (int)(x1 - x0);
    b.h = (int)(y1 - y0);
    return b;
}

// Callback I/O so the driver never assumes the image or the output fit in
// host memory. A source fills a box of vertex values (row-major, x-fastest,
// with the given row pitch in elements); a sink receives a box of doubled-
// lattice cell bytes. Boxes are always inside the domain. Return false to
// abort the pipeline.
typedef std::function<bool(int64_t x0, int64_t y0, int w, int h, float* dst, int pitch)>
    TileSource2D;
typedef std::function<bool(int64_t cx0, int64_t cy0, int w, int h, const uint8_t* src,
                           int pitch)>
    ByteSink2D;

// In-core adapters over host arrays (X*Y floats; (2X-1)*(2Y-1) bytes).
inline TileSource2D HostImageSource2D(const float* values, int64_t X) {
    return [values, X](int64_t x0, int64_t y0, int w, int h, float* dst, int pitch) {
        for (int r = 0; r < h; r++) {
            const float* s = values + (y0 + r) * X + x0;
            float* d = dst + (int64_t)r * pitch;
            for (int c = 0; c < w; c++) d[c] = s[c];
        }
        return true;
    };
}

inline ByteSink2D HostGradSink2D(uint8_t* grad, int64_t X) {
    const int64_t W = 2 * X - 1;
    return [grad, W](int64_t cx0, int64_t cy0, int w, int h, const uint8_t* src, int pitch) {
        for (int r = 0; r < h; r++) {
            uint8_t* d = grad + (cy0 + r) * W + cx0;
            const uint8_t* s = src + (int64_t)r * pitch;
            for (int c = 0; c < w; c++) d[c] = s[c];
        }
        return true;
    };
}

struct TiledOptions2D {
    int slots = 3;           // resident tile batches (pinned + device buffers)
    int tilesPerLaunch = 1;  // tiles staged into one launch (blockIdx.z)
    int block_shape = 0;     // 0 = 16x16, 1 = 32x8, 2 = 8x8 (as the batched path)
};

struct TileStats2D {
    int64_t tiles = 0;
    int64_t launches = 0;
    int64_t bytesH2D = 0;
    int64_t bytesD2H = 0;
    float h2d_ms = 0.f;      // per-batch event windows, summed (overlapped legs
    float kernel_ms = 0.f;   //  may add up to more than wall_ms)
    float d2h_ms = 0.f;
    float wall_ms = 0.f;
    int64_t deviceFreeMin = -1;  // cudaMemGetInfo sampled once per batch
    int64_t deviceFreeMax = -1;
    int64_t deviceBytesAllocated = 0;  // this call's own device allocations
    int64_t hostPinnedBytes = 0;       // this call's own pinned allocations
};

} // namespace gpu
} // namespace GInt

#endif
