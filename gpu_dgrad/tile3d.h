#ifndef DGRAD_TILE3D_H
#define DGRAD_TILE3D_H

// 3D tile geometry for the tiled GPU lower-star gradient (dgrad3d_*). Plain
// C++: no CUDA types, includable by any host TU.
//
// Same ownership rule as tile2d.h, per axis: tile (tx,ty,tz) owns the grid
// vertices [ox, ox+tw) x [oy, oy+th) x [oz, oz+td), and a doubled-lattice cell
// c is owned by the tile owning vertex floor(c/2). Every cell lies in exactly
// one lower star, so every cell is written exactly once.
//
// Host-side bookkeeping (tile origins, boxes, output offsets) is int64. The
// device kernel sees only int32 tile-local quantities; per-axis coordinates
// must be below 2^30 (checked by the drivers), which bounds the volume far
// beyond anything that fits in memory.
//
// The tile size T may differ per axis (slab tiles, one tile for in-core runs).
// Tile size, residency and tiles per launch are scheduling choices and never
// change the output: tiled == in-core byte for byte.

#include <cstdint>
#include <cstring>
#include <functional>

namespace GInt {
namespace gpu {

struct TileGrid3D {
    int64_t X = 0, Y = 0, Z = 0;   // volume size in vertices
    int Tx = 0, Ty = 0, Tz = 0;    // nominal tile size per axis (vertices)
    int ntx = 0, nty = 0, ntz = 0; // tiles per axis
    int64_t numTiles = 0;
};

struct TileDesc3D {
    int tx = 0, ty = 0, tz = 0;
    int64_t ox = 0, oy = 0, oz = 0;  // first owned vertex
    int tw = 0, th = 0, td = 0;      // owned vertices per axis
};

// [x0, x0+w) x [y0, y0+h) x [z0, z0+d) in some lattice.
struct Box3D {
    int64_t x0 = 0, y0 = 0, z0 = 0;
    int w = 0, h = 0, d = 0;
    int64_t Count() const { return (int64_t)w * h * d; }
};

inline TileGrid3D MakeTileGrid3D(int64_t X, int64_t Y, int64_t Z, int Tx, int Ty, int Tz) {
    TileGrid3D g;
    g.X = X;
    g.Y = Y;
    g.Z = Z;
    // a tile never needs to be larger than the volume
    g.Tx = Tx < 1 ? 1 : (Tx > X && X > 0 ? (int)X : Tx);
    g.Ty = Ty < 1 ? 1 : (Ty > Y && Y > 0 ? (int)Y : Ty);
    g.Tz = Tz < 1 ? 1 : (Tz > Z && Z > 0 ? (int)Z : Tz);
    g.ntx = (int)((X + g.Tx - 1) / g.Tx);
    g.nty = (int)((Y + g.Ty - 1) / g.Ty);
    g.ntz = (int)((Z + g.Tz - 1) / g.Tz);
    g.numTiles = (int64_t)g.ntx * g.nty * g.ntz;
    return g;
}
inline TileGrid3D MakeTileGrid3D(int64_t X, int64_t Y, int64_t Z, int T) {
    return MakeTileGrid3D(X, Y, Z, T, T, T);
}

// tileId = tx + ntx*(ty + nty*tz)
inline TileDesc3D TileAt3D(const TileGrid3D& g, int64_t tileId) {
    TileDesc3D d;
    d.tx = (int)(tileId % g.ntx);
    d.ty = (int)((tileId / g.ntx) % g.nty);
    d.tz = (int)(tileId / ((int64_t)g.ntx * g.nty));
    d.ox = (int64_t)d.tx * g.Tx;
    d.oy = (int64_t)d.ty * g.Ty;
    d.oz = (int64_t)d.tz * g.Tz;
    d.tw = (int)((g.X - d.ox) < g.Tx ? (g.X - d.ox) : g.Tx);
    d.th = (int)((g.Y - d.oy) < g.Ty ? (g.Y - d.oy) : g.Ty);
    d.td = (int)((g.Z - d.oz) < g.Tz ? (g.Z - d.oz) : g.Tz);
    return d;
}

namespace tile3d_detail {
inline void Clip(int64_t lo, int64_t hi, int64_t n, int64_t& x0, int& w) {
    if (lo < 0) lo = 0;
    if (hi > n) hi = n;
    x0 = lo;
    w = (int)(hi - lo);
}
} // namespace tile3d_detail

// Owned doubled-lattice cell box: [2o, 2(o+t)) per axis clipped to 2N-1.
inline Box3D OwnedCellBox3D(const TileGrid3D& g, const TileDesc3D& d) {
    Box3D b;
    tile3d_detail::Clip(2 * d.ox, 2 * (d.ox + d.tw), 2 * g.X - 1, b.x0, b.w);
    tile3d_detail::Clip(2 * d.oy, 2 * (d.oy + d.th), 2 * g.Y - 1, b.y0, b.h);
    tile3d_detail::Clip(2 * d.oz, 2 * (d.oz + d.td), 2 * g.Z - 1, b.z0, b.d);
    return b;
}

// Vertices whose lower stars can contain an owned cell: an owned cell's
// coordinate c lies in [2o, 2(o+t)), so its vertices (floor(c/2), ceil(c/2))
// lie in [o, o+t]. One vertex beyond the owned box on the high side only.
inline Box3D ComputedVertexBox3D(const TileGrid3D& g, const TileDesc3D& d) {
    Box3D b;
    tile3d_detail::Clip(d.ox, d.ox + d.tw + 1, g.X, b.x0, b.w);
    tile3d_detail::Clip(d.oy, d.oy + d.th + 1, g.Y, b.y0, b.h);
    tile3d_detail::Clip(d.oz, d.oz + d.td + 1, g.Z, b.z0, b.d);
    return b;
}

// Function values the kernel reads: computed vertices + their 1-ring.
inline Box3D StagedVertexBox3D(const TileGrid3D& g, const TileDesc3D& d) {
    const Box3D c = ComputedVertexBox3D(g, d);
    Box3D b;
    tile3d_detail::Clip(c.x0 - 1, c.x0 + c.w + 1, g.X, b.x0, b.w);
    tile3d_detail::Clip(c.y0 - 1, c.y0 + c.h + 1, g.Y, b.y0, b.h);
    tile3d_detail::Clip(c.z0 - 1, c.z0 + c.d + 1, g.Z, b.z0, b.d);
    return b;
}

// Cells of the lower stars of the computed vertices (the restriction-class
// bytes the constrained kernel reads): the last computed vertex v = c0+cw-1
// reaches cell 2v+1, so [2c0-1, 2(c0+cw)) clipped.
inline Box3D ClassCellBox3D(const TileGrid3D& g, const TileDesc3D& d) {
    const Box3D c = ComputedVertexBox3D(g, d);
    Box3D b;
    tile3d_detail::Clip(2 * c.x0 - 1, 2 * (c.x0 + c.w), 2 * g.X - 1, b.x0, b.w);
    tile3d_detail::Clip(2 * c.y0 - 1, 2 * (c.y0 + c.h), 2 * g.Y - 1, b.y0, b.h);
    tile3d_detail::Clip(2 * c.z0 - 1, 2 * (c.z0 + c.d), 2 * g.Z - 1, b.z0, b.d);
    return b;
}

// Callback I/O: a source fills a box of vertex values, x fastest, with the
// given row pitch and slice pitch (elements); a sink receives a box of cell
// bytes. A byte source supplies a box of cell bytes (restriction classes).
// Boxes are always inside the domain. Return false to abort.
typedef std::function<bool(const Box3D& box, float* dst, int64_t pitch, int64_t slicePitch)>
    TileSource3D;
typedef std::function<bool(const Box3D& box, uint8_t* dst, int64_t pitch, int64_t slicePitch)>
    ByteSource3D;
typedef std::function<bool(const Box3D& box, const uint8_t* src, int64_t pitch,
                           int64_t slicePitch)>
    ByteSink3D;

// In-core adapters over host arrays (X*Y*Z floats; (2X-1)(2Y-1)(2Z-1) bytes).
inline TileSource3D HostVolumeSource3D(const float* values, int64_t X, int64_t Y) {
    return [values, X, Y](const Box3D& b, float* dst, int64_t pitch, int64_t spitch) {
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch,
                            values + ((b.z0 + z) * Y + (b.y0 + y)) * X + b.x0,
                            sizeof(float) * (size_t)b.w);
        return true;
    };
}
inline ByteSource3D HostCellSource3D(const uint8_t* cells, int64_t X, int64_t Y) {
    const int64_t W = 2 * X - 1, H = 2 * Y - 1;
    return [cells, W, H](const Box3D& b, uint8_t* dst, int64_t pitch, int64_t spitch) {
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch,
                            cells + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0, (size_t)b.w);
        return true;
    };
}
inline ByteSink3D HostGradSink3D(uint8_t* grad, int64_t X, int64_t Y) {
    const int64_t W = 2 * X - 1, H = 2 * Y - 1;
    return [grad, W, H](const Box3D& b, const uint8_t* src, int64_t pitch, int64_t spitch) {
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(grad + ((b.z0 + z) * H + (b.y0 + y)) * W + b.x0,
                            src + z * spitch + y * pitch, (size_t)b.w);
        return true;
    };
}

struct TiledOptions3D {
    int slots = 3;           // resident tile batches (pinned + device buffers)
    int tilesPerLaunch = 1;  // tiles staged into one launch
    int block_shape = 0;     // 0 = 32x4x2, 1 = 16x8x2, 2 = 8x8x4
    int hostThreads = 0;     // threads for the in-core adapters' copies (0 = auto)
};

struct TileStats3D {
    int64_t tiles = 0;
    int64_t launches = 0;
    int64_t bytesH2D = 0;
    int64_t bytesD2H = 0;
    float h2d_ms = 0.f;      // per-batch event windows, summed (overlapped
    float kernel_ms = 0.f;   //  legs may add up to more than wall_ms)
    float d2h_ms = 0.f;
    float source_ms = 0.f;   // host time inside the source callbacks
    float sink_ms = 0.f;     // host time inside the sink callbacks
    float wall_ms = 0.f;
    int64_t deviceFreeMin = -1;        // cudaMemGetInfo sampled once per batch
    int64_t deviceFreeMax = -1;
    int64_t deviceBytesAllocated = 0;  // this call's own device allocations
    int64_t hostPinnedBytes = 0;       // this call's own pinned allocations
};

} // namespace gpu
} // namespace GInt

#endif
