// Stage A4: tiled driver for the numeric separatrix lines (lines2d.cuh).
//
// Two phases, each a set of rounds over tiles. Round 0 of a phase visits
// every tile: a seed kernel replicates the CPU's depth-10 BFS from each
// critical edge the tile owns (descending: facet vertices along the V-path
// down; ascending: cofacet quads along the V-path up) into a device seed
// list with a seniority key per seed (saddle value, seed index, saddle
// hash); a walk kernel then runs one thread per seed. Walks that leave the
// owned box record their key (atomicMin) on the tile's outer ring for the
// exit target; the host maps records to the owning neighbour tile's pending
// list, and later rounds visit only tiles with pending resumes.
//
// Phase 1 claims ownership of every lattice point a walk passes (host-side
// owner arrays, min-merged per tile per round); phase 2 marks the class
// bits, each walk stopping at the first point a more senior walk owns.
// Marks go to the sink per tile per round as a box one cell wider than the
// owned cell box, to be OR'd by the sink. The device footprint is slots x
// tilesPerLaunch tiles; the owner keys (8 B per lattice point, two
// lattices) are the one allocation proportional to the image in this stage
// (a tile-streamed state; in-core here).

#include "cls/cls2d_internal.cuh"
#include "common/gpu_transfer.h"
#include "lines2d.cuh"
#include "tile2d.h"

#include <chrono>
#include <cstdio>
#include <vector>

namespace GInt {
namespace gpu {

namespace {

using namespace lines2d;

__constant__ LineParams c_lp;

struct SeedEntry { int32_t x, y; Key key; };   // global lattice coords + seniority key

struct LineWork {
    LineTile t;
    SeedEntry* seedsV; int* nSeedsV; int capV;
    SeedEntry* seedsQ; int* nSeedsQ; int capQ;
    int* walkEnds;   // [5]
    int* dropped;    // [1]
    int doAsc, doDsc;
    int ocx0, ocy0, ocw, och;   // owned cell box (edges to seed from)
};

__device__ __forceinline__ void pushSeed(SeedEntry* list, int* n, int cap, int x, int y, Key key, int* dropped) {
    const int k = atomicAdd(n, 1);
    if (k < cap) { list[k].x = x; list[k].y = y; list[k].key = key; }
    else atomicAdd(dropped, 1);
}

// monotone float -> uint32
__device__ __forceinline__ unsigned orderedFloat(float f) {
    unsigned u = __float_as_uint(f);
    return (u & 0x80000000u) ? ~u : (u | 0x80000000u);
}
__device__ __forceinline__ unsigned hash32(unsigned a) {
    a ^= a >> 16; a *= 0x7feb352du; a ^= a >> 15; a *= 0x846ca68bu; a ^= a >> 16;
    return a;
}

// The CPU BFS (gi_discrete_gradient_computer.h:237-271 / :382-416): from
// the saddle edge, pop up to bfsDepth cells; for each facet (DSC) / cofacet
// (ASC) that is critical -> seed; else -> seed, mark it and its pair
// visited, push the pair. Replicated on the doubled lattice with the pair
// codes. Seed k of saddle s gets key (value order, k, hash(s)).
__global__ void seed_kernel(const LineWork* __restrict__ works) {
    const LineWork w = works[blockIdx.y];
    const LineTile& t = w.t;
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= w.ocw * w.och) return;
    const int cx = w.ocx0 + i % w.ocw, cy = w.ocy0 + i / w.ocw;
    if (((cx ^ cy) & 1) == 0) return;          // not an edge
    if (gcode(t, cx, cy) != 7u) return;        // not critical
    // saddle value = max of its two vertices (lower-star function)
    float sv;
    if (cx & 1) { const float a = fval(t, (cx - 1) / 2, cy / 2), b = fval(t, (cx + 1) / 2, cy / 2); sv = a > b ? a : b; }
    else        { const float a = fval(t, cx / 2, (cy - 1) / 2), b = fval(t, cx / 2, (cy + 1) / 2); sv = a > b ? a : b; }
    const unsigned ov = orderedFloat(sv);
    const unsigned h = hash32((unsigned)(cx + cy * (2 * c_lp.X - 1))) & 0x3FFFFFEu;
    for (int cls = 0; cls < 2; cls++) {
        const bool asc = cls == 0;
        if (asc && !w.doAsc) continue;
        if (!asc && !w.doDsc) continue;
        // ASC processes the highest saddle first (smaller key = more senior)
        const Key hi = asc ? (Key)(0xFFFFFFFFu - ov) : (Key)ov;
        int qx[32], qy[32]; int qh = 0, qt = 0;
        int vx[64], vy[64]; int nv = 0;
        auto visited = [&](int x, int y) { for (int k = 0; k < nv; k++) if (vx[k] == x && vy[k] == y) return true; return false; };
        auto visit = [&](int x, int y) { if (nv < 64) { vx[nv] = x; vy[nv] = y; nv++; } };
        qx[qt] = cx; qy[qt] = cy; qt++;
        int counter = c_lp.bfsDepth;
        int seedIdx = 0;
        while (qh < qt && counter >= 0) {
            const int px = qx[qh], py = qy[qh]; qh++;
            visit(px, py);
            const int dimp = (px & 1) + (py & 1);
            int nx[2], ny[2], nn = 0;
            if (dimp == 1) {
                if (!asc) {   // edge -> its two vertices
                    if (px & 1) { nx[0] = px - 1; ny[0] = py; nx[1] = px + 1; ny[1] = py; }
                    else        { nx[0] = px; ny[0] = py - 1; nx[1] = px; ny[1] = py + 1; }
                } else {      // edge -> its (up to) two quads
                    if (px & 1) { nx[0] = px; ny[0] = py - 1; nx[1] = px; ny[1] = py + 1; }
                    else        { nx[0] = px - 1; ny[0] = py; nx[1] = px + 1; ny[1] = py; }
                }
                nn = 2;
            }
            for (int k = 0; k < nn; k++) {
                const int tx = nx[k], ty = ny[k];
                if (tx < 0 || ty < 0 || tx > 2 * c_lp.X - 2 || ty > 2 * c_lp.Y - 2) continue;
                if (!gradInBox(t, tx, ty)) continue;   // BFS clamped to the staged box
                if (visited(tx, ty)) continue;
                const unsigned code = gcode(t, tx, ty);
                const Key key = (hi << 32) | ((Key)(seedIdx < 63 ? seedIdx : 63) << 26) | (Key)h;
                seedIdx++;
                if (asc) pushSeed(w.seedsQ, w.nSeedsQ, w.capQ, (tx - 1) / 2, (ty - 1) / 2, key, w.dropped);
                else pushSeed(w.seedsV, w.nSeedsV, w.capV, tx / 2, ty / 2, key, w.dropped);
                visit(tx, ty);
                if (code == 7u) continue;
                const int ppx = tx + c_lp.dx[code], ppy = ty + c_lp.dy[code];   // its pair (an edge)
                visit(ppx, ppy);
                if (qt < 32) { qx[qt] = ppx; qy[qt] = ppy; qt++; }
            }
            counter--;
        }
    }
    // the saddle cell itself gets both classes it seeds (CPU: set_label(sad_id))
    if (c_lp.phase == 2) {
        if (w.doAsc) markCell(t, cx, cy, kClassAsc);
        if (w.doDsc) markCell(t, cx, cy, kClassDsc);
    }
}

__global__ void walk_kernel(const LineWork* __restrict__ works, int asc) {
    const LineWork w = works[blockIdx.y];
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = asc ? min(*w.nSeedsQ, w.capQ) : min(*w.nSeedsV, w.capV);
    if (i >= n) return;
    const SeedEntry s = asc ? w.seedsQ[i] : w.seedsV[i];
    const int e = asc ? walk<true>(w.t, c_lp, s.x, s.y, s.key) : walk<false>(w.t, c_lp, s.x, s.y, s.key);
    atomicAdd(w.walkEnds + e, 1);
}

__global__ void fill_keys_kernel(Key* p, int n, Key v) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) p[i] = v;
}

#define LINE_CUDA_CHECK(call)                                                   \
    do {                                                                        \
        cudaError_t err__ = (call);                                             \
        if (err__ != cudaSuccess) {                                             \
            std::fprintf(stderr, "CUDA error %s at %s:%d: %s\n",                \
                         cudaGetErrorName(err__), __FILE__, __LINE__,           \
                         cudaGetErrorString(err__));                            \
            return false;                                                       \
        }                                                                       \
    } while (0)

struct Pending { std::vector<SeedEntry> v, q; };

static inline size_t round8(size_t b) { return (b + 7) & ~(size_t)7; }

} // namespace

bool ComputeLineClasses2DTiled(const TileSource2D& values,
                               const ByteSource2D& gradSource,
                               const TileGrid2D& grid,
                               const Dgrad2DTables& tables,
                               const ByteSink2D& sinkOr,
                               const LineOptions2D& options,
                               LineStats2D* stats) {
    if (grid.X < 3 || grid.Y < 3 || grid.T < 1 || grid.numTiles < 1) return false;
    const auto t0 = std::chrono::steady_clock::now();
    LineStats2D st;
    st.tiles = grid.numTiles;
    const int64_t X = grid.X, Y = grid.Y;

    LineParams lp{};
    lp.X = (int)X; lp.Y = (int)Y;
    lp.hopCells = options.hopCellsLine < 1 ? 1 : options.hopCellsLine;
    lp.iterCap = options.iterCap < 1 ? 1 : options.iterCap;
    lp.maxHops = options.maxHops < 1 ? 1 : options.maxHops;
    lp.walkIter = options.walkIter;
    lp.bfsDepth = options.bfsDepth;
    lp.errThresh = options.errorThreshold;
    lp.gradThresh = options.gradThreshold;
    lp.dx[tables.dir_px] = 1;  lp.dy[tables.dir_px] = 0;
    lp.dx[tables.dir_mx] = -1; lp.dy[tables.dir_mx] = 0;
    lp.dx[tables.dir_py] = 0;  lp.dy[tables.dir_py] = 1;
    lp.dx[tables.dir_my] = 0;  lp.dy[tables.dir_my] = -1;

    const int T = grid.T;
    // residency never exceeds the tile count (a one-tile image must not pay
    // for 12 resident tiles of buffers)
    int tpl = options.tilesPerLaunch < 1 ? 1 : options.tilesPerLaunch > 64 ? 64 : options.tilesPerLaunch;
    if ((int64_t)tpl > grid.numTiles) tpl = (int)grid.numTiles;
    int nslots = options.slots < 1 ? 1 : options.slots;
    const int64_t batches = (grid.numTiles + tpl - 1) / tpl;
    if ((int64_t)nslots > batches) nslots = (int)batches;
    // BFS seeds per saddle per class: each of the bfsDepth+1 pops pushes at
    // most 2 -> an exact per-tile bound from the critical edges the host
    // counts while staging the gradient. Seed capacity must never drop a
    // seed: which seed drops would depend on atomic order, and the result
    // would stop being a function of the gradient.
    const int seedsPerSaddle = 2 * (options.bfsDepth + 1) + 2;
    // pending resumes per tile are bounded by the neighbours' outer rings
    const int pendCap = 8 * (T + 2) + 8;
    const size_t valsB = (size_t)(T + 4) * (T + 4) * sizeof(float);
    const size_t gradB = round8((size_t)(2 * T + 8) * (2 * T + 8));
    const size_t ownB = (size_t)(T + 2) * (T + 2) * sizeof(Key);   // one lattice box + ring
    const size_t pendB = (size_t)pendCap * 2 * sizeof(SeedEntry);
    const size_t capIn = (valsB + gradB + 2 * ownB + pendB) * (size_t)tpl;
    const int RoutMax = RingCount2D(T + 2, T + 2);
    const size_t clsB = round8((size_t)(2 * T + 2) * (2 * T + 2));
    const size_t flagsB = (size_t)2 * RoutMax * sizeof(Key);
    const size_t capOut = (clsB + flagsB + 2 * ownB) * (size_t)tpl;

    // host owner arrays (vertex lattice X*Y, quad lattice (X-1)*(Y-1))
    const int64_t QX = X - 1, QY = Y - 1;
    std::vector<Key> ownV((size_t)(X * Y), LINES_NO_KEY), ownQ((size_t)(QX * QY), LINES_NO_KEY);
    st.hostOwnerBytes = (int64_t)((ownV.size() + ownQ.size()) * sizeof(Key));

    TransferRing ring;
    if (!ring.Create(nslots, capIn, capOut)) return false;
    LineWork* d_works = nullptr;
    // seed lists: one grow-only device buffer per slot, sized per batch to
    // the exact bound (retire waits for the slot, so growing is safe)
    std::vector<SeedEntry*> d_seeds((size_t)nslots, nullptr);
    std::vector<size_t> seedAlloc((size_t)nslots, 0);
    int* d_counts = nullptr;
    int* d_stats = nullptr;
    LINE_CUDA_CHECK(cudaMalloc(&d_works, (size_t)nslots * tpl * sizeof(LineWork)));
    LINE_CUDA_CHECK(cudaMalloc(&d_counts, (size_t)nslots * tpl * 2 * sizeof(int)));
    LINE_CUDA_CHECK(cudaMalloc(&d_stats, 8 * sizeof(int)));
    LINE_CUDA_CHECK(cudaMemset(d_stats, 0, 8 * sizeof(int)));
    st.deviceBytesAllocated = (int64_t)ring.DeviceBytes() +
                              (int64_t)((size_t)nslots * tpl * (sizeof(LineWork) + 2 * sizeof(int)));
    st.hostPinnedBytes = (int64_t)ring.PinnedBytes();
    int64_t seedBytesMax = 0;

    bool ok = true;
    for (int phase = 1; phase <= 2 && ok; phase++) {
        lp.phase = phase;
        LINE_CUDA_CHECK(cudaMemcpyToSymbol(c_lp, &lp, sizeof(lp)));
        LINE_CUDA_CHECK(cudaMemset(d_stats, 0, 8 * sizeof(int)));
        std::vector<Pending> pending((size_t)grid.numTiles);
        std::vector<int64_t> todo;
        for (int64_t i = 0; i < grid.numTiles; i++) todo.push_back(i);

        for (int round = 0; round < options.maxRounds && !todo.empty() && ok; round++) {
            st.phaseRounds[phase - 1]++;
            std::vector<Pending> nextPending((size_t)grid.numTiles);
            std::vector<std::vector<LineWork>> works((size_t)nslots);
            std::vector<int> counts((size_t)nslots, 0);
            std::vector<std::vector<size_t>> seedOff((size_t)nslots);
            std::vector<std::vector<int>> seedCounts((size_t)nslots);

            auto retire = [&](int slot) -> bool {
                if (counts[(size_t)slot] == 0) return true;
                if (!ring.Retire(slot, nullptr, nullptr, nullptr)) return false;
                const unsigned char* h_out = (const unsigned char*)ring.Slot(slot).h_out;
                size_t off = 0;
                for (int i = 0; i < counts[(size_t)slot]; i++) {
                    const LineWork& w = works[(size_t)slot][(size_t)i];
                    const LineTile& lt = w.t;
                    // class box (phase 2 only)
                    if (phase == 2 && !sinkOr(lt.ccx0, lt.ccy0, lt.cw, lt.ch, h_out + off, lt.cpitch)) return false;
                    off += clsB;
                    // resume keys -> pending seeds of the owning tiles
                    const Key* fl = (const Key*)(h_out + off);
                    off += flagsB;
                    for (int cls = 0; cls < 2; cls++) {
                        const bool asc = cls == 0;
                        const int bw = asc ? lt.qw : lt.tw, bh = asc ? lt.qh : lt.th;
                        const int Rout = RingCount2D(bw + 2, bh + 2);
                        const Key* f = fl + (asc ? RoutMax : 0);
                        const int64_t NX = asc ? QX : X, NY = asc ? QY : Y;
                        for (int r = 0; r < Rout; r++) {
                            if (f[r] == LINES_NO_KEY) continue;
                            int x, y;
                            RingToLocal2D(r, bw + 2, bh + 2, x, y);
                            const int64_t gx = lt.ox + x - 1, gy = lt.oy + y - 1;
                            if (gx < 0 || gy < 0 || gx >= NX || gy >= NY) continue;
                            const int64_t nt = (gx / T) + (gy / T) * (int64_t)grid.ntx;
                            SeedEntry e; e.x = (int32_t)gx; e.y = (int32_t)gy; e.key = f[r];
                            if (asc) nextPending[(size_t)nt].q.push_back(e);
                            else nextPending[(size_t)nt].v.push_back(e);
                            st.resumes++;
                        }
                    }
                    // owner boxes: min-merge into the host arrays (phase 1)
                    const Key* ov = (const Key*)(h_out + off);
                    off += ownB;
                    const Key* oq = (const Key*)(h_out + off);
                    off += ownB;
                    if (phase == 1) {
                        for (int r = 0; r < lt.ovh; r++)
                            for (int c = 0; c < lt.ovw; c++) {
                                Key& h = ownV[(size_t)((lt.ovy0 + r) * X + lt.ovx0 + c)];
                                const Key d = ov[(size_t)r * lt.ovpitch + c];
                                if (d < h) h = d;
                            }
                        for (int r = 0; r < lt.oqh; r++)
                            for (int c = 0; c < lt.oqw; c++) {
                                Key& h = ownQ[(size_t)((lt.oqy0 + r) * QX + lt.oqx0 + c)];
                                const Key d = oq[(size_t)r * lt.oqpitch + c];
                                if (d < h) h = d;
                            }
                    }
                }
                counts[(size_t)slot] = 0;
                return true;
            };

            size_t next = 0;
            int slot = 0;
            while (next < todo.size() && ok) {
                if (!retire(slot)) { ok = false; break; }
                const int cnt = (int)((todo.size() - next) < (size_t)tpl ? (todo.size() - next) : (size_t)tpl);
                std::vector<LineWork>& wk = works[(size_t)slot];
                wk.assign((size_t)cnt, LineWork{});
                seedOff[(size_t)slot].assign((size_t)cnt, 0);
                seedCounts[(size_t)slot].assign((size_t)cnt * 2, 0);
                TransferSlot& s = ring.Slot(slot);
                unsigned char* h_in = (unsigned char*)s.h_in;
                const unsigned char* d_in = (const unsigned char*)s.d_in;
                size_t inOff = 0, outOff = 0;
                int maxCells = 0, maxSeeds = 1;
                std::vector<size_t> ownInOff((size_t)cnt * 2, 0);
                for (int i = 0; i < cnt && ok; i++) {
                    const int64_t tid = todo[next + i];
                    const TileDesc2D d = TileAt(grid, tid);
                    LineWork& w = wk[(size_t)i];
                    LineTile& lt = w.t;
                    lt.ox = (int)d.ox; lt.oy = (int)d.oy; lt.tw = d.tw; lt.th = d.th;
                    const Box2D qb = LatticeBox2D(d, QX, QY);
                    lt.qw = qb.w; lt.qh = qb.h;
                    // values, halo 2
                    {
                        const int64_t x0 = d.ox - 2 < 0 ? 0 : d.ox - 2, y0 = d.oy - 2 < 0 ? 0 : d.oy - 2;
                        const int64_t x1 = d.ox + d.tw + 2 > X ? X : d.ox + d.tw + 2;
                        const int64_t y1 = d.oy + d.th + 2 > Y ? Y : d.oy + d.th + 2;
                        const int wv = (int)(x1 - x0), hv = (int)(y1 - y0);
                        if (!values(x0, y0, wv, hv, (float*)(h_in + inOff), wv)) { ok = false; break; }
                        lt.vals = (const float*)(d_in + inOff);
                        lt.vpitch = wv; lt.sx0 = (int)x0; lt.sy0 = (int)y0; lt.svx = wv; lt.svy = hv;
                        inOff += (size_t)wv * hv * sizeof(float);
                    }
                    // gradient bytes: owned cell box + 4 cells halo
                    const Box2D ob = OwnedCellBox2D(grid, d);
                    {
                        const int64_t x0 = ob.x0 - 4 < 0 ? 0 : ob.x0 - 4, y0 = ob.y0 - 4 < 0 ? 0 : ob.y0 - 4;
                        const int64_t x1 = ob.x0 + ob.w + 4 > 2 * X - 1 ? 2 * X - 1 : ob.x0 + ob.w + 4;
                        const int64_t y1 = ob.y0 + ob.h + 4 > 2 * Y - 1 ? 2 * Y - 1 : ob.y0 + ob.h + 4;
                        const int wg = (int)(x1 - x0), hg = (int)(y1 - y0);
                        if (!gradSource(x0, y0, wg, hg, h_in + inOff, wg)) { ok = false; break; }
                        lt.grad = d_in + inOff;
                        lt.gpitch = wg; lt.gcx0 = (int)x0; lt.gcy0 = (int)y0; lt.gw = wg; lt.gh = hg;
                        w.ocx0 = (int)ob.x0; w.ocy0 = (int)ob.y0; w.ocw = ob.w; w.och = ob.h;
                        // exact seed bound: critical edges of the owned cell box
                        int crit = 0;
                        if (round == 0) {
                            const unsigned char* gb = h_in + inOff;
                            for (int r = 0; r < ob.h; r++) {
                                const int64_t cy = ob.y0 + r;
                                const unsigned char* row = gb + (size_t)(cy - y0) * wg + (ob.x0 - x0);
                                for (int c = 0; c < ob.w; c++) {
                                    const int64_t cx = ob.x0 + c;
                                    if (((cx ^ cy) & 1) && ((row[c] >> 2) & 7) == 7) crit++;
                                }
                            }
                        }
                        w.capV = w.capQ = crit * seedsPerSaddle;
                        inOff += round8((size_t)wg * hg);
                    }
                    // owner boxes (+1 ring), staged from the host arrays
                    {
                        const int64_t x0 = d.ox - 1 < 0 ? 0 : d.ox - 1, y0 = d.oy - 1 < 0 ? 0 : d.oy - 1;
                        const int64_t x1 = d.ox + d.tw + 1 > X ? X : d.ox + d.tw + 1;
                        const int64_t y1 = d.oy + d.th + 1 > Y ? Y : d.oy + d.th + 1;
                        lt.ovx0 = (int)x0; lt.ovy0 = (int)y0; lt.ovw = (int)(x1 - x0); lt.ovh = (int)(y1 - y0); lt.ovpitch = lt.ovw;
                        Key* hk = (Key*)(h_in + inOff);
                        for (int r = 0; r < lt.ovh; r++)
                            for (int c = 0; c < lt.ovw; c++) hk[(size_t)r * lt.ovw + c] = ownV[(size_t)((y0 + r) * X + x0 + c)];
                        ownInOff[(size_t)i * 2] = inOff;
                        inOff += ownB;
                        const int64_t qx0 = qb.x0 - 1 < 0 ? 0 : qb.x0 - 1, qy0 = qb.y0 - 1 < 0 ? 0 : qb.y0 - 1;
                        const int64_t qx1 = qb.x0 + qb.w + 1 > QX ? QX : qb.x0 + qb.w + 1;
                        const int64_t qy1 = qb.y0 + qb.h + 1 > QY ? QY : qb.y0 + qb.h + 1;
                        lt.oqx0 = (int)qx0; lt.oqy0 = (int)qy0; lt.oqw = (int)(qx1 > qx0 ? qx1 - qx0 : 0); lt.oqh = (int)(qy1 > qy0 ? qy1 - qy0 : 0); lt.oqpitch = lt.oqw;
                        hk = (Key*)(h_in + inOff);
                        for (int r = 0; r < lt.oqh; r++)
                            for (int c = 0; c < lt.oqw; c++) hk[(size_t)r * lt.oqw + c] = ownQ[(size_t)((qy0 + r) * QX + qx0 + c)];
                        ownInOff[(size_t)i * 2 + 1] = inOff;
                        inOff += ownB;
                    }
                    // outputs: class box, resume keys, owner boxes (working copies)
                    {
                        const int64_t cx0 = ob.x0 - 1 < 0 ? 0 : ob.x0 - 1, cy0 = ob.y0 - 1 < 0 ? 0 : ob.y0 - 1;
                        const int64_t cx1 = ob.x0 + ob.w + 1 > 2 * X - 1 ? 2 * X - 1 : ob.x0 + ob.w + 1;
                        const int64_t cy1 = ob.y0 + ob.h + 1 > 2 * Y - 1 ? 2 * Y - 1 : ob.y0 + ob.h + 1;
                        lt.ccx0 = (int)cx0; lt.ccy0 = (int)cy0; lt.cw = (int)(cx1 - cx0); lt.ch = (int)(cy1 - cy0);
                        lt.cpitch = lt.cw;
                        lt.cls = (unsigned char*)s.d_out + outOff;
                        outOff += clsB;
                        lt.resumeV = (Key*)((unsigned char*)s.d_out + outOff);
                        lt.resumeQ = lt.resumeV + RoutMax;
                        lt.RoutV = RingCount2D(d.tw + 2, d.th + 2);
                        lt.RoutQ = RingCount2D(qb.w + 2, qb.h + 2);
                        outOff += flagsB;
                        lt.ownV = (Key*)((unsigned char*)s.d_out + outOff);
                        outOff += ownB;
                        lt.ownQ = (Key*)((unsigned char*)s.d_out + outOff);
                        outOff += ownB;
                    }
                    // seed lists: pending resumes staged through the pinned input
                    {
                        w.nSeedsV = d_counts + ((size_t)slot * tpl + (size_t)i) * 2;
                        w.nSeedsQ = w.nSeedsV + 1;
                        const Pending& p = pending[(size_t)tid];
                        int nv = (int)p.v.size(), nq = (int)p.q.size();
                        if (nv > pendCap) { st.droppedSeeds += nv - pendCap; nv = pendCap; }
                        if (nq > pendCap) { st.droppedSeeds += nq - pendCap; nq = pendCap; }
                        // capacity = BFS bound (round 0) + pending resumes
                        w.capV += nv; w.capQ += nq;
                        if (w.capV < 1) w.capV = 1;
                        if (w.capQ < 1) w.capQ = 1;
                        SeedEntry* hv = (SeedEntry*)(h_in + inOff);
                        for (int k = 0; k < nv; k++) hv[k] = p.v[(size_t)k];
                        for (int k = 0; k < nq; k++) hv[nv + k] = p.q[(size_t)k];
                        seedOff[(size_t)slot][(size_t)i] = inOff;
                        seedCounts[(size_t)slot][(size_t)i * 2] = nv;
                        seedCounts[(size_t)slot][(size_t)i * 2 + 1] = nq;
                        inOff += (size_t)(nv + nq) * sizeof(SeedEntry);
                        w.walkEnds = d_stats; w.dropped = d_stats + 5;
                        w.doAsc = options.ascending ? 1 : 0; w.doDsc = options.descending ? 1 : 0;
                        const int ns = nv > nq ? nv : nq;
                        if (ns > maxSeeds) maxSeeds = ns;
                    }
                    if (w.ocw * w.och > maxCells) maxCells = w.ocw * w.och;
                }
                if (!ok) break;
                counts[(size_t)slot] = cnt;

                // seed buffer of this slot: grow to the batch's exact bound
                {
                    size_t need = 0;
                    for (int i = 0; i < cnt; i++) need += (size_t)wk[(size_t)i].capV + (size_t)wk[(size_t)i].capQ;
                    if (need > seedAlloc[(size_t)slot]) {
                        if (d_seeds[(size_t)slot]) cudaFree(d_seeds[(size_t)slot]);
                        d_seeds[(size_t)slot] = nullptr;
                        LINE_CUDA_CHECK(cudaMalloc(&d_seeds[(size_t)slot], need * sizeof(SeedEntry)));
                        seedAlloc[(size_t)slot] = need;
                        if ((int64_t)(need * sizeof(SeedEntry)) > seedBytesMax) seedBytesMax = (int64_t)(need * sizeof(SeedEntry));
                    }
                    size_t off = 0;
                    for (int i = 0; i < cnt; i++) {
                        LineWork& w = wk[(size_t)i];
                        w.seedsV = d_seeds[(size_t)slot] + off; off += (size_t)w.capV;
                        w.seedsQ = d_seeds[(size_t)slot] + off; off += (size_t)w.capQ;
                        const int c = w.capV > w.capQ ? w.capV : w.capQ;
                        if (c > maxSeeds) maxSeeds = c;
                    }
                }

                LineWork* dw = d_works + (size_t)slot * tpl;
                LINE_CUDA_CHECK(cudaMemcpyAsync(dw, wk.data(), (size_t)cnt * sizeof(LineWork),
                                                cudaMemcpyHostToDevice, ring.StreamH2D()));
                if (!ring.BeginUpload(slot, inOff)) { ok = false; break; }
                cudaStream_t ks = ring.StreamKernel();
                LINE_CUDA_CHECK(cudaMemsetAsync(s.d_out, 0, outOff, ks));
                const int TPB = 128;
                for (int i = 0; i < cnt; i++) {
                    const LineWork& w = wk[(size_t)i];
                    // resume keys start at NO_KEY; owner working copies from the staged input
                    fill_keys_kernel<<<(2 * RoutMax + TPB - 1) / TPB, TPB, 0, ks>>>(w.t.resumeV, 2 * RoutMax, LINES_NO_KEY);
                    LINE_CUDA_CHECK(cudaMemcpyAsync(w.t.ownV, d_in + ownInOff[(size_t)i * 2], ownB, cudaMemcpyDeviceToDevice, ks));
                    LINE_CUDA_CHECK(cudaMemcpyAsync(w.t.ownQ, d_in + ownInOff[(size_t)i * 2 + 1], ownB, cudaMemcpyDeviceToDevice, ks));
                    const int nv = seedCounts[(size_t)slot][(size_t)i * 2], nq = seedCounts[(size_t)slot][(size_t)i * 2 + 1];
                    const size_t so = seedOff[(size_t)slot][(size_t)i];
                    if (nv) LINE_CUDA_CHECK(cudaMemcpyAsync(w.seedsV, d_in + so, (size_t)nv * sizeof(SeedEntry), cudaMemcpyDeviceToDevice, ks));
                    if (nq) LINE_CUDA_CHECK(cudaMemcpyAsync(w.seedsQ, d_in + so + (size_t)nv * sizeof(SeedEntry), (size_t)nq * sizeof(SeedEntry), cudaMemcpyDeviceToDevice, ks));
                    LINE_CUDA_CHECK(cudaMemcpyAsync(w.nSeedsV, &seedCounts[(size_t)slot][(size_t)i * 2], 2 * sizeof(int), cudaMemcpyHostToDevice, ks));
                }
                LINE_CUDA_CHECK(cudaGetLastError());
                if (round == 0) {
                    dim3 g1((unsigned)((maxCells + TPB - 1) / TPB), (unsigned)cnt);
                    seed_kernel<<<g1, TPB, 0, ks>>>(dw);
                    LINE_CUDA_CHECK(cudaGetLastError());
                    // counts live on the device; the grid covers the exact
                    // bound and threads beyond the count return
                }
                dim3 g2((unsigned)((maxSeeds + TPB - 1) / TPB), (unsigned)cnt);
                if (options.descending) { walk_kernel<<<g2, TPB, 0, ks>>>(dw, 0); LINE_CUDA_CHECK(cudaGetLastError()); }
                if (options.ascending) { walk_kernel<<<g2, TPB, 0, ks>>>(dw, 1); LINE_CUDA_CHECK(cudaGetLastError()); }
                if (!ring.EndComputeAndDownload(slot, outOff)) { ok = false; break; }
                next += (size_t)cnt;
                slot = (slot + 1) % nslots;
            }
            for (int i = 0; i < nslots && ok; i++)
                if (!retire((slot + i) % nslots)) ok = false;
            if (!ok) break;
            pending.swap(nextPending);
            todo.clear();
            for (int64_t i = 0; i < grid.numTiles; i++)
                if (!pending[(size_t)i].v.empty() || !pending[(size_t)i].q.empty()) todo.push_back(i);
        }
        if (ok && !todo.empty())
            std::fprintf(stderr, "ComputeLineClasses2DTiled: maxRounds reached in phase %d with %zu tiles pending\n",
                         phase, todo.size());
        int h[8];
        if (cudaMemcpy(h, d_stats, sizeof(h), cudaMemcpyDeviceToHost) == cudaSuccess) {
            if (phase == 2) for (int k = 0; k < 5; k++) st.walkEnds[k] = h[k];
            st.droppedSeeds += h[5];
        }
    }
    st.rounds = st.phaseRounds[0] + st.phaseRounds[1];
    st.deviceBytesAllocated += seedBytesMax * nslots;
    cudaFree(d_works);
    for (int i = 0; i < nslots; i++) if (d_seeds[(size_t)i]) cudaFree(d_seeds[(size_t)i]);
    cudaFree(d_counts);
    cudaFree(d_stats);
    ring.Destroy();
    st.wall_ms = (float)std::chrono::duration<double, std::milli>(
                     std::chrono::steady_clock::now() - t0).count();
    if (stats) *stats = st;
    return ok;
}

} // namespace gpu
} // namespace GInt
