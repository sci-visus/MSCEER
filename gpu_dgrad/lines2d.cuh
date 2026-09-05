#ifndef DGRAD_LINES2D_CUH
#define DGRAD_LINES2D_CUH

// Stage A4: device line walks - the numeric 1-manifolds (separatrices) of
// the 2D accurate gradient, digitized into the doubled lattice as class
// marks (1 = ascending line: quad-edge-quad chain toward a maximum, 2 =
// descending line: vertex-edge-vertex chain toward a minimum), the cells the
// CPU's DigitizingNumericStreamlineIntegrator2dASC/DSC paint.
//
// A walk is a chain of numeric hops (numeric2d.cuh semantics: the CPU's
// adaptive 2-step Euler on the bilinear gradient) of hopCellsLine voxels
// each, snapping to a lattice point after every hop (descending: the
// nearest vertex of the end point, accepted only if strictly below the
// start under the (value, index) order, else the discrete V-path step;
// ascending: the containing quad under the quad order); each accepted Euler
// segment is digitized with the CPU's plane-crossing rules.
//
// Thinness the CPU gets from HIT_PREASSIGNED (a line stops when it touches
// a cell an earlier line painted - earlier in thread timing) is reproduced
// DETERMINISTICALLY in two phases: phase 1 walks every seed and records at
// each lattice point it passes the most SENIOR walk (atomicMin of a 64-bit
// key: saddle value first - the CPU's processing order - then seed index,
// then a saddle hash), phase 2 walks again and marks the path only until
// the first lattice point owned by a more senior walk, marking that point
// too so the junior line joins the senior one. Keys make the result a
// function of the gradient alone.
//
// A walk ends at a critical cell of its lattice, at the domain boundary
// band (where the CPU line integrators stop too), when it merges into a
// senior line (phase 2), on leaving the tile's owned box (a resume record -
// the most senior key exiting through that target - is set on the tile's
// outer ring; the host re-seeds the neighbouring tile next round), when its
// Euler budget (walkIter, the CPU's 500 per streamline) is spent, or at
// maxHops.

#include "numeric2d.cuh"

#include <cuda_runtime.h>
#include <cstdint>

namespace GInt {
namespace gpu {
namespace lines2d {

enum { kClassAsc = 1, kClassDsc = 2 };
typedef unsigned long long Key;
#define LINES_NO_KEY 0xFFFFFFFFFFFFFFFFull

struct LineTile {
    // staged values, owned vertex box + 2 halo (clipped)
    const float* vals; int vpitch, sx0, sy0, svx, svy;
    // staged gradient bytes: owned cell box + 4 cells halo (clipped)
    const unsigned char* grad; int gpitch, gcx0, gcy0, gw, gh;
    // owned vertex box and owned quad box
    int ox, oy, tw, th;      // vertices
    int qw, qh;              // quads: [ox, ox+qw) x [oy, oy+qh)
    // class marks: owned cell box + 1 cell ring (clipped); byte OR
    unsigned char* cls; int cpitch, ccx0, ccy0, cw, ch;
    // owner keys per lattice point of the owned boxes (+1 ring so a walk can
    // test the point it exits to): vertex lattice and quad lattice
    Key* ownV; int ovpitch, ovx0, ovy0, ovw, ovh;
    Key* ownQ; int oqpitch, oqx0, oqy0, oqw, oqh;
    // resume keys (atomicMin) on the outer ring of the vertex box (descending)
    // and of the quad box (ascending)
    Key* resumeV; int RoutV;
    Key* resumeQ; int RoutQ;
};

struct LineParams {
    int X, Y;
    int hopCells, iterCap, maxHops;
    int walkIter;            // total Euler iterations per walk (CPU: 500 per streamline); <=0 = unlimited
    int bfsDepth;            // seed BFS pops per saddle (CPU: 10)
    int phase;               // 1 = ownership, 2 = marking
    float errThresh, gradThresh;
    int dx[8], dy[8];
};

enum WalkEnd { kEndTerminal = 0, kEndBoundary, kEndExit, kEndCap, kEndMerged, kEndCount };

__device__ __forceinline__ int ringIndexDev(int x, int y, int w, int h) {
    if (y == 0) return x;
    if (y == h - 1) return w + x;
    if (x == 0) return 2 * w + (y - 1);
    if (x == w - 1) return 2 * w + (h - 2) + (y - 1);
    return -1;
}

__device__ __forceinline__ float fval(const LineTile& t, int gx, int gy) {
    return t.vals[(gy - t.sy0) * t.vpitch + (gx - t.sx0)];
}
__device__ __forceinline__ bool gradInBox(const LineTile& t, int cx, int cy) {
    return cx >= t.gcx0 && cy >= t.gcy0 && cx < t.gcx0 + t.gw && cy < t.gcy0 + t.gh;
}
__device__ __forceinline__ unsigned gcode(const LineTile& t, int cx, int cy) {
    return (t.grad[(cy - t.gcy0) * t.gpitch + (cx - t.gcx0)] >> 2) & 7u;
}
// byte OR through the containing 32-bit word
__device__ __forceinline__ void markCell(const LineTile& t, int cx, int cy, unsigned cls) {
    if (cx < t.ccx0 || cy < t.ccy0 || cx >= t.ccx0 + t.cw || cy >= t.ccy0 + t.ch) return;
    unsigned char* p = t.cls + (size_t)(cy - t.ccy0) * t.cpitch + (cx - t.ccx0);
    const size_t addr = (size_t)p;
    unsigned* w = (unsigned*)(addr & ~(size_t)3);
    const unsigned shift = (unsigned)(addr & 3) * 8u;
    atomicOr(w, cls << shift);
}

// Ownership at a lattice point. Phase 1: claim (atomicMin). Phase 2: true
// if a more senior walk owns the point (the caller marks it and stops).
template <bool ASC>
__device__ __forceinline__ bool ownerStop(const LineTile& t, const LineParams& P, int lx, int ly, Key key) {
    Key* arr; int pitch, x0, y0, w, h;
    if (ASC) { arr = t.ownQ; pitch = t.oqpitch; x0 = t.oqx0; y0 = t.oqy0; w = t.oqw; h = t.oqh; }
    else     { arr = t.ownV; pitch = t.ovpitch; x0 = t.ovx0; y0 = t.ovy0; w = t.ovw; h = t.ovh; }
    if (lx < x0 || ly < y0 || lx >= x0 + w || ly >= y0 + h) return false;
    Key* p = arr + (size_t)(ly - y0) * pitch + (lx - x0);
    if (P.phase == 1) { atomicMin(p, key); return false; }
    return *p < key;
}

// mark a lattice point (vertex or quad) and test ownership; returns false
// when the walk must stop here (merged into a senior line)
template <bool ASC>
__device__ __forceinline__ bool visitPoint(const LineTile& t, const LineParams& P, int lx, int ly, Key key) {
    if (P.phase == 2) {
        if (ASC) markCell(t, 2 * lx + 1, 2 * ly + 1, kClassAsc); else markCell(t, 2 * lx, 2 * ly, kClassDsc);
    }
    return !ownerStop<ASC>(t, P, lx, ly, key);
}

__device__ __forceinline__ float2 vertGrad(const LineTile& t, const LineParams& P, int gx, int gy) {
    float2 g;
    g.x = (gx > 0 && gx < P.X - 1) ? 0.5f * (fval(t, gx + 1, gy) - fval(t, gx - 1, gy)) : 0.f;
    g.y = (gy > 0 && gy < P.Y - 1) ? 0.5f * (fval(t, gx, gy + 1) - fval(t, gx, gy - 1)) : 0.f;
    return g;
}

// Digitize the segment a -> b as a lattice chain: mark (phase 2) and claim /
// test ownership at every lattice point reached. Returns false to stop.
// ASC: quad-edge-quad chain (unshifted); DSC: vertex-edge-vertex chain
// (coordinates shifted by +0.5 so floor() names the nearest vertex). A
// segment is at most ~1.5 cells long, so each axis is crossed at most once.
template <bool ASC>
__device__ __forceinline__ bool digitize(const LineTile& t, const LineParams& P, float2 a, float2 b, Key key) {
    if (!ASC) { a.x += 0.5f; a.y += 0.5f; b.x += 0.5f; b.y += 0.5f; }
    int cx = (int)floorf(a.x), cy = (int)floorf(a.y);   // current lattice point
    const int ex = (int)floorf(b.x), ey = (int)floorf(b.y);
    float tx = 2.f, ty = 2.f;
    int sx = 0, sy = 0;
    if (ex != cx) { sx = ex > cx ? 1 : -1; const float plane = (float)(sx > 0 ? cx + 1 : cx); tx = (plane - a.x) / (b.x - a.x); }
    if (ey != cy) { sy = ey > cy ? 1 : -1; const float plane = (float)(sy > 0 ? cy + 1 : cy); ty = (plane - a.y) / (b.y - a.y); }
    const unsigned cls = ASC ? kClassAsc : kClassDsc;
    for (int k = 0; k < 2; k++) {
        if (sx != 0 && (sy == 0 || tx <= ty)) {
            if (P.phase == 2) markCell(t, ASC ? 2 * cx + 1 + sx : 2 * cx + sx, ASC ? 2 * cy + 1 : 2 * cy, cls);
            cx += sx;
            if (!visitPoint<ASC>(t, P, cx, cy, key)) return false;
            sx = 0; tx = 2.f;
        } else if (sy != 0) {
            if (P.phase == 2) markCell(t, ASC ? 2 * cx + 1 : 2 * cx, ASC ? 2 * cy + 1 + sy : 2 * cy + sy, cls);
            cy += sy;
            if (!visitPoint<ASC>(t, P, cx, cy, key)) return false;
            sy = 0; ty = 2.f;
        }
    }
    return true;
}

// quad total order: (max corner value, max corner index, quad index)
__device__ __forceinline__ void quadKey(const LineTile& t, const LineParams& P, int qx, int qy,
                                        float& f, long long& vi) {
    f = fval(t, qx, qy); vi = (long long)qx + (long long)qy * P.X;
    const int cs[3][2] = {{qx + 1, qy}, {qx, qy + 1}, {qx + 1, qy + 1}};
    for (int k = 0; k < 3; k++) {
        const float fc = fval(t, cs[k][0], cs[k][1]);
        const long long ic = (long long)cs[k][0] + (long long)cs[k][1] * P.X;
        if (!numeric2d::vlessI(fc, ic, f, vi)) { f = fc; vi = ic; }
    }
}

enum HopResult { kHopNext = 0, kHopMerged = 1 };

// One hop of a line walk from lattice point (lx, ly). Returns kHopNext with
// the successor lattice point in nx, ny, or kHopMerged when the digitized
// path touched a senior line (phase 2). Charges the walk budget.
template <bool ASC>
__device__ int lineHop(const LineTile& t, const LineParams& P, int lx, int ly, int& nx, int& ny,
                       int& budget, Key key) {
    const float sgn = ASC ? 1.f : -1.f;
    float2 p = ASC ? make_float2(lx + 0.5f, ly + 0.5f) : make_float2((float)lx, (float)ly);
    int bx = lx, by = ly;
    if (bx > P.X - 2) bx = P.X - 2;
    if (by > P.Y - 2) by = P.Y - 2;
    const int ccx = ASC ? 2 * lx + 1 : 2 * lx, ccy = ASC ? 2 * ly + 1 : 2 * ly;
    const unsigned code = gcode(t, ccx, ccy);
    const int fbx = lx + P.dx[code], fby = ly + P.dy[code];
    const unsigned cls = ASC ? kClassAsc : kClassDsc;
    int iters = P.iterCap < budget ? P.iterCap : budget;
    const int itersStart = iters;
    // discrete V-path step: mark the shared edge, then the neighbour
    auto fallback = [&]() {
        budget -= itersStart - (iters < 0 ? 0 : iters);
        if (P.phase == 2) markCell(t, ccx + P.dx[code], ccy + P.dy[code], cls);
        nx = fbx; ny = fby;
        return (int)kHopNext;
    };
    float f0; long long i0;
    if (ASC) quadKey(t, P, lx, ly, f0, i0);
    else { f0 = fval(t, lx, ly); i0 = (long long)lx + (long long)ly * P.X; }

    int cells = 0;
    bool snapNow = false;
    while (!snapNow) {
        float2 g[4];
        g[0] = vertGrad(t, P, bx, by);
        g[1] = vertGrad(t, P, bx + 1, by);
        g[2] = vertGrad(t, P, bx, by + 1);
        g[3] = vertGrad(t, P, bx + 1, by + 1);
        float2 old = p;
        float2 grad = numeric2d::bilerp(p, bx, by, g);
        float mag = sqrtf(grad.x * grad.x + grad.y * grad.y);
        if (mag < P.gradThresh) return fallback();
        float step = 0.5f / mag;
        bool outOfVoxel = false;
        while (iters-- > 0) {
            if (old.x != p.x || old.y != p.y) {
                if ((int)p.x != (int)old.x || (int)p.y != (int)old.y) { outOfVoxel = true; break; }
                grad = numeric2d::bilerp(p, bx, by, g);
                mag = sqrtf(grad.x * grad.x + grad.y * grad.y);
                if (mag * step > 0.5f) step = 0.5f / mag;
                old = p;
            }
            if (mag < P.gradThresh) return fallback();
            const float2 euler = make_float2(grad.x * (sgn * step), grad.y * (sgn * step));
            const float2 half = make_float2(grad.x * (sgn * 0.5f * step), grad.y * (sgn * 0.5f * step));
            const float2 hp = make_float2(p.x + half.x, p.y + half.y);
            const float2 hg = numeric2d::bilerp(hp, bx, by, g);
            const float2 two = make_float2(hg.x * (sgn * 0.5f * step) + half.x,
                                           hg.y * (sgn * 0.5f * step) + half.y);
            const float ex = two.x - euler.x, ey = two.y - euler.y;
            if (ex * ex + ey * ey < P.errThresh) {
                const float2 prev = p;
                p.x += two.x; p.y += two.y;
                if (p.x < 0.f) p.x = 0.f;
                if (p.y < 0.f) p.y = 0.f;
                if (p.x > (float)(P.X - 1)) p.x = (float)(P.X - 1);
                if (p.y > (float)(P.Y - 1)) p.y = (float)(P.Y - 1);
                if (!digitize<ASC>(t, P, prev, p, key)) {
                    budget -= itersStart - (iters < 0 ? 0 : iters);
                    return (int)kHopMerged;
                }
                const float lim = 0.5f / mag;
                step = step * 1.8f < lim ? step * 1.8f : lim;
                if (p.x < (float)t.ox || p.y < (float)t.oy ||
                    p.x >= (float)(t.ox + t.tw) || p.y >= (float)(t.oy + t.th)) { snapNow = true; break; }
            } else {
                step *= 0.5f;
            }
        }
        if (snapNow) break;
        if (!outOfVoxel) break;   // iteration cap: snap from here
        bx = (int)p.x; by = (int)p.y;
        if (bx > P.X - 2) bx = P.X - 2;
        if (by > P.Y - 2) by = P.Y - 2;
        cells++;
        if (bx <= 0 || by <= 0 || bx >= P.X - 2 || by >= P.Y - 2) break;   // boundary band: snap
        if (cells >= P.hopCells) break;
    }
    // snap
    int sxq, syq;
    if (ASC) {
        sxq = bx; syq = by;   // the quad containing p = the digitizer's current quad
        float fq; long long iq;
        quadKey(t, P, sxq, syq, fq, iq);
        const bool greater = numeric2d::vlessI(f0, i0, fq, iq) ||
                             (fq == f0 && iq == i0 && ((long long)sxq + (long long)syq * P.X) > ((long long)lx + (long long)ly * P.X));
        if (!greater) return fallback();
    } else {
        sxq = (int)floorf(p.x + 0.5f); syq = (int)floorf(p.y + 0.5f);   // nearest vertex = digitizer's current vertex
        if (sxq > P.X - 1) sxq = P.X - 1;
        if (syq > P.Y - 1) syq = P.Y - 1;
        const float fs = fval(t, sxq, syq);
        const long long is = (long long)sxq + (long long)syq * P.X;
        if (!numeric2d::vlessI(fs, is, f0, i0)) return fallback();
    }
    budget -= itersStart - (iters < 0 ? 0 : iters);
    nx = sxq; ny = syq;
    return (int)kHopNext;
}

// Full walk from a seed lattice point. Returns the WalkEnd event.
template <bool ASC>
__device__ int walk(const LineTile& t, const LineParams& P, int lx, int ly, Key key) {
    const int bw = ASC ? t.qw : t.tw, bh = ASC ? t.qh : t.th;
    int budget = P.walkIter > 0 ? P.walkIter : 0x3FFFFFFF;
    for (int h = 0; h < P.maxHops; h++) {
        const int rx = lx - t.ox, ry = ly - t.oy;
        if (rx < 0 || ry < 0 || rx >= bw || ry >= bh) {
            // hand over to the neighbour tile: the most senior key wins the slot
            const int r = ringIndexDev(rx + 1, ry + 1, bw + 2, bh + 2);
            if (r >= 0) atomicMin(ASC ? t.resumeQ + r : t.resumeV + r, key);
            return kEndExit;
        }
        if (!visitPoint<ASC>(t, P, lx, ly, key)) return kEndMerged;
        const int ccx = ASC ? 2 * lx + 1 : 2 * lx, ccy = ASC ? 2 * ly + 1 : 2 * ly;
        if (gcode(t, ccx, ccy) == 7u) return kEndTerminal;
        const int NX = ASC ? P.X - 1 : P.X, NY = ASC ? P.Y - 1 : P.Y;
        if (lx < 2 || ly < 2 || lx > NX - 3 || ly > NY - 3) return kEndBoundary;
        if (budget <= 0) return kEndCap;
        int nx, ny;
        if (lineHop<ASC>(t, P, lx, ly, nx, ny, budget, key) == kHopMerged) return kEndMerged;
        lx = nx; ly = ny;
    }
    return kEndCap;
}

} // namespace lines2d
} // namespace gpu
} // namespace GInt

#endif
