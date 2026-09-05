#ifndef DGRAD_NUMERIC2D_CUH
#define DGRAD_NUMERIC2D_CUH

// Stage A2: device numeric streamline hop for the 2D vertex lattice
// (descending flow toward minima). Shared by the region successor (A2) and
// the line walks (A4).
//
// Ported line for line from the CPU:
//   gradient at a vertex = central differences GradientFromImage(p, 1) with
//     the stencil reduced to nothing (component 0) ON the domain edge
//     (Gather1DNeighborhood clamps the kernel to DistToBoundary);
//   gradient at a point  = bilinear interpolation of the 4 corner gradients
//     of the containing voxel (BiLinInterpGrad with a FIXED voxel base);
//   integrator           = AdaptiveEulerAdvector2D::AdvectThroughVoxelNoCheck:
//     step <= 0.5 cell, one Euler step vs two half steps (the second half
//     step re-interpolates with the SAME voxel base - the CPU's "incorrect
//     method", reproduced on purpose so agreement is measurable), accept if
//     the squared displacement difference < errorThreshold, grow x1.8 /
//     shrink x0.5, StepMultiplier = -1 for descent.
//
// What is NOT the CPU's: termination. A hop ends when
//   (a) it enters a voxel whose 4 corners are all certain with one living
//       id -> TERMINAL(living id);
//   (b) it leaves the tile's owned vertex box -> snap from the voxel it was
//       in (its corners are at most one step outside the box);
//   (c) the current voxel touches the domain boundary -> snap (the CPU
//       switches to its clamped near-boundary advector there);
//   (d) |g| < gradThreshold -> discrete fallback (plateaus are the base
//       gradient's business: index tie-breaking);
//   (e) hopCells voxels traversed -> snap;
//   (f) iterCap Euler iterations -> snap.
// Snap = the lowest corner of the voxel under the (value, index) total
// order, accepted only if strictly below the start vertex - otherwise the
// discrete V-path step v + 2*pairdir (strictly lower by construction). That
// guard is what makes the successor field acyclic for ANY vector field.

#include <cuda_runtime.h>
#include <cstdint>

namespace GInt {
namespace gpu {
namespace numeric2d {

struct TileIn {
    const float* vals;         // staged values (owned box + 2 halo, clipped)
    int vpitch, sx0, sy0, svx, svy;
    const int32_t* base;       // compact base labels (owned box + 1 halo, clipped)
    int bpitch, bx0, by0, bw, bh;
    const unsigned char* grad; // owned cell box, origin (2ox, 2oy)
    int gpitch;
    int ox, oy, tw, th;        // owned vertex box
};

struct Params {
    int X, Y;
    int hopCells, iterCap;
    float errThresh, gradThresh;
    const int32_t* remap;      // [nBase] -> living id
    const float* sStar;        // [nLiving] lowest boundary saddle value
    int nBase, nLiving;
    int dx[8], dy[8];          // pair code -> lattice step (code 7 = critical)
};

enum HopEvent {
    kEvCertain = 0,       // decided by f(x) < s*(living)
    kEvMinimum,           // critical vertex: TERMINAL(remap[base])
    kEvCornerCertain,     // entered a voxel with 4 certain corners
    kEvBoxExit,           // left the owned box -> snap
    kEvBoundary,          // voxel touches the domain boundary -> snap
    kEvLowGrad,           // |g| < threshold -> fallback
    kEvCellCap,           // hopCells reached -> snap
    kEvIterCap,           // iterCap reached -> snap
    kEvGuardFallback,     // snap corner not strictly lower -> V-path step
    kEvStartBoundary,     // vertex within 2 of the domain edge -> fallback
    kEvCount
};

__device__ __forceinline__ float fval(const TileIn& t, int gx, int gy) {
    return t.vals[(gy - t.sy0) * t.vpitch + (gx - t.sx0)];
}
__device__ __forceinline__ int32_t baseOf(const TileIn& t, int gx, int gy) {
    return t.base[(gy - t.by0) * t.bpitch + (gx - t.bx0)];
}
__device__ __forceinline__ int livingOf(const TileIn& t, const Params& P, int gx, int gy) {
    const int32_t b = baseOf(t, gx, gy);
    return (b < 0 || b >= P.nBase) ? -1 : P.remap[b];
}
__device__ __forceinline__ bool isCertain(const TileIn& t, const Params& P, int gx, int gy,
                                          int& living) {
    living = livingOf(t, P, gx, gy);
    return living >= 0 && fval(t, gx, gy) < P.sStar[living];
}
__device__ __forceinline__ float2 vertGrad(const TileIn& t, const Params& P, int gx, int gy) {
    float2 g;
    g.x = (gx > 0 && gx < P.X - 1) ? 0.5f * (fval(t, gx + 1, gy) - fval(t, gx - 1, gy)) : 0.f;
    g.y = (gy > 0 && gy < P.Y - 1) ? 0.5f * (fval(t, gx, gy + 1) - fval(t, gx, gy - 1)) : 0.f;
    return g;
}
__device__ __forceinline__ float2 lerp2(float2 a, float2 b, float f) {
    return make_float2(a.x + (b.x - a.x) * f, a.y + (b.y - a.y) * f);
}
__device__ __forceinline__ float2 bilerp(float2 p, int bx, int by, const float2 (&g)[4]) {
    const float fx = p.x - (float)bx, fy = p.y - (float)by;
    const float2 x0 = lerp2(g[0], g[1], fx);
    const float2 x1 = lerp2(g[2], g[3], fx);
    return lerp2(x0, x1, fy);
}
// (value, index) total order: a < b
__device__ __forceinline__ bool vlessI(float fa, long long ia, float fb, long long ib) {
    if (fa < fb) return true;
    if (fb < fa) return false;
    return ia < ib;
}

// Successor decision for the start vertex x = (gx, gy). Returns the global
// lattice index of the successor vertex, or -1 with *living set when the
// vertex is a terminal (certain / minimum / corner-certain). ev receives the
// terminating event. Caller converts the successor to a word.
__device__ inline long long hop(const TileIn& t, const Params& P, int gx, int gy,
                                int& living, int& ev) {
    const long long X = P.X;
    const long long xi = gx + gy * X;
    const float fx0 = fval(t, gx, gy);

    // decide first
    if (isCertain(t, P, gx, gy, living)) { ev = kEvCertain; return -1; }
    const unsigned code = (t.grad[(2 * (gy - t.oy)) * t.gpitch + 2 * (gx - t.ox)] >> 2) & 7u;
    if (code == 7u) { ev = kEvMinimum; return -1; }   // living already = remap[base]

    // discrete fallback target (always strictly lower under the total order)
    const long long fb = (long long)(gx + P.dx[code]) + (long long)(gy + P.dy[code]) * X;

    if (gx < 2 || gy < 2 || gx > P.X - 3 || gy > P.Y - 3) { ev = kEvStartBoundary; return fb; }

    float2 p = make_float2((float)gx, (float)gy);
    int bx = gx, by = gy;   // current voxel base
    if (bx > P.X - 2) bx = P.X - 2;
    if (by > P.Y - 2) by = P.Y - 2;
    int cells = 0;
    int iters = P.iterCap;
    int snapEv = -1;

    while (snapEv < 0) {
        // ---- AdvectThroughVoxelNoCheck on voxel (bx, by) ----
        float2 g[4];
        g[0] = vertGrad(t, P, bx, by);
        g[1] = vertGrad(t, P, bx + 1, by);
        g[2] = vertGrad(t, P, bx, by + 1);
        g[3] = vertGrad(t, P, bx + 1, by + 1);
        float2 old = p;
        float2 grad = bilerp(p, bx, by, g);
        float mag = sqrtf(grad.x * grad.x + grad.y * grad.y);
        if (mag < P.gradThresh) { ev = kEvLowGrad; return fb; }
        float step = 0.5f / mag;
        bool outOfVoxel = false;
        while (iters-- > 0) {
            if (old.x != p.x || old.y != p.y) {
                if ((int)p.x != (int)old.x || (int)p.y != (int)old.y) { outOfVoxel = true; break; }
                grad = bilerp(p, bx, by, g);
                mag = sqrtf(grad.x * grad.x + grad.y * grad.y);
                if (mag * step > 0.5f) step = 0.5f / mag;
                old = p;
            }
            if (mag < P.gradThresh) { ev = kEvLowGrad; return fb; }
            const float2 euler = make_float2(grad.x * (-step), grad.y * (-step));
            const float2 half = make_float2(grad.x * (-0.5f * step), grad.y * (-0.5f * step));
            const float2 hp = make_float2(p.x + half.x, p.y + half.y);
            const float2 hg = bilerp(hp, bx, by, g);
            const float2 two = make_float2(hg.x * (-0.5f * step) + half.x,
                                           hg.y * (-0.5f * step) + half.y);
            const float ex = two.x - euler.x, ey = two.y - euler.y;
            if (ex * ex + ey * ey < P.errThresh) {
                p.x += two.x;
                p.y += two.y;
                const float lim = 0.5f / mag;
                step = step * 1.8f < lim ? step * 1.8f : lim;
                // left the owned box? snap from the voxel we are in
                if (p.x < (float)t.ox || p.y < (float)t.oy ||
                    p.x >= (float)(t.ox + t.tw) || p.y >= (float)(t.oy + t.th)) {
                    snapEv = kEvBoxExit;
                    break;
                }
            } else {
                step *= 0.5f;
            }
        }
        if (snapEv >= 0) break;
        if (!outOfVoxel) { snapEv = kEvIterCap; break; }
        // ---- entered a new voxel ----
        bx = (int)p.x;
        by = (int)p.y;
        if (bx > P.X - 2) bx = P.X - 2;
        if (by > P.Y - 2) by = P.Y - 2;
        cells++;
        if (bx <= 0 || by <= 0 || bx >= P.X - 2 || by >= P.Y - 2) { snapEv = kEvBoundary; break; }
        int l0, l1, l2, l3;
        if (isCertain(t, P, bx, by, l0) && isCertain(t, P, bx + 1, by, l1) &&
            isCertain(t, P, bx, by + 1, l2) && isCertain(t, P, bx + 1, by + 1, l3) &&
            l0 == l1 && l0 == l2 && l0 == l3) {
            living = l0;
            ev = kEvCornerCertain;
            return -1;
        }
        if (cells >= P.hopCells) { snapEv = kEvCellCap; break; }
    }

    // ---- snap: lowest corner of voxel (bx, by), guarded ----
    long long best = -1;
    float bestf = 0.f;
    for (int j = 0; j < 2; j++)
        for (int i = 0; i < 2; i++) {
            const int cx = bx + i, cy = by + j;
            const float fc = fval(t, cx, cy);
            const long long ci = (long long)cx + (long long)cy * X;
            if (best < 0 || vlessI(fc, ci, bestf, best)) { best = ci; bestf = fc; }
        }
    if (vlessI(bestf, best, fx0, xi)) { ev = snapEv; return best; }
    ev = kEvGuardFallback;
    return fb;
}

} // namespace numeric2d
} // namespace gpu
} // namespace GInt

#endif
