#ifndef DGRAD3D_INTERNAL_CUH
#define DGRAD3D_INTERNAL_CUH

// Internal contract between the 3D lower-star kernel TU and its drivers.
// CUDA types allowed; never included by a host-only TU.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>

#include <functional>
#include <vector>

namespace GInt {
namespace gpu {

// One unit of work for the 3D kernel. Per-axis coordinates are global
// (int32; the drivers require every axis < 2^30); every linear index the
// kernel forms is tile-local int32 (staged box, owned cell box, class box).
struct DeviceTileView3D {
    const float* vals;           // staged values, x fastest
    int vpy, vpz;                // row / slice pitch of vals (elements)
    int sx0, sy0, sz0;           // global vertex coords of vals[0]
    int svx, svy, svz;           // staged box size
    int cx0, cy0, cz0;           // first computed vertex
    int cw, ch, cd;              // computed box size
    int ox, oy, oz;              // first owned vertex
    int tw, th, td;              // owned vertices per axis
    int X, Y, Z;                 // volume size (vertices)
    unsigned char* grad;         // owned cell box, x fastest
    int gpy, gpz;                // row / slice pitch of grad (bytes)
    int gcx0, gcy0, gcz0;        // global cell coords of grad[0]
    const unsigned char* cls;    // restriction classes (null = unconstrained)
    int kpy, kpz;                // row / slice pitch of cls
    int kcx0, kcy0, kcz0;        // global cell coords of cls[0]
    int keyClamp;                // max subset key (CPU: LMAXDIM*LMAX_LABEL-1)
    unsigned char dircode[8];    // GradBitfield pair code per raw direction 2*axis+neg
};

// Uploads the topology tables to the kernel TU's constant memory. Returns
// false if the tables break a kernel assumption (e.g. > 4 window cofacets).
bool UploadTables3D(const Dgrad3DTables& tables);

// Launches the kernel over nviews device views (blockIdx.z = view*nbz + bz),
// gridding the largest computed box. interior = true runs the interior kernel
// (views of domain-interior vertices, no classes; see SplitViews3D), false the
// boundary / label kernel (any view). Asynchronous on 'stream'.
// block_shape 0-2: 32x4x2, 16x8x2, 8x8x4; kShapeFace{X,Y,Z}: 1x32x8, 32x1x8,
// 32x8x1 (boundary kernel only), for views one vertex thick across that axis.
enum { kShapeFaceX = 16, kShapeFaceY = 17, kShapeFaceZ = 18 };
bool LaunchDgrad3DViews(const DeviceTileView3D* d_views, int nviews, int maxCw, int maxCh,
                        int maxCd, int block_shape, cudaStream_t stream, bool interior = false);

// A batch of views, regrouped for the two kernels: group 0 = the interior
// boxes (interior kernel), groups 1-3 = the domain-face boxes across x, y and
// z (boundary kernel, each with a block one vertex thick across its axis),
// group 4 = constrained views, not split: the classes are matched at every
// vertex, so all of them go to the boundary / label kernel whole.
constexpr int kViewGroups3D = 5;
struct ViewGroups3D {
    std::vector<DeviceTileView3D> views;   // group after group
    int offset[kViewGroups3D] = {0, 0, 0, 0, 0}, count[kViewGroups3D] = {0, 0, 0, 0, 0};
    int maxC[kViewGroups3D][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

// Splits each view's computed box into the domain interior (every axis in
// [1, N-2]) and the shell of domain-face vertices: x faces over the full y / z
// extents, y faces over the interior x range, z faces over the interior x and
// y ranges, so every computed vertex is in exactly one box. Ownership, staged
// values and outputs are unchanged, so the result is byte for byte the same.
inline ViewGroups3D SplitViews3D(const DeviceTileView3D* in, int n) {
    ViewGroups3D G;
    std::vector<DeviceTileView3D> grp[kViewGroups3D];
    for (int k = 0; k < n; k++) {
        const DeviceTileView3D& v = in[k];
        if (v.cls != nullptr) {   // labels: every vertex through the label kernel
            grp[4].push_back(v);
            continue;
        }
        const int c0[3] = {v.cx0, v.cy0, v.cz0}, cn[3] = {v.cw, v.ch, v.cd}, N[3] = {v.X, v.Y, v.Z};
        int ilo[3], ihi[3];
        for (int a = 0; a < 3; a++) {
            ilo[a] = c0[a] > 1 ? c0[a] : 1;
            ihi[a] = c0[a] + cn[a] < N[a] - 1 ? c0[a] + cn[a] : N[a] - 1;
            if (ihi[a] < ilo[a]) ilo[a] = ihi[a] = c0[a] + cn[a];   // no interior on this axis
            if (ilo[a] < c0[a]) ilo[a] = c0[a];
            if (ihi[a] < ilo[a]) ihi[a] = ilo[a];
        }
        auto box = [&](int g, int x0, int x1, int y0, int y1, int z0, int z1) {
            if (x1 <= x0 || y1 <= y0 || z1 <= z0) return;
            DeviceTileView3D u = v;
            u.cx0 = x0; u.cw = x1 - x0;
            u.cy0 = y0; u.ch = y1 - y0;
            u.cz0 = z0; u.cd = z1 - z0;
            grp[g].push_back(u);
        };
        const int e[3] = {c0[0] + cn[0], c0[1] + cn[1], c0[2] + cn[2]};
        box(0, ilo[0], ihi[0], ilo[1], ihi[1], ilo[2], ihi[2]);                  // interior
        box(1, c0[0], ilo[0], c0[1], e[1], c0[2], e[2]);                          // x faces
        box(1, ihi[0], e[0], c0[1], e[1], c0[2], e[2]);
        box(2, ilo[0], ihi[0], c0[1], ilo[1], c0[2], e[2]);                       // y faces
        box(2, ilo[0], ihi[0], ihi[1], e[1], c0[2], e[2]);
        box(3, ilo[0], ihi[0], ilo[1], ihi[1], c0[2], ilo[2]);                    // z faces
        box(3, ilo[0], ihi[0], ilo[1], ihi[1], ihi[2], e[2]);
    }
    for (int g = 0; g < kViewGroups3D; g++) {
        G.offset[g] = (int)G.views.size();
        G.count[g] = (int)grp[g].size();
        for (const DeviceTileView3D& u : grp[g]) {
            if (u.cw > G.maxC[g][0]) G.maxC[g][0] = u.cw;
            if (u.ch > G.maxC[g][1]) G.maxC[g][1] = u.ch;
            if (u.cd > G.maxC[g][2]) G.maxC[g][2] = u.cd;
            G.views.push_back(u);
        }
    }
    return G;
}

// Launches the groups of a split batch whose views sit at d_views (in the
// order of G.views): the interior kernel for group 0, the boundary / label
// kernel for groups 1-4 (face-shaped blocks for the face groups).
inline bool LaunchDgrad3DGroups(const DeviceTileView3D* d_views, const ViewGroups3D& G, int block_shape,
                                cudaStream_t stream) {
    const int shape[kViewGroups3D] = {block_shape, kShapeFaceX, kShapeFaceY, kShapeFaceZ, block_shape};
    for (int g = 0; g < kViewGroups3D; g++) {
        if (G.count[g] == 0) continue;
        if (!LaunchDgrad3DViews(d_views + G.offset[g], G.count[g], G.maxC[g][0], G.maxC[g][1], G.maxC[g][2],
                                shape[g], stream, g == 0))
            return false;
    }
    return true;
}

// Runs f(i) for i in [0, n) on up to `threads` threads of the library's
// persistent host worker pool (dgrad3d_tiled.cu); HostCopyThreads picks a
// thread count for a copy of `bytes` (1 below 256 KB).
bool HostParallelFor(int n, int threads, const std::function<void(int)>& f);
int HostCopyThreads(int64_t bytes);

// Packs gradient bytes two cells per byte for the download: nibble =
// pair code (bits 2-4) | assigned (bit 0) << 3; cell 2k in the low nibble of
// byte k. Lossless for every byte the lower-star kernels write (flag and ldir
// are 0). nOut (packed bytes) must be a multiple of 8. Asynchronous.
bool LaunchPackNibbles(const unsigned char* src, unsigned char* dst, size_t nOut, cudaStream_t stream);

// Inverse of the nibble code, on the host.
inline unsigned char UnpackNibble(unsigned v) {
    return (unsigned char)(((v >> 3) & 1u) | ((v & 7u) << 2));
}

// The GradBitfield pair codes TopologicalRegularGrid3D::Compress6NeighborOffsetToByte
// returns for the six unit steps of an X*Y*Z grid (index 2*axis + (step<0)).
// It matches |diff| against the offsets 1, MX, MX*MY in that order, so on a
// degenerate axis (MX == 1 or MY == 1) two axes share a code; reproduced here.
inline void DirCodes3D(int64_t X, int64_t Y, unsigned char out[8]) {
    const int64_t MX = 2 * X - 1, MY = 2 * Y - 1;
    const int64_t offs[3] = {1, MX, MX * MY};
    for (int a = 0; a < 3; a++) {
        int i = 0;
        while (i < 3 && offs[i] != offs[a]) i++;
        out[2 * a] = (unsigned char)(2 * i);
        out[2 * a + 1] = (unsigned char)(2 * i + 1);
    }
    out[6] = 7;
    out[7] = 7;
}

} // namespace gpu
} // namespace GInt

#endif
