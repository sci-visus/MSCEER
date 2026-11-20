#ifndef DGRAD3D_INTERNAL_CUH
#define DGRAD3D_INTERNAL_CUH

// Internal contract between the 3D lower-star kernel TU and its drivers.
// CUDA types allowed; never included by a host-only TU.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>

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
// gridding the largest computed box. Asynchronous on 'stream'.
bool LaunchDgrad3DViews(const DeviceTileView3D* d_views, int nviews, int maxCw, int maxCh,
                        int maxCd, int block_shape, cudaStream_t stream);

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
