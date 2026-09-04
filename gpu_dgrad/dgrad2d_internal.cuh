#ifndef DGRAD2D_INTERNAL_CUH
#define DGRAD2D_INTERNAL_CUH

// Internal (library-private) contract between the 2D kernel TU and its
// drivers. CUDA types allowed here; never included by a host-only TU.

#include "dgrad_gpu_api.h"

#include <cuda_runtime.h>

namespace GInt {
namespace gpu {

// One unit of work for the fused kernel: a box of vertices to compute, the
// staged function values covering their 1-ring, and the owned cell box that
// receives the gradient bytes. The in-core batched path is the special case
// "one view = one whole slice".
struct DeviceTileView2D {
    const float* vals;   // staged values, row-major, x-fastest
    int pitch;           // row pitch of vals (elements)
    int sx0, sy0;        // global vertex coords of vals[0]
    int svx, svy;        // staged box size
    int cx0, cy0;        // first computed vertex (global)
    int cw, ch;          // computed box size
    int ox, oy;          // first owned vertex (global)
    int tw, th;          // owned vertices per axis
    int X, Y;            // image size (for validity + boundary classes)
    unsigned char* grad; // owned cell box, row-major
    int gpitch;          // row pitch of grad (bytes)
    int gcx0, gcy0;      // global cell coords of grad[0]
};

// Uploads the topology tables to the kernel TU's constant memory (once per
// driver call; cheap). Synchronous.
bool UploadTables2D(const Dgrad2DTables& tables);

// Launches the fused full-domain kernel over nviews views (device array),
// gridding the max computed box; each block picks its view by blockIdx.z.
// Returns false on a launch error. Asynchronous on 'stream'.
bool LaunchFused2DViews(const DeviceTileView2D* d_views, int nviews,
                        int maxComputedW, int maxComputedH, int block_shape,
                        cudaStream_t stream);

} // namespace gpu
} // namespace GInt

#endif
