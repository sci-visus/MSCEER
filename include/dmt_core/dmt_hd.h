#ifndef DMT_HD_H
#define DMT_HD_H

// Device-safe core of the discrete Morse theory kernels.
//
// Every header under include/dmt_core/ obeys the same rules so that the exact
// same source compiles as host code and as CUDA device code:
//   - includes only <cstdint> and other dmt_core headers (never gi_*.h);
//   - no heap allocation, STL containers, exceptions, recursion or virtual calls;
//   - fixed-size state (at most 27 slots for a 3x3x3 star) and integer masks;
//   - every function is marked DMT_HD.
// Host-side drivers, memo structures and analysis live in include/gi_dmt_*.h.

#include <cstdint>

#if defined(__CUDACC__)
#define DMT_HD __host__ __device__ __forceinline__
#define DMT_UNROLL _Pragma("unroll")
#else
#define DMT_HD inline
#define DMT_UNROLL
#endif

namespace dmt {

	// One byte per cell, identical to GInt::GradBitfield
	// (gi_discrete_gradient_labeling.h) and to the gpu_dgrad output contract:
	//   bit 0     assigned
	//   bit 1     flag
	//   bits 2-4  pair direction 0..5 (+x,-x,+y,-y,+z,-z), 7 = critical
	//   bits 5-7  ldir
	constexpr uint8_t kPairCritical = 7;

	DMT_HD uint8_t MakeGradByte(uint8_t dir) { return (uint8_t)(1u | ((dir & 7u) << 2)); }
	DMT_HD uint8_t GradPairCode(uint8_t b) { return (uint8_t)((b >> 2) & 7u); }
	DMT_HD bool GradAssigned(uint8_t b) { return (b & 1u) != 0; }
	DMT_HD bool GradCritical(uint8_t b) { return GradPairCode(b) == kPairCritical; }

} // namespace dmt

#endif
