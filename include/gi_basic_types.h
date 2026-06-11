#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

// Type definitions for mesh indexing and computation
// These are macros for backwards compatibility, but could be converted to using/typedef
#define INDEX_TYPE long long		// type to index elements of a mesh - need to support big meshes! - keep this signed to make arithmetic consistent
#define FLOATTYPE float				// type used for computing numerical integration
#define INT_TYPE int				// regular old ints
#define DIM_TYPE unsigned char		// used to query the dimension of a cell - we usually have values between 0-3
#define BYTE_TYPE unsigned char		// just a regular old byte
#define BOUNDARY_TYPE unsigned char	// used to query the boundary classification of a cell - usually a small integer
#define ASSIGNED_TYPE unsigned char	// used to query if a cell has been assigned - usually 0, 1 values

// Parallel execution backend selection
// Only one should be active at a time. If none are defined, defaults to single-threaded.
// GINT_USE_OPENMP - Use OpenMP for CPU parallelism
// GINT_USE_SYCL   - Use SYCL for GPU offloading (mutually exclusive with OpenMP macros)
// GINT_USE_STDPAR - Use C++17 parallel algorithms (mutually exclusive with OpenMP macros)

// OpenMP compatibility layer
// Only provide OpenMP stubs if:
// 1. OpenMP is available (_OPENMP defined by compiler), OR
// 2. We're explicitly using OpenMP backend (GINT_USE_OPENMP defined)
// Don't provide stubs if using SYCL or stdpar to avoid conflicts
#if defined(_OPENMP) && !defined(GINT_USE_SYCL) && !defined(GINT_USE_STDPAR)
#include <omp.h>
#define GINT_OPENMP_ENABLED 1
#elif !defined(GINT_USE_SYCL) && !defined(GINT_USE_STDPAR)
// Provide single-threaded stubs only when not using alternative parallel backends
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define GINT_OPENMP_ENABLED 0
#endif

// Marker for device-compatible code (used in SYCL/CUDA contexts)
#ifdef GINT_USE_SYCL
#define GINT_DEVICE_FUNC // Could add SYCL-specific attributes here
#else
#define GINT_DEVICE_FUNC
#endif


enum DestType : BYTE_TYPE {
	BACKGROUND,
	UNASSIGNED,
	ASSIGNED,
	CERTAIN_TERMINAL,
	CERTAIN_NONTERMINAL
};
  
#endif
