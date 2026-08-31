# GPU discrete gradient path.
#
# Structure (established in Stage 0, filled in by later stages):
#   dgrad_gpu_api.h    - host-facing interface; plain C++ header, includable by any
#                        existing tool without a CUDA compiler
#   dgrad_gpu_api.cu   - host-side launcher implementations
#   (later)            - standalone device-side kernel headers; these must NEVER
#                        include the gi_* template stack (see the cuda-gradient
#                        master plan, Appendix A.5) - the CPU and GPU paths share
#                        only generated topology tables and a byte-exact output
#                        contract with include/gi_discrete_gradient_labeling.h

if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 12.8)
    message(WARNING
        "CUDA ${CMAKE_CUDA_COMPILER_VERSION} cannot generate native sm_120 (Blackwell) code. "
        "With the Visual Studio generator select the toolkit via -T cuda=13.0; "
        "with Ninja/Makefiles pass -DCMAKE_CUDA_COMPILER=<path to CUDA >= 12.8 nvcc>.")
endif()

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Host-facing GPU library. Deliberately NOT linked against GInt: device code and
# its launchers stay independent of the gi_* headers.
add_library(gpu_dgrad STATIC
    dgrad_gpu_api.cu
    dgrad2d_kernel.cu
    dgrad2d_labeling.cu
    dgrad2d_label_ctx.cu
    dgrad_gpu_api.h
    dgrad2d_tables.h
)
target_include_directories(gpu_dgrad PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Carry the CUDA runtime on the INTERFACE: CMake only injects it automatically
# for targets configured where the CUDA language is enabled, so a pure-C++
# consumer in an embedding project (e.g. MSSeg linking gpu_dgrad transitively
# through msc_2d_lib) would otherwise hit LNK2019 on cudaMalloc & co. Static
# cudart is kept deliberately so executables run without the CUDA bin on PATH.
find_package(CUDAToolkit REQUIRED)
target_link_libraries(gpu_dgrad PUBLIC CUDA::cudart_static)

# Smoke test / benchmark driver: device query, trivial kernel round-trip,
# pinned H2D/D2H and device-to-device bandwidth measurement.
add_executable(gpu_dgrad_smoke smoke_test.cu)
target_link_libraries(gpu_dgrad_smoke PRIVATE gpu_dgrad)

# Stage 1 validation harness: CPU 2D reference vs GPU kernel, bit-exact.
# Host-only C++ (compiled by MSVC, not nvcc) - it may include the gi_* stack.
add_executable(dgrad2d_validate dgrad2d_validate.cxx)
target_link_libraries(dgrad2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(dgrad2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# The Stage 4 parity test (msc2d_gpu_parity.cxx) is added from the ROOT
# CMakeLists after all subdirectories, because it needs msc_2d_lib, which is
# configured after this directory.
