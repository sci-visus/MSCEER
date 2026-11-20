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
    dgrad2d_tiled.cu
    successor2d.cu
    lines2d.cu
    gpu_transfer.cu
    dmt3d_gradient.cu
    dgrad3d_kernel.cu
    dgrad3d_tiled.cu
    dmt3d_tiled.cu
    dgrad_gpu_api.h
    dgrad2d_tables.h
    tile2d.h
    dgrad2d_internal.cuh
    gpu_transfer.h
    frontier2d.h
    numeric2d.cuh
    lines2d.cuh
    tile3d.h
    dgrad3d_tables.h
    dgrad3d_internal.cuh
)
target_include_directories(gpu_dgrad PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Host code of the .cu files uses OpenMP for the in-core adapters' copies
# (the project requires OpenMP outside Apple, see cmake/Common.cmake).
if(OpenMP_CXX_FOUND)
    target_compile_options(gpu_dgrad PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=${OpenMP_CXX_FLAGS}>)
    target_link_libraries(gpu_dgrad PUBLIC OpenMP::OpenMP_CXX)
endif()

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

# Stage A0 gate: tiled == in-core == CPU, byte-exact, plus the constant-
# footprint run on a procedural image larger than device memory.
add_executable(dgrad2d_tiled_validate dgrad2d_tiled_validate.cxx)
target_link_libraries(dgrad2d_tiled_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(dgrad2d_tiled_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A1 gate: tiled successor resolver == in-core flow terminals.
add_executable(succ2d_validate succ2d_validate.cxx)
target_link_libraries(succ2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(succ2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A2 gate: tiled numeric region labeling (acyclicity, decrease,
# tiled == in-core, agreement vs host-integrated basins).
add_executable(numeric2d_validate numeric2d_validate.cxx)
target_link_libraries(numeric2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(numeric2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A4 gate: tiled numeric separatrix lines.
add_executable(lines2d_validate lines2d_validate.cxx)
target_link_libraries(lines2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(lines2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A5 gate: constrained fused kernel == CPU constrained MyRobinsNoalloc.
add_executable(dgrad2d_constrained_validate dgrad2d_constrained_validate.cxx)
target_link_libraries(dgrad2d_constrained_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(dgrad2d_constrained_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# 3D discrete Morse theory rules (PLS, PLS-BR, DMF): GPU drivers ==
# CPU raw drivers of include/gi_dmt_*.h, byte-exact, plus timings.
add_executable(dmt3d_validate dmt3d_validate.cxx)
target_link_libraries(dmt3d_validate PRIVATE GInt gpu_dgrad)
target_compile_features(dmt3d_validate PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(dmt3d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(dmt3d_validate PRIVATE /bigobj)
endif()

# 3D lower-star (Robins) gradient: GPU (in-core and tiled, unconstrained and
# constrained) == CPU SlidingWindowRobinsNoalloc::ComputePairing_sliding /
# MyRobinsNoalloc, byte-exact; tiled == in-core; fixed-footprint large run.
add_executable(dgrad3d_validate dgrad3d_validate.cxx)
target_link_libraries(dgrad3d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(dgrad3d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(dgrad3d_validate PRIVATE /bigobj)
endif()

# Tiled DMT drivers (multi-pass DMF, PLS, PLS-BR): tiled == in-core == CPU,
# halo sweep, fixed-footprint large run.
add_executable(dmt3d_tiled_validate dmt3d_tiled_validate.cxx)
target_link_libraries(dmt3d_tiled_validate PRIVATE GInt gpu_dgrad)
target_compile_features(dmt3d_tiled_validate PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(dmt3d_tiled_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(dmt3d_tiled_validate PRIVATE /bigobj)
endif()

# GPU-only benchmark of the tiled 3D gradient paths (Robins, PLS, PLS-BR, DMF).
add_executable(gpu3d_bench gpu3d_bench.cxx)
target_link_libraries(gpu3d_bench PRIVATE GInt gpu_dgrad)
target_compile_features(gpu3d_bench PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(gpu3d_bench PRIVATE OpenMP::OpenMP_CXX)
endif()

# The Stage 4 parity test (msc2d_gpu_parity.cxx) is added from the ROOT
# CMakeLists after all subdirectories, because it needs msc_2d_lib, which is
# configured after this directory.
