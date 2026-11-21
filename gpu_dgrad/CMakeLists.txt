# GPU discrete gradient path.
#
# Structure:
#   api/         - host-facing interface (dgrad_gpu_api.h and the tile types); plain
#                  C++ headers, includable by any tool without a CUDA compiler
#   common/      - device query, transfer ring, frontier helpers
#   cls/         - lower-star gradient kernels (2D, 3D), their tiled drivers and the
#                  topology-table builders
#   dmt3d/       - PLS, PLS-BR and DMF drivers (in-core, tiled, streamed)
#   labeling2d/  - 2D flow terminals, numeric regions and lines
#   tests/       - validators, smoke test, benchmark
#   docs/        - README and the CPU-to-GPU guide
#
# Device-side code must NEVER include the gi_* template stack: the CPU and GPU
# paths share only generated topology tables and a byte-exact output contract
# with include/gi_discrete_gradient_labeling.h.

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
    common/gpu_device.cu
    cls/cls2d_kernel.cu
    labeling2d/labeling2d.cu
    labeling2d/label_ctx2d.cu
    cls/cls2d_tiled.cu
    labeling2d/successor2d.cu
    labeling2d/lines2d.cu
    common/gpu_transfer.cu
    dmt3d/dmt3d_gradient.cu
    cls/cls3d_kernel.cu
    cls/cls3d_tiled.cu
    dmt3d/dmt3d_tiled.cu
    api/dgrad_gpu_api.h
    cls/cls2d_tables.h
    api/tile2d.h
    cls/cls2d_internal.cuh
    common/gpu_transfer.h
    common/frontier2d.h
    labeling2d/numeric2d.cuh
    labeling2d/lines2d.cuh
    api/tile3d.h
    cls/cls3d_tables.h
    cls/cls3d_internal.cuh
)
# api/ is the public include folder ("dgrad_gpu_api.h"). The library root is on
# the path too, so the other headers are included by subfolder, for example
# "cls/cls2d_tables.h".
target_include_directories(gpu_dgrad PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/api
    ${CMAKE_CURRENT_SOURCE_DIR})

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
add_executable(gpu_dgrad_smoke tests/smoke_test.cu)
target_link_libraries(gpu_dgrad_smoke PRIVATE gpu_dgrad)

# Stage 1 validation harness: CPU 2D reference vs GPU kernel, bit-exact.
# Host-only C++ (compiled by MSVC, not nvcc) - it may include the gi_* stack.
add_executable(cls2d_validate tests/cls2d_validate.cxx)
target_link_libraries(cls2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(cls2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A0 gate: tiled == in-core == CPU, byte-exact, plus the constant-
# footprint run on a procedural image larger than device memory.
add_executable(cls2d_tiled_validate tests/cls2d_tiled_validate.cxx)
target_link_libraries(cls2d_tiled_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(cls2d_tiled_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A1 gate: tiled successor resolver == in-core flow terminals.
add_executable(succ2d_validate tests/succ2d_validate.cxx)
target_link_libraries(succ2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(succ2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A2 gate: tiled numeric region labeling (acyclicity, decrease,
# tiled == in-core, agreement vs host-integrated basins).
add_executable(numeric2d_validate tests/numeric2d_validate.cxx)
target_link_libraries(numeric2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(numeric2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A4 gate: tiled numeric separatrix lines.
add_executable(lines2d_validate tests/lines2d_validate.cxx)
target_link_libraries(lines2d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(lines2d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# Stage A5 gate: constrained fused kernel == CPU constrained MyRobinsNoalloc.
add_executable(cls2d_constrained_validate tests/cls2d_constrained_validate.cxx)
target_link_libraries(cls2d_constrained_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(cls2d_constrained_validate PRIVATE OpenMP::OpenMP_CXX)
endif()

# 3D discrete Morse theory rules (PLS, PLS-BR, DMF): GPU drivers ==
# CPU raw drivers of include/gi_dmt_*.h, byte-exact, plus timings.
add_executable(dmt3d_validate tests/dmt3d_validate.cxx)
target_link_libraries(dmt3d_validate PRIVATE GInt gpu_dgrad)
target_compile_features(dmt3d_validate PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(dmt3d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(dmt3d_validate PRIVATE /bigobj)
endif()

# 3D lower-star (CLS) gradient: GPU (in-core and tiled, unconstrained and
# constrained) == CPU SlidingWindowRobinsNoalloc::ComputePairing_sliding /
# MyRobinsNoalloc, byte-exact; tiled == in-core; fixed-footprint large run.
add_executable(cls3d_validate tests/cls3d_validate.cxx)
target_link_libraries(cls3d_validate PRIVATE GInt gpu_dgrad)
if(OpenMP_FOUND)
    target_link_libraries(cls3d_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(cls3d_validate PRIVATE /bigobj)
endif()

# Tiled DMT drivers (multi-pass DMF, PLS, PLS-BR): tiled == in-core == CPU,
# halo sweep, fixed-footprint large run.
add_executable(dmt3d_tiled_validate tests/dmt3d_tiled_validate.cxx)
target_link_libraries(dmt3d_tiled_validate PRIVATE GInt gpu_dgrad)
target_compile_features(dmt3d_tiled_validate PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(dmt3d_tiled_validate PRIVATE OpenMP::OpenMP_CXX)
endif()
if(MSVC)
    target_compile_options(dmt3d_tiled_validate PRIVATE /bigobj)
endif()

# GPU-only benchmark of the tiled 3D gradient paths (Robins, PLS, PLS-BR, DMF).
add_executable(gpu3d_bench tests/gpu3d_bench.cxx)
target_link_libraries(gpu3d_bench PRIVATE GInt gpu_dgrad)
target_compile_features(gpu3d_bench PRIVATE cxx_std_17)
if(OpenMP_FOUND)
    target_link_libraries(gpu3d_bench PRIVATE OpenMP::OpenMP_CXX)
endif()

# The Stage 4 parity test (tests/msc2d_gpu_parity.cxx) is added from the ROOT
# CMakeLists after all subdirectories, because it needs msc_2d_lib, which is
# configured after this directory.
