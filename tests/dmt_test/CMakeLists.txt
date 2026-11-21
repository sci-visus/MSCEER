# Discrete Morse theory kernel tests.
#
# Each test is a standalone executable whose exit code is 0 on success, the
# same convention as libs/msc_3d_lib/*_smoke and gpu_dgrad/*_validate.

# Every include/dmt_core header, compiled with ONLY include/dmt_core on the
# include path. The directory-level include_directories() from
# cmake/Common.cmake is replaced here, so a core header that includes a gi_*
# header fails to compile.
add_executable(dmt_core_selfcontained dmt_core_selfcontained.cxx)
set_property(TARGET dmt_core_selfcontained PROPERTY
    INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/dmt_core")
target_compile_features(dmt_core_selfcontained PRIVATE cxx_std_17)
set_target_properties(dmt_core_selfcontained PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Host tests link GInt (static offset tables of the topological grid).
foreach(t loopcheck_smoke dmt_tables_validate dmt_smoke dmt_morse_smoke dmt_cycles_smoke gradient_methods_smoke)
    add_executable(${t} ${t}.cxx)
    target_link_libraries(${t} PRIVATE GInt)
    if(OpenMP_FOUND)
        target_link_libraries(${t} PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_compile_features(${t} PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(${t} PRIVATE /bigobj)
    endif()
    set_target_properties(${t} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endforeach()

# The same core compiled by nvcc and run on the GPU, compared byte for byte with
# the host-compiled core. Built only when the CUDA language is enabled.
if(GPU_DGRAD_ENABLED)
    add_executable(dmt_core_nvcc_check dmt_core_nvcc_check.cu)
    set_property(TARGET dmt_core_nvcc_check PROPERTY
        INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/dmt_core")
    set_target_properties(dmt_core_nvcc_check PROPERTIES
        CUDA_STANDARD 17
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()
