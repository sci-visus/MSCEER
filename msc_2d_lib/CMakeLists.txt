# msc_2d_lib pure C++ library

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/msc_2d_lib_smoke.cxx")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/msc_2d_partitioned_smoke.cxx")

add_library(msc_2d_lib ${SOURCES} ${HEADERS})

target_include_directories(msc_2d_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

# Require C++17 (see src/CMakeLists.txt); inherited via GInt too, but explicit
# here so msc_2d_lib's own sources build at the right standard regardless.
target_compile_features(msc_2d_lib PUBLIC cxx_std_17)

target_link_libraries(msc_2d_lib
    PUBLIC GInt
)

if(OpenMP_FOUND)
    target_link_libraries(msc_2d_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

# GPU gradient producer (ComputeOptions::useGpuGradient) when built with
# GPU_DGRAD_ENABLED; without it the option falls back to the CPU path.
if(TARGET gpu_dgrad)
    target_link_libraries(msc_2d_lib PUBLIC gpu_dgrad)
    target_compile_definitions(msc_2d_lib PRIVATE MSC2D_HAS_GPU_DGRAD)
endif()

add_executable(msc_2d_lib_smoke "${CMAKE_CURRENT_SOURCE_DIR}/msc_2d_lib_smoke.cxx")
target_link_libraries(msc_2d_lib_smoke PUBLIC msc_2d_lib)

set_target_properties(msc_2d_lib_smoke PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

add_executable(msc_2d_partitioned_smoke "${CMAKE_CURRENT_SOURCE_DIR}/msc_2d_partitioned_smoke.cxx")
target_link_libraries(msc_2d_partitioned_smoke PUBLIC msc_2d_lib)

set_target_properties(msc_2d_partitioned_smoke PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
