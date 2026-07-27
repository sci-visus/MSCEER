# msc_3d_lib pure C++ library

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_lib_smoke.cxx")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_partitioned_smoke.cxx")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_partitioned_benchmark.cxx")

add_library(msc_3d_lib ${SOURCES} ${HEADERS})

target_include_directories(msc_3d_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(msc_3d_lib
    PUBLIC GInt
)

if(OpenMP_FOUND)
    target_link_libraries(msc_3d_lib PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(msc_3d_lib_smoke "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_lib_smoke.cxx")
target_link_libraries(msc_3d_lib_smoke PUBLIC msc_3d_lib)

set_target_properties(msc_3d_lib_smoke PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

add_executable(msc_3d_partitioned_smoke "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_partitioned_smoke.cxx")
target_link_libraries(msc_3d_partitioned_smoke PUBLIC msc_3d_lib)

set_target_properties(msc_3d_partitioned_smoke PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

add_executable(msc_3d_partitioned_benchmark "${CMAKE_CURRENT_SOURCE_DIR}/msc_3d_partitioned_benchmark.cxx")
target_link_libraries(msc_3d_partitioned_benchmark PUBLIC msc_3d_lib)

set_target_properties(msc_3d_partitioned_benchmark PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
