# Core GInt library

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")

# Create the library
add_library(GInt ${SOURCES})

# Set public include directories
target_include_directories(GInt PUBLIC
    $<BUILD_INTERFACE:${GInt_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# GInt's headers use C++14+ features (e.g. gi_vectors.h constexpr operator[],
# which is implicitly const under C++11). Require C++17 as a usage requirement
# so both GInt and anything linking it compile at a sufficient standard,
# overriding the C++11 default from cmake/Common.cmake. Portable to MSVC.
target_compile_features(GInt PUBLIC cxx_std_17)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(GInt PUBLIC OpenMP::OpenMP_CXX)
endif()
                     