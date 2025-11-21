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

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(GInt PUBLIC OpenMP::OpenMP_CXX)
endif()
                     