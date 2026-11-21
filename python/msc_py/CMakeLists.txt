# msc_py Python module

# Gather source and header files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add pybind11 module
pybind11_add_module(msc_py MODULE ${SOURCES} ${HEADERS})

# Set output directory
set_target_properties(msc_py PROPERTIES
    OUTPUT_NAME "msc_py"
    LIBRARY_OUTPUT_NAME "msc_py"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Link with necessary libraries
target_link_libraries(msc_py
    PUBLIC GInt
    PUBLIC msc_2d_lib
    PRIVATE pybind11::module
)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(msc_py PUBLIC OpenMP::OpenMP_CXX)
endif()
