# msc_test executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(msc_test ${SOURCES} ${HEADERS})
else()
    add_executable(msc_test ${SOURCES})
endif()

# Set output directory
set_target_properties(msc_test PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

# Link with GInt library
target_link_libraries(msc_test PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(msc_test PUBLIC OpenMP::OpenMP_CXX)
endif()
