# PathCompressGradIntegrator executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(PathCompressGradIntegrator ${SOURCES} ${HEADERS})
else()
    add_executable(PathCompressGradIntegrator ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(PathCompressGradIntegrator PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(PathCompressGradIntegrator PUBLIC OpenMP::OpenMP_CXX)
endif()
