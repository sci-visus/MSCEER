# generateisosurface executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(generateisosurface ${SOURCES} ${HEADERS})
else()
    add_executable(generateisosurface ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(generateisosurface PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(generateisosurface PUBLIC OpenMP::OpenMP_CXX)
endif()
