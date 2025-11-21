# sparseblockedcompute executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(sparseblockedcompute ${SOURCES} ${HEADERS})
else()
    add_executable(sparseblockedcompute ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(sparseblockedcompute PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(sparseblockedcompute PUBLIC OpenMP::OpenMP_CXX)
endif()
