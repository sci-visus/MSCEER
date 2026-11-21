# make_const_small_grid executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(make_const_small_grid ${SOURCES} ${HEADERS})
else()
    add_executable(make_const_small_grid ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(make_const_small_grid PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(make_const_small_grid PUBLIC OpenMP::OpenMP_CXX)
endif()
