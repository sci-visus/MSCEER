# streamlines2d executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(streamlines2d ${SOURCES} ${HEADERS})
else()
    add_executable(streamlines2d ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(streamlines2d PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(streamlines2d PUBLIC OpenMP::OpenMP_CXX)
endif()
