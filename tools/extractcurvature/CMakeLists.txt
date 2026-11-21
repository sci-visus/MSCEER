# extractcurvature executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(extractcurvature ${SOURCES} ${HEADERS})
else()
    add_executable(extractcurvature ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(extractcurvature PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(extractcurvature PUBLIC OpenMP::OpenMP_CXX)
endif()
