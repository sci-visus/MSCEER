# extractgraph executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(extractgraph ${SOURCES} ${HEADERS})
else()
    add_executable(extractgraph ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(extractgraph PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(extractgraph PUBLIC OpenMP::OpenMP_CXX)
endif()
