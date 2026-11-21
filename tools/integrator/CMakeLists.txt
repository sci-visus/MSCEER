# GradIntegrator executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(GradIntegrator ${SOURCES} ${HEADERS})
else()
    add_executable(GradIntegrator ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(GradIntegrator PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(GradIntegrator PUBLIC OpenMP::OpenMP_CXX)
endif()
