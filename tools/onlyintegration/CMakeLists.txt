# OnlyIntegration executable

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(OnlyIntegration ${SOURCES} ${HEADERS})
else()
    add_executable(OnlyIntegration ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(OnlyIntegration PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(OnlyIntegration PUBLIC OpenMP::OpenMP_CXX)
endif()
