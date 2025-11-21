# extraccentroids executable (with VTK support)

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(extraccentroids ${SOURCES} ${HEADERS})
else()
    add_executable(extraccentroids ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(extraccentroids PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(extraccentroids PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link VTK if enabled
if(VTK_SUPPORT_ENABLED)
    if(TARGET VTK::VTK)
        target_link_libraries(extraccentroids PUBLIC VTK::VTK)
    elseif(VTK_LIBRARIES)
        target_link_libraries(extraccentroids PUBLIC ${VTK_LIBRARIES})
    endif()
endif()
