# Common CMake utilities and configuration for GInt project

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Configure OpenMP
if(APPLE)
    find_package(OpenMP QUIET)
    if(OpenMP_FOUND)
        set(OpenMP_EXE_LINKER_FLAGS "-fopenmp")
    else()
        message(WARNING "OpenMP not found on Apple. Some features may not work correctly.")
    endif()
else()
    find_package(OpenMP REQUIRED)
endif()

# Add OpenMP to compiler flags if available
if(OpenMP_FOUND)
    if(APPLE)
        # On Apple, OpenMP_CXX_FLAGS might be empty, add manually
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    else()
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

# Common include directory
include_directories(${GInt_SOURCE_DIR}/include)

# Function to create a standard executable that links with GInt library
function(add_gint_executable target_name)
    # Get source files from current directory
    file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
    file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
    
    # Add headers to sources if they exist
    if(HEADERS)
        add_executable(${target_name} ${SOURCES} ${HEADERS})
    else()
        add_executable(${target_name} ${SOURCES})
    endif()
    
    # Link with GInt library
    target_link_libraries(${target_name} PUBLIC GInt)
    
    # Link OpenMP if available
    if(OpenMP_FOUND)
        target_link_libraries(${target_name} PUBLIC OpenMP::OpenMP_CXX)
    endif()
    
    # Link VTK if enabled
    if(VTK_SUPPORT_ENABLED AND TARGET VTK::VTK)
        target_link_libraries(${target_name} PUBLIC VTK::VTK)
    elseif(VTK_SUPPORT_ENABLED AND VTK_LIBRARIES)
        target_link_libraries(${target_name} PUBLIC ${VTK_LIBRARIES})
    endif()
endfunction()

