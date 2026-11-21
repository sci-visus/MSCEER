# CMakeLists.txt for test_convolution
# 
# Build options:
#   ENABLE_OPENMP  - Enable OpenMP for CPU parallelization (default: ON)
#   ENABLE_SYCL    - Enable SYCL for GPU offloading (default: OFF)
#   ENABLE_STDPAR  - Enable C++17 stdpar for GPU offloading (default: OFF)
#   ENABLE_DEBUG   - Enable debug output (GINT_DEBUG) (default: OFF)
#
# Example build:
#   mkdir build && cd build
#   cmake -DENABLE_OPENMP=ON -DENABLE_SYCL=OFF ..
#   cmake --build .

cmake_minimum_required(VERSION 3.14)

# Set compiler early if SYCL is requested and compiler is specified
# This allows us to use Intel DPC++ compiler before project() call
if(DEFINED ENV{ENABLE_SYCL_EARLY} OR DEFINED CMAKE_CXX_COMPILER)
    # If CMAKE_CXX_COMPILER is already set (e.g., from command line), use it
    if(CMAKE_CXX_COMPILER)
        string(TOLOWER "${CMAKE_CXX_COMPILER}" COMPILER_CHECK)
        if(COMPILER_CHECK MATCHES "icpx" OR COMPILER_CHECK MATCHES "dpcpp")
            set(CMAKE_CXX_COMPILER "${CMAKE_CXX_COMPILER}" CACHE STRING "C++ compiler" FORCE)
            # Force CMake to recognize Intel DPC++ as Clang-based
            set(CMAKE_CXX_COMPILER_ID "Clang" CACHE STRING "C++ compiler ID" FORCE)
            set(CMAKE_CXX_COMPILER_FORCED TRUE)
        endif()
    endif()
endif()

project(test_convolution CXX)

# C++17 required for stdpar and modern features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build options
option(ENABLE_OPENMP "Enable OpenMP for CPU parallelization" ON)
option(ENABLE_SYCL "Enable SYCL for GPU offloading" OFF)
option(ENABLE_STDPAR "Enable C++17 stdpar for GPU offloading" OFF)
option(ENABLE_DEBUG "Enable debug output (GINT_DEBUG)" OFF)
option(PREFER_ADAPTIVECPP "Prefer AdaptiveCpp over Intel DPC++ when both are available" OFF)

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

# Source files
set(SOURCES
    test_convolution.cxx
)

# Headers (for IDE support)
set(HEADERS
    block_decomposition.h
    convolution_kernel.h
    execution_strategy.h
    execution_strategy_cpu.h
    execution_strategy_sycl.h
    execution_strategy_stdpar.h
)

# Create executable
add_executable(test_convolution ${SOURCES} ${HEADERS})

# Compile definitions
if(ENABLE_DEBUG)
    target_compile_definitions(test_convolution PRIVATE GINT_DEBUG)
endif()

# OpenMP support
if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(test_convolution PUBLIC OpenMP::OpenMP_CXX)
        message(STATUS "OpenMP enabled")
    else()
        message(WARNING "OpenMP requested but not found")
    endif()
endif()

# SYCL support
if(ENABLE_SYCL)
    message(STATUS "=== Configuring SYCL Support ===")
    
    set(SYCL_IMPLEMENTATION_FOUND FALSE)
    
    # Try AdaptiveCpp first (if preferred or if no other SYCL found)
    if(PREFER_ADAPTIVECPP OR NOT SYCL_IMPLEMENTATION_FOUND)
        # AdaptiveCpp requires CMake 3.22+ and uses CONFIG mode
        if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.22")
            find_package(AdaptiveCpp CONFIG QUIET)
            if(AdaptiveCpp_FOUND)
                message(STATUS "Found AdaptiveCpp")
                target_compile_definitions(test_convolution PRIVATE GINT_USE_SYCL)
                
                # AdaptiveCpp uses add_sycl_to_target to mark sources for SYCL compilation
                # This must be called after the target is created
                add_sycl_to_target(
                    TARGET test_convolution
                    SOURCES ${SOURCES}
                )
                
                message(STATUS "SYCL enabled via AdaptiveCpp")
                message(STATUS "  AdaptiveCpp version: ${AdaptiveCpp_VERSION}")
                set(SYCL_IMPLEMENTATION_FOUND TRUE)
            endif()
        else()
            message(STATUS "AdaptiveCpp requires CMake 3.22+, current version: ${CMAKE_VERSION}")
        endif()
    endif()
    
    # Try to find SYCL via find_package (works with Intel oneAPI, hipSYCL, etc.)
    if(NOT SYCL_IMPLEMENTATION_FOUND)
        find_package(SYCL QUIET)
        
        if(SYCL_FOUND)
            target_compile_definitions(test_convolution PRIVATE GINT_USE_SYCL)
            target_link_libraries(test_convolution PUBLIC SYCL::SYCL)
            message(STATUS "SYCL enabled via find_package")
            set(SYCL_IMPLEMENTATION_FOUND TRUE)
        endif()
    endif()
    
    if(NOT SYCL_IMPLEMENTATION_FOUND)
        # Manual detection for Intel oneAPI DPC++
        message(STATUS "SYCL package not found, attempting manual detection...")
        
        # Check if we're already using a SYCL-capable compiler
        set(SYCL_COMPILER_FOUND FALSE)
        
        # Check for Intel DPC++ (icpx or dpcpp)
        # Use execute_process to get Windows Program Files paths (avoids CMake variable name issues)
        if(WIN32)
            execute_process(
                COMMAND powershell -Command "[Environment]::GetEnvironmentVariable('ProgramFiles(x86)')"
                OUTPUT_VARIABLE PROGRAM_FILES_X86
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            execute_process(
                COMMAND powershell -Command "[Environment]::GetEnvironmentVariable('ProgramFiles')"
                OUTPUT_VARIABLE PROGRAM_FILES
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            if(NOT PROGRAM_FILES_X86)
                set(PROGRAM_FILES_X86 "C:/Program Files (x86)")
            endif()
            if(NOT PROGRAM_FILES)
                set(PROGRAM_FILES "C:/Program Files")
            endif()
        endif()
        
        # First check if compiler is already set to a SYCL compiler
        # Note: CMAKE_CXX_COMPILER might be cached from previous runs
        get_property(CACHED_COMPILER CACHE CMAKE_CXX_COMPILER PROPERTY VALUE)
        if(CACHED_COMPILER)
            string(TOLOWER "${CACHED_COMPILER}" COMPILER_LOWER)
            if(COMPILER_LOWER MATCHES "icpx" OR COMPILER_LOWER MATCHES "dpcpp")
                set(ICPX_COMPILER "${CACHED_COMPILER}")
                message(STATUS "Compiler already set to SYCL compiler: ${CACHED_COMPILER}")
            endif()
        endif()
        # Also check CMAKE_CXX_COMPILER variable (after project() call)
        if(CMAKE_CXX_COMPILER AND NOT ICPX_COMPILER)
            string(TOLOWER "${CMAKE_CXX_COMPILER}" COMPILER_LOWER)
            if(COMPILER_LOWER MATCHES "icpx" OR COMPILER_LOWER MATCHES "dpcpp")
                set(ICPX_COMPILER "${CMAKE_CXX_COMPILER}")
                message(STATUS "Compiler variable indicates SYCL compiler: ${CMAKE_CXX_COMPILER}")
            endif()
        endif()
        
        # If not found, search for it
        if(NOT ICPX_COMPILER)
            find_program(ICPX_COMPILER 
                NAMES icpx dpcpp icpx.exe dpcpp.exe
                PATHS
                    "$ENV{ONEAPI_ROOT}/compiler/latest/windows/bin"
                    "$ENV{ONEAPI_ROOT}/compiler/latest/linux/bin"
                    "$ENV{ONEAPI_ROOT}/compiler/latest/bin/compiler"
                    "$ENV{ONEAPI_ROOT}/compiler/*/bin/compiler"
                    "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/windows/bin"
                    "C:/Program Files/Intel/oneAPI/compiler/latest/windows/bin"
                DOC "Intel DPC++ compiler"
            )
            
            # Try Windows-specific paths if not found, including version-specific paths
            if(NOT ICPX_COMPILER AND WIN32)
                # Search for version-specific directories (2025.3, 2024.0, etc.)
                file(GLOB COMPILER_VERSIONS "${PROGRAM_FILES_X86}/Intel/oneAPI/compiler/*")
                file(GLOB COMPILER_VERSIONS2 "${PROGRAM_FILES}/Intel/oneAPI/compiler/*")
                list(APPEND COMPILER_VERSIONS ${COMPILER_VERSIONS2})
                
                foreach(COMPILER_DIR ${COMPILER_VERSIONS})
                    if(IS_DIRECTORY "${COMPILER_DIR}")
                        find_program(ICPX_COMPILER 
                            NAMES icpx dpcpp icpx.exe dpcpp.exe
                            PATHS
                                "${COMPILER_DIR}/bin/compiler"
                                "${COMPILER_DIR}/windows/bin"
                            NO_DEFAULT_PATH
                        )
                        if(ICPX_COMPILER)
                            break()
                        endif()
                    endif()
                endforeach()
            endif()
        endif()
        
        if(ICPX_COMPILER)
            message(STATUS "Found Intel DPC++ compiler: ${ICPX_COMPILER}")
            set(CMAKE_CXX_COMPILER "${ICPX_COMPILER}" CACHE FILEPATH "CXX compiler" FORCE)
            # Force CMake to re-detect the compiler
            unset(CMAKE_CXX_COMPILER_ID CACHE)
            set(SYCL_COMPILER_FOUND TRUE)
        else()
            # Check if current compiler supports SYCL
            # Also check by looking at the compiler path/name
            if(CMAKE_CXX_COMPILER)
                string(TOLOWER "${CMAKE_CXX_COMPILER}" COMPILER_LOWER)
                if(COMPILER_LOWER MATCHES "icpx" OR COMPILER_LOWER MATCHES "dpcpp")
                    set(SYCL_COMPILER_FOUND TRUE)
                    message(STATUS "Compiler path indicates SYCL support: ${CMAKE_CXX_COMPILER}")
                elseif(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" OR 
                       CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
                    # Try to detect if it's a SYCL-capable Clang
                    execute_process(
                        COMMAND ${CMAKE_CXX_COMPILER} --version
                        OUTPUT_VARIABLE COMPILER_VERSION
                        ERROR_VARIABLE COMPILER_ERROR
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                    )
                    if(COMPILER_VERSION MATCHES "DPC" OR COMPILER_VERSION MATCHES "oneAPI" OR COMPILER_VERSION MATCHES "SYCL")
                        set(SYCL_COMPILER_FOUND TRUE)
                        message(STATUS "Current compiler appears to support SYCL: ${CMAKE_CXX_COMPILER}")
                    endif()
                endif()
            endif()
        endif()
        
        if(SYCL_COMPILER_FOUND)
            target_compile_definitions(test_convolution PRIVATE GINT_USE_SYCL)
            
            # Intel DPC++ compiler flags for GPU offloading
            # Use JIT compilation by default (no AOT flags) - works without ocloc tool
            # For AOT compilation, uncomment the -fsycl-targets and -Xsycl-target-backend lines
            target_compile_options(test_convolution PRIVATE 
                -fsycl
                # -fsycl-targets=spir64_gen  # AOT: For Intel GPU (requires ocloc tool)
                # -Xsycl-target-backend=spir64_gen "-device gen12"  # AOT: For Intel integrated/discrete GPUs
            )
            
            # Link options
            target_link_options(test_convolution PRIVATE 
                -fsycl
                # -fsycl-targets=spir64_gen  # AOT: Uncomment if using AOT compilation
                # -Xsycl-target-backend=spir64_gen "-device gen12"  # AOT: Uncomment if using AOT compilation
            )
            
            message(STATUS "SYCL enabled with Intel DPC++ for GPU offloading (JIT mode)")
            message(STATUS "  Compiler: ${CMAKE_CXX_COMPILER}")
            message(STATUS "  Mode: Just-In-Time (JIT) compilation")
            message(STATUS "  Note: For AOT compilation, enable -fsycl-targets flags (requires ocloc tool)")
        else()
            # Try hipSYCL (alternative SYCL implementation)
            find_package(hipSYCL QUIET)
            if(hipSYCL_FOUND)
                target_compile_definitions(test_convolution PRIVATE GINT_USE_SYCL)
                target_compile_options(test_convolution PRIVATE ${hipSYCL_CXX_FLAGS})
                target_link_libraries(test_convolution PRIVATE ${hipSYCL_LIBRARIES})
                message(STATUS "SYCL enabled via hipSYCL")
            else()
                message(WARNING "SYCL requested but no SYCL implementation found.")
                message(WARNING "  To enable SYCL, either:")
                message(WARNING "    1. Install AdaptiveCpp and set CMAKE_PREFIX_PATH or ACPP_ROOT")
                message(WARNING "    2. Install Intel oneAPI Base Toolkit")
                message(WARNING "    3. Set CMAKE_CXX_COMPILER to icpx or dpcpp")
                message(WARNING "    4. Set ONEAPI_ROOT environment variable")
                message(WARNING "    5. Install hipSYCL or another SYCL implementation")
                message(WARNING "  SYCL execution will not be available at runtime.")
            endif()
        endif()
    endif()
endif()

# C++17 stdpar support
if(ENABLE_STDPAR)
    target_compile_definitions(test_convolution PRIVATE GINT_USE_STDPAR)
    
    # MSVC: Built-in support via Concurrency Runtime - no TBB needed
    # GCC: Typically needs TBB as backend
    # Intel oneAPI: Has its own backend
    # NVIDIA HPC SDK: Use -stdpar=gpu flag
    
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        message(STATUS "Stdpar enabled (MSVC built-in Concurrency Runtime)")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # GCC typically needs TBB for parallel algorithms
        find_package(TBB QUIET)
        if(TBB_FOUND)
            target_link_libraries(test_convolution PRIVATE TBB::tbb)
            message(STATUS "Stdpar enabled with TBB backend")
        else()
            message(WARNING "GCC stdpar typically requires TBB. Install TBB or parallel algorithms may run sequentially.")
        endif()
    else()
        message(STATUS "Stdpar enabled (backend depends on compiler)")
    endif()
endif()

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(test_convolution PRIVATE -Wall -Wextra -pedantic)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(test_convolution PRIVATE -O3 -march=native)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(test_convolution PRIVATE /W4)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(test_convolution PRIVATE /O2)
    endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(test_convolution PRIVATE -Wall -Wextra -pedantic)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        target_compile_options(test_convolution PRIVATE -O3 -march=native)
    endif()
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "=== test_convolution Configuration ===")
message(STATUS "  OpenMP: ${ENABLE_OPENMP}")
message(STATUS "  SYCL:   ${ENABLE_SYCL}")
message(STATUS "  Stdpar: ${ENABLE_STDPAR}")
message(STATUS "  Debug:  ${ENABLE_DEBUG}")
message(STATUS "")
