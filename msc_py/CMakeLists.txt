CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

SET(CMAKE_CXX_STANDARD 14)  # Or 17 if possible
ENABLE_LANGUAGE(CXX)

if(APPLE)
  find_package(OpenMP)
  if(OpenMP_FOUND)
    SET(OpenMP_EXE_LINKER_FLAGS "-fopenmp")
  else()
    message(WARNING "OpenMP not found! Compile will proceed without OpenMP support.")
  endif()
else()
  FIND_PACKAGE(OpenMP REQUIRED)
endif()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

# Gather source and header files
FILE(GLOB INCLUDES ${GInt_SOURCE_DIR}/include/*.h *.h)
FILE(GLOB SOURCES ${GInt_SOURCE_DIR}/msc_py/*.cxx)

SET(HEADERS ${INCLUDES})
INCLUDE_DIRECTORIES(${GInt_SOURCE_DIR}/include/)

# Add pybind11 module
pybind11_add_module(msc_py MODULE ${SOURCES} ${HEADERS})

# Set the output name for the shared library
set_target_properties(msc_py PROPERTIES OUTPUT_NAME "msc_py")
set_target_properties(msc_py PROPERTIES LIBRARY_OUTPUT_NAME "msc_py")
set_target_properties(msc_py PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set_target_properties(msc_py PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Link with necessary libraries
TARGET_LINK_LIBRARIES(msc_py LINK_PUBLIC GInt)
target_link_libraries(msc_py PRIVATE pybind11::module)  # Consider using pybind11::embed only if embedding Python in C++
