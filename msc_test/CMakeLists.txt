CMAKE_MINIMUM_REQUIRED(VERSION 3.14)

SET(CMAKE_CXX_STANDARD 11)

ENABLE_LANGUAGE(CXX)


if(APPLE)
  SET(OpenMP_C_FLAGS "")
  SET(OpenMP_CXX_FLAGS "")
  SET(OpenMP_EXE_LINKER_FLAGS "-fopenmp")
else()
  FIND_PACKAGE( OpenMP REQUIRED )
endif()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

FILE(GLOB INCLUDES ${GInt_SOURCE_DIR}/include/*.h *.h)

SET(HEADERS ${INCLUDES})
INCLUDE_DIRECTORIES(${GInt_SOURCE_DIR}/include/)

FILE(GLOB SOURCES ${GInt_SOURCE_DIR}/msc_test/*.cxx)


ADD_EXECUTABLE(msc_test ${SOURCES} ${HEADERS})

# Set the output name for the shared library
set_target_properties(msc_test PROPERTIES OUTPUT_NAME "msc_test")
set_target_properties(msc_test PROPERTIES LIBRARY_OUTPUT_NAME "msc_test")
set_target_properties(msc_test PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set_target_properties(msc_test PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
# Specify any additional libraries that the module depends on

  
TARGET_LINK_LIBRARIES(msc_test LINK_PUBLIC GInt)


                             
 
