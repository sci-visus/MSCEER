CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_STANDARD 11)

ENABLE_LANGUAGE(CXX)

if(VTK_SUPPORT_ENABLED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()

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

FILE(GLOB SOURCES ${GInt_SOURCE_DIR}/extractfilteredlinesoup/*.cxx)

ADD_EXECUTABLE(extractfilteredlinesoup ${HEADERS} ${SOURCES})
  
TARGET_LINK_LIBRARIES(extractfilteredlinesoup LINK_PUBLIC GInt)

if(APPLE)
  set_target_properties(extractfilteredlinesoup PROPERTIES XCODE_ATTRIBUTE_CC "${APPLE_CLANG_COMPILER}")
  set_target_properties(extractfilteredlinesoup PROPERTIES XCODE_ATTRIBUTE_CXX "${APPLE_CLANG_COMPILER}")
endif()
                             
 
