CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_STANDARD 11)

ENABLE_LANGUAGE(CXX)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(APPLE)
    SET(OpenMP_C_FLAGS "")
    SET(OpenMP_CXX_FLAGS "")
    #SET(OpenMP_EXE_LINKER_FLAGS "-fopenmp")
else()
    FIND_PACKAGE( OpenMP REQUIRED )
endif()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

FILE(GLOB INCLUDES ${GInt_SOURCE_DIR}/include/*.h)

SET(HEADERS ${INCLUDES})
INCLUDE_DIRECTORIES(${GInt_SOURCE_DIR}/include/)

FILE(GLOB SOURCES ${GInt_SOURCE_DIR}/extractmeshfunction/*.cxx)

ADD_EXECUTABLE(extractmeshfunction ${HEADERS} ${SOURCES})
  
TARGET_LINK_LIBRARIES(extractmeshfunction LINK_PUBLIC GInt)
                             
    