#CMakeLists.txt
cmake_minimum_required (VERSION 3.8)

set(CMAKE_VERBOSE_MAKEFILE 1)

#add_executable(Tutorial tutorial.cpp)
set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/fibonacci.h)
add_executable(Fibonacci fibonacci.cpp ${HEADER_FILES})
target_link_libraries(Fibonacci PRIVATE testlibA)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

FIND_PACKAGE(PythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(CMAKE_SWIG_FLAGS "-Wall")

SET(SWIG_MODULE_fibonacci_EXTRA_DEPS 
  Fibonacci
  testlibA
)
SET_SOURCE_FILES_PROPERTIES(fibonacci.i PROPERTIES CPLUSPLUS ON)
SET_PROPERTY(SOURCE fibonacci.i PROPERTY SWIG_FLAGS "-builtin")

SWIG_ADD_LIBRARY(fibonacci TYPE MODULE LANGUAGE python SOURCES fibonacci.i fibonacci.cpp)
SWIG_LINK_LIBRARIES(fibonacci ${PYTHON_LIBRARIES} testlibA)


install(TARGETS ${SWIG_MODULE_fibonacci_REAL_NAME} DESTINATION test_install)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fibonacci.py DESTINATION test_install) 