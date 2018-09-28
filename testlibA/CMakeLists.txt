#CMakeLists.txt
cmake_minimum_required (VERSION 3.8)

message(STATUS "current dir" ${CMAKE_CURRENT_SOURCE_DIR})

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/testlibA.h)

add_library(testlibA STATIC 
	src/testlibA.cpp ${HEADER_FILES}
)


target_include_directories(testlibA 
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

install(TARGETS testlibA DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION include)
