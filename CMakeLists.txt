#CMakeLists.txt
cmake_minimum_required (VERSION 3.12)

project (MSCEER VERSION 1.0
              DESCRIPTION "Morse-Smale Complex Extraction Exploration and Reasonign"
              LANGUAGES CXX)
				  
add_subdirectory(MSCEER)
add_subdirectory(DiscreteGradComputeCL)

add_subdirectory(testlibA)
add_subdirectory(test)