CMAKE_MINIMUM_REQUIRED(VERSION 3.12)



FIND_PACKAGE( OpenMP REQUIRED )

ADD_EXECUTABLE(OnDemandAccurate2D ondemandaccurate2d.cxx)
  
if(OpenMP_CXX_FOUND)
    target_link_libraries(OnDemandAccurate2D PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(OnDemandAccurate2D PUBLIC msceer)
target_compile_features(OnDemandAccurate2D PUBLIC cxx_std_11)
set_target_properties(OnDemandAccurate2D PROPERTIES CXX_EXTENSIONS OFF)

                              