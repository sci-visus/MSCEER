CMAKE_MINIMUM_REQUIRED(VERSION 3.12)



FIND_PACKAGE( OpenMP REQUIRED )

ADD_EXECUTABLE(OnDemandAccurate3D ondemandaccurate.cxx)
  
if(OpenMP_CXX_FOUND)
    target_link_libraries(OnDemandAccurate3D PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(OnDemandAccurate3D PUBLIC msceer)
target_compile_features(OnDemandAccurate3D PUBLIC cxx_std_11)
set_target_properties(OnDemandAccurate3D PROPERTIES CXX_EXTENSIONS OFF)

                              