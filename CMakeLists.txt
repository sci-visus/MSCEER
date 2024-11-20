CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(GInt)

set(CMAKE_BUILD_TYPE Release)

add_definitions(-std=c++11)
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()



ADD_SUBDIRECTORY(src)

option(VTK_SUPPORT_ENABLED "Enable VTK" FALSE)
option(VTK_ENABLED "Enable VTK" FALSE)
if(VTK_SUPPORT_ENABLED)
  find_package(VTK REQUIRED)
  #include(${VTK_USE_FILE})
  ADD_DEFINITIONS(-DVTK_ENABLED)
endif()

if(CMAKE_GENERATOR STREQUAL Xcode)
  set(APPLE_CLANG_COMPILER /usr/local/opt/llvm/bin/clang)
  set(CMAKE_CXX_FLAGS "-fopenmp -std=c++11")
  link_directories(/usr/local/Cellar/llvm/5.0.1/lib)
endif()

#option(FIND_OPENMP_ENABLED "Enable FindOpenMP, not necessary on Apple" OFF)

option(INTEGRATION_ENABLED "Build integration" FALSE)
if(INTEGRATION_ENABLED)
  ADD_SUBDIRECTORY(integrator)
endif()

option(ONLY_INTEGRATOR_ENABLED "Build only integration" FALSE)
if(ONLY_INTEGRATOR_ENABLED)
  ADD_SUBDIRECTORY(onlyintegration)
endif()

option(EXTRACT_MESH_FUNCTION_ENABLED "Build extract mesh function" FALSE)
if(EXTRACT_MESH_FUNCTION_ENABLED)
  ADD_SUBDIRECTORY(extractmeshfunction)
endif()

option(ON_DEMAND_ACCURATE_ENABLED "Build ondemand accurate" FALSE)
if(ON_DEMAND_ACCURATE_ENABLED)
  ADD_SUBDIRECTORY(ondemandaccurate)
  ADD_SUBDIRECTORY(ondemandaccurate2d)
endif()

option(STREMLINES_2D_ENABLED "Build 2d streamline integrator" FALSE)
if(STREMLINES_2D_ENABLED)
  ADD_SUBDIRECTORY(streamlines2d)
endif()

option(SPARSEBLOCKEDCOMPUTE_ENABLED "Build sparse blocked compute integrator" FALSE)
if(SPARSEBLOCKEDCOMPUTE_ENABLED)
  ADD_SUBDIRECTORY(sparseblockedcompute)
endif()

option(ALL_INTEGRATION_MSC_ENABLED "Build all integration MSC" FALSE)
if(ALL_INTEGRATION_MSC_ENABLED)
  ADD_SUBDIRECTORY(allintegrationmorsecomplex)
endif()

option(CONVERGENT_ENABLED "Build convergent" FALSE)
if(CONVERGENT_ENABLED)
  ADD_SUBDIRECTORY(convergent)
endif()

option(EXTRACT_CENTROIDS_ENABLED "Build extract ligaments centroids" FALSE)
if(EXTRACT_CENTROIDS_ENABLED)
  ADD_SUBDIRECTORY(extractcentroids)
endif()

option(ACCURATE_3D_MORSE_LABELING "Compute accurate 3d morse labeling" FALSE)
if(ACCURATE_3D_MORSE_LABELING)
  ADD_SUBDIRECTORY(accurate3dmorselabeling)
endif()

option(LABELING_TO_STRATIFIED_COMPLEX "Turn grid labeling into stratified complex" FALSE)
if(LABELING_TO_STRATIFIED_COMPLEX)
  ADD_SUBDIRECTORY(labelingtostratifiedcomplex)
endif()

option(PRUNE_LIGAMENTS "prune away boundary ligaments" FALSE)
if(PRUNE_LIGAMENTS)
  ADD_SUBDIRECTORY(pruneligaments)
endif()

option(EXTRACT_GRAPH_ENABLED "Build extract graph" FALSE)
if(EXTRACT_GRAPH_ENABLED)
  ADD_SUBDIRECTORY(extractgraph)
endif()

option(EXTRACT_MSC_ENABLED "Build extract msc" FALSE)
if(EXTRACT_MSC_ENABLED)
  ADD_SUBDIRECTORY(extractmsc)
endif()

option(EXTRACT_FILTERED_LINE_SOUP_ENABLED "Build extract filtered line soup" FALSE)
if(EXTRACT_FILTERED_LINE_SOUP_ENABLED)
  ADD_SUBDIRECTORY(extractfilteredlinesoup)
endif()

option(EXTRACT_METRICS_ENABLED "Build extract metrics" FALSE)
if(EXTRACT_METRICS_ENABLED)
  ADD_SUBDIRECTORY(extractmetrics)
endif()

option(EXTRACT_SIMPLIFIED_MORSE "Build extract simplified Morse labeling" FALSE)
if(EXTRACT_SIMPLIFIED_MORSE)
  ADD_SUBDIRECTORY(extractsimplifiedmorse)
endif()

option(EXTRACT_CROSS_SECTION_ENABLED "Build extract cross section" FALSE)
if(EXTRACT_CROSS_SECTION_ENABLED)
  ADD_SUBDIRECTORY(extractcrosssection)
endif()

option(EXTRACT_MULTI_CROSS_SECTION_ENABLED "Build extract multiple cross sections" FALSE)
if(EXTRACT_MULTI_CROSS_SECTION_ENABLED)
  ADD_SUBDIRECTORY(extractmulticrosssections)
endif()

option(EXTRACT_CURVATURE_ENABLED "Build extract curvature" FALSE)
if(EXTRACT_CURVATURE_ENABLED)
  ADD_SUBDIRECTORY(extractcurvature)
endif()

option(LABEL_CONN_COMP_ENABLED "Label connected components" FALSE)
if(LABEL_CONN_COMP_ENABLED)
  ADD_SUBDIRECTORY(labelconnectedcomponents)
endif()

option(INTENSITY_FIELD_CONTOUR_VIEWER_ENABLED "View Contours of For Various Thresholds" FALSE)
if(INTENSITY_FIELD_CONTOUR_VIEWER_ENABLED)
  ADD_SUBDIRECTORY(intensityFieldContourView)
endif()

option(EXTRACT_CONN_COMP_ENABLED "Extract connected components" FALSE)
if(EXTRACT_CONN_COMP_ENABLED)
  ADD_SUBDIRECTORY(extractconnectedcomponents)
endif()

option(EXTRACT_LABELED_COMP_ENABLED "Extract labeled components" FALSE)
if(EXTRACT_LABELED_COMP_ENABLED)
  ADD_SUBDIRECTORY(extractlabeledcomponents)
endif()

option(FILTERED_LINESOUP_PY "Compute grad and extract linesoup" FALSE)
if(FILTERED_LINESOUP_PY)

  add_subdirectory(extern/pybind11)
  ADD_SUBDIRECTORY(filteredlinesoup_py)
  ADD_SUBDIRECTORY(filteredlinesoup)
endif()

option(MSC_PY "Python interface to compute dgrad and msc" FALSE)
if(MSC_PY)
  add_subdirectory(extern/pybind11)
  ADD_SUBDIRECTORY(msc_py)
  add_subdirectory(msc_test)
endif()

option(STEEPEST_ENABLED "Build steepest" FALSE)
if(STEEPEST_ENABLED)
  ADD_SUBDIRECTORY(steepest)
endif()

option(STEEPEST_LSTAR_ENABLED "Build steepest lstar" FALSE)
if(STEEPEST_LSTAR_ENABLED)
  ADD_SUBDIRECTORY(steepest_lstar)
endif()

option(STEEPEST_LSTAR_TILED_ENABLED "Build steepest lstar tiled" FALSE)
if(STEEPEST_LSTAR_TILED_ENABLED)
  ADD_SUBDIRECTORY(steepest_lstar_tiled)
endif()

option(MAKE_CONST_SMALL_GRID_ENABLED "Build make const small grid" FALSE)
if(MAKE_CONST_SMALL_GRID_ENABLED)
  ADD_SUBDIRECTORY(make_const_small_grid)
endif()

option(MORSE_HIERARCHY_2D_ENABLED "Build Morse hierarchy 2D" FALSE)
if(MORSE_HIERARCHY_2D_ENABLED)
  ADD_SUBDIRECTORY(morsehierarchy2d)
endif()

option(EXTRACT_2D_RIDGE_GRAPH_ENABLED "Extract 2D geometric ridge graph" FALSE)
if(EXTRACT_2D_RIDGE_GRAPH_ENABLED)
  ADD_SUBDIRECTORY(extract2dridgegraph)
endif()

option(EXTRACT_2D_REGION_GRAPH_ENABLED "Extract 2D geometric region graph" FALSE)
if(EXTRACT_2D_REGION_GRAPH_ENABLED)
  ADD_SUBDIRECTORY(extract2dregiongraph)
endif()

option(PATH_COMPRESS_INTEGRATOR_ENABLED "Build compress integrator" FALSE)
if(PATH_COMPRESS_INTEGRATOR_ENABLED)
  ADD_SUBDIRECTORY(pathcompressintegrator)
endif()

option(MIN_NETWORK_ENABLED "Build min network" FALSE)
if(MIN_NETWORK_ENABLED)
  ADD_SUBDIRECTORY(buildminnetwork)
endif()

option(TRIANGLE_TEST_ENABLED "Build triangle test" FALSE)
if(TRIANGLE_TEST_ENABLED)
  ADD_SUBDIRECTORY(triangletest)
endif()

option(GENERATE_ISOSURFACE_ENABLED "Generate isosurface util" FALSE)
if(GENERATE_ISOSURFACE_ENABLED)
  ADD_SUBDIRECTORY(generateisosurface)
endif()

option(EXPERIMENTAL_ENABLED "Build experimental" FALSE)
if(EXPERIMENTAL_ENABLED)
  ADD_SUBDIRECTORY(experimental)
endif()  
