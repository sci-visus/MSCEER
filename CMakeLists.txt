cmake_minimum_required(VERSION 3.14)
project(GInt LANGUAGES CXX)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Export compile commands for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include common utilities
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(Common)

# Build the core library first
add_subdirectory(src)

# VTK support
option(VTK_SUPPORT_ENABLED "Enable VTK support" FALSE)
if(VTK_SUPPORT_ENABLED)
    find_package(VTK REQUIRED)
    add_compile_definitions(VTK_ENABLED)
endif()

# Integration tools
option(INTEGRATION_ENABLED "Build integration" FALSE)
if(INTEGRATION_ENABLED)
    add_subdirectory(integrator)
endif()

option(ONLY_INTEGRATOR_ENABLED "Build only integration" FALSE)
if(ONLY_INTEGRATOR_ENABLED)
    add_subdirectory(onlyintegration)
endif()

option(EXTRACT_MESH_FUNCTION_ENABLED "Build extract mesh function" FALSE)
if(EXTRACT_MESH_FUNCTION_ENABLED)
    add_subdirectory(extractmeshfunction)
endif()

option(ON_DEMAND_ACCURATE_ENABLED "Build ondemand accurate" FALSE)
if(ON_DEMAND_ACCURATE_ENABLED)
    add_subdirectory(ondemandaccurate)
    add_subdirectory(ondemandaccurate2d)
endif()

option(STREMLINES_2D_ENABLED "Build 2d streamline integrator" FALSE)
if(STREMLINES_2D_ENABLED)
    add_subdirectory(streamlines2d)
endif()

option(SPARSEBLOCKEDCOMPUTE_ENABLED "Build sparse blocked compute integrator" FALSE)
if(SPARSEBLOCKEDCOMPUTE_ENABLED)
    add_subdirectory(sparseblockedcompute)
endif()

option(ALL_INTEGRATION_MSC_ENABLED "Build all integration MSC" FALSE)
if(ALL_INTEGRATION_MSC_ENABLED)
    add_subdirectory(allintegrationmorsecomplex)
endif()

option(CONVERGENT_ENABLED "Build convergent" FALSE)
if(CONVERGENT_ENABLED)
    add_subdirectory(convergent)
endif()

# Extraction tools
option(EXTRACT_CENTROIDS_ENABLED "Build extract ligaments centroids" FALSE)
if(EXTRACT_CENTROIDS_ENABLED)
    add_subdirectory(extractcentroids)
endif()

option(ACCURATE_3D_MORSE_LABELING "Compute accurate 3d morse labeling" FALSE)
if(ACCURATE_3D_MORSE_LABELING)
    add_subdirectory(accurate3dmorselabeling)
endif()

option(LABELING_TO_STRATIFIED_COMPLEX "Turn grid labeling into stratified complex" FALSE)
if(LABELING_TO_STRATIFIED_COMPLEX)
    add_subdirectory(labelingtostratifiedcomplex)
endif()

option(PRUNE_LIGAMENTS "Prune away boundary ligaments" FALSE)
if(PRUNE_LIGAMENTS)
    add_subdirectory(pruneligaments)
endif()

option(EXTRACT_GRAPH_ENABLED "Build extract graph" FALSE)
if(EXTRACT_GRAPH_ENABLED)
    add_subdirectory(extractgraph)
endif()

option(EXTRACT_MSC_ENABLED "Build extract msc" FALSE)
if(EXTRACT_MSC_ENABLED)
    add_subdirectory(extractmsc)
endif()

option(EXTRACT_FILTERED_LINE_SOUP_ENABLED "Build extract filtered line soup" FALSE)
if(EXTRACT_FILTERED_LINE_SOUP_ENABLED)
    add_subdirectory(extractfilteredlinesoup)
endif()

option(EXTRACT_METRICS_ENABLED "Build extract metrics" FALSE)
if(EXTRACT_METRICS_ENABLED)
    add_subdirectory(extractmetrics)
endif()

option(EXTRACT_SIMPLIFIED_MORSE "Build extract simplified Morse labeling" FALSE)
if(EXTRACT_SIMPLIFIED_MORSE)
    add_subdirectory(extractsimplifiedmorse)
endif()

option(EXTRACT_CROSS_SECTION_ENABLED "Build extract cross section" FALSE)
if(EXTRACT_CROSS_SECTION_ENABLED)
    add_subdirectory(extractcrosssection)
endif()

option(EXTRACT_MULTI_CROSS_SECTION_ENABLED "Build extract multiple cross sections" FALSE)
if(EXTRACT_MULTI_CROSS_SECTION_ENABLED)
    add_subdirectory(extractmulticrosssections)
endif()

option(EXTRACT_CURVATURE_ENABLED "Build extract curvature" FALSE)
if(EXTRACT_CURVATURE_ENABLED)
    add_subdirectory(extractcurvature)
endif()

option(LABEL_CONN_COMP_ENABLED "Label connected components" FALSE)
if(LABEL_CONN_COMP_ENABLED)
    add_subdirectory(labelconnectedcomponents)
endif()

option(INTENSITY_FIELD_CONTOUR_VIEWER_ENABLED "View Contours of For Various Thresholds" FALSE)
if(INTENSITY_FIELD_CONTOUR_VIEWER_ENABLED)
    add_subdirectory(intensityFieldContourView)
endif()

option(EXTRACT_CONN_COMP_ENABLED "Extract connected components" FALSE)
if(EXTRACT_CONN_COMP_ENABLED)
    add_subdirectory(extractconnectedcomponents)
endif()

option(EXTRACT_LABELED_COMP_ENABLED "Extract labeled components" FALSE)
if(EXTRACT_LABELED_COMP_ENABLED)
    add_subdirectory(extractlabeledcomponents)
endif()

# Python bindings
option(FILTERED_LINESOUP_PY "Compute grad and extract linesoup" FALSE)
if(FILTERED_LINESOUP_PY)
    add_subdirectory(extern/pybind11)
    add_subdirectory(filteredlinesoup_py)
    add_subdirectory(filteredlinesoup)
endif()

option(MSC_PY "Python interface to compute dgrad and msc" FALSE)
if(MSC_PY)
    add_subdirectory(extern/pybind11)
    add_subdirectory(msc_py)
    add_subdirectory(msc_test)
endif()

option(MSC_2D_LIB "Pure C++ 2D Morse-Smale complex library" FALSE)
if(MSC_2D_LIB)
    add_subdirectory(msc_2d_lib)
endif()

# Steepest descent tools
option(STEEPEST_ENABLED "Build steepest" FALSE)
if(STEEPEST_ENABLED)
    add_subdirectory(steepest)
endif()

option(STEEPEST_LSTAR_ENABLED "Build steepest lstar" FALSE)
if(STEEPEST_LSTAR_ENABLED)
    add_subdirectory(steepest_lstar)
endif()

option(STEEPEST_LSTAR_TILED_ENABLED "Build steepest lstar tiled" FALSE)
if(STEEPEST_LSTAR_TILED_ENABLED)
    add_subdirectory(steepest_lstar_tiled)
endif()

option(MAKE_CONST_SMALL_GRID_ENABLED "Build make const small grid" FALSE)
if(MAKE_CONST_SMALL_GRID_ENABLED)
    add_subdirectory(make_const_small_grid)
endif()

option(MORSE_HIERARCHY_2D_ENABLED "Build Morse hierarchy 2D" FALSE)
if(MORSE_HIERARCHY_2D_ENABLED)
    add_subdirectory(morsehierarchy2d)
endif()

option(EXTRACT_2D_RIDGE_GRAPH_ENABLED "Extract 2D geometric ridge graph" FALSE)
if(EXTRACT_2D_RIDGE_GRAPH_ENABLED)
    add_subdirectory(extract2dridgegraph)
endif()

option(EXTRACT_2D_REGION_GRAPH_ENABLED "Extract 2D geometric region graph" FALSE)
if(EXTRACT_2D_REGION_GRAPH_ENABLED)
    add_subdirectory(extract2dregiongraph)
endif()

option(PATH_COMPRESS_INTEGRATOR_ENABLED "Build compress integrator" FALSE)
if(PATH_COMPRESS_INTEGRATOR_ENABLED)
    add_subdirectory(pathcompressintegrator)
endif()

option(MIN_NETWORK_ENABLED "Build min network" FALSE)
if(MIN_NETWORK_ENABLED)
    add_subdirectory(buildminnetwork)
endif()

option(TRIANGLE_TEST_ENABLED "Build triangle test" FALSE)
if(TRIANGLE_TEST_ENABLED)
    add_subdirectory(triangletest)
endif()

option(GENERATE_ISOSURFACE_ENABLED "Generate isosurface util" FALSE)
if(GENERATE_ISOSURFACE_ENABLED)
    add_subdirectory(generateisosurface)
endif()

option(EXPERIMENTAL_ENABLED "Build experimental" FALSE)
if(EXPERIMENTAL_ENABLED)
    add_subdirectory(experimental)
endif()
