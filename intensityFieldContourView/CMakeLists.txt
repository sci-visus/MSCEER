# intensityFieldContourView executable (with VTK support)

# VTK support with specific components
if(VTK_SUPPORT_ENABLED)
    find_package(VTK REQUIRED)
    find_package(VTK COMPONENTS
        vtkChartsCore
        vtkCommonColor
        vtkCommonComputationalGeometry
        vtkCommonCore
        vtkCommonDataModel
        vtkCommonExecutionModel
        vtkCommonMath
        vtkCommonMisc
        vtkCommonSystem
        vtkCommonTransforms
        vtkDICOMParser
        vtkDomainsChemistry
        vtkFiltersAMR
        vtkFiltersCore
        vtkFiltersExtraction
        vtkFiltersFlowPaths
        vtkFiltersGeneral
        vtkFiltersGeneric
        vtkFiltersGeometry
        vtkFiltersHybrid
        vtkFiltersHyperTree
        vtkFiltersImaging
        vtkFiltersModeling
        vtkFiltersParallel
        vtkFiltersParallelImaging
        vtkFiltersPoints
        vtkFiltersProgrammable
        vtkFiltersSMP
        vtkFiltersSelection
        vtkFiltersSources
        vtkFiltersStatistics
        vtkFiltersTexture
        vtkFiltersTopology
        vtkFiltersVerdict
        vtkGeovisCore
        vtkIOAMR
        vtkIOCore
        vtkIOEnSight
        vtkIOExodus
        vtkIOExport
        vtkIOGeometry
        vtkIOImage
        vtkIOImport
        vtkIOInfovis
        vtkIOLSDyna
        vtkIOLegacy
        vtkIOMINC
        vtkIOMovie
        vtkIOPLY
        vtkIOParallel
        vtkIOParallelXML
        vtkIOSQL
        vtkIOSegY
        vtkIOTecplotTable
        vtkIOVideo
        vtkIOXML
        vtkIOXMLParser
        vtkImagingColor
        vtkImagingCore
        vtkImagingFourier
        vtkImagingGeneral
        vtkImagingHybrid
        vtkImagingMath
        vtkImagingMorphological
        vtkImagingSources
        vtkImagingStatistics
        vtkImagingStencil
        vtkInfovisCore
        vtkInfovisLayout
        vtkInteractionImage
        vtkInteractionStyle
        vtkInteractionWidgets
        vtkParallelCore
        vtkRenderingAnnotation
        vtkRenderingContext2D
        vtkRenderingContextOpenGL2
        vtkRenderingCore
        vtkRenderingFreeType
        vtkRenderingGL2PSOpenGL2
        vtkRenderingImage
        vtkRenderingLOD
        vtkRenderingLabel
        vtkRenderingOpenGL2
        vtkRenderingSceneGraph
        vtkRenderingVolume
        vtkRenderingVolumeOpenGL2
        vtkTestingCore
        vtkTestingRendering
        vtkViewsContext2D
        vtkViewsCore
        vtkViewsInfovis
        vtkWrappingTools
        vtkdiy2
        vtkexodusII
        vtkexpat
        vtkfreetype
        vtkgl2ps
        vtkglew
        vtkhdf5
        vtkjpeg
        vtkjsoncpp
        vtkkwiml
        vtklibharu
        vtklibxml2
        vtklz4
        vtknetcdf
        vtkpng
        vtksqlite
        vtksys
        vtktiff
        vtkverdict
        vtkzlib
        QUIET)
    if(NOT VTK_FOUND)
        message("Skipping ImageWriter: ${VTK_NOT_FOUND_MESSAGE}")
        return()
    endif()
    message(STATUS "VTK_VERSION: ${VTK_VERSION}")
    if(VTK_VERSION VERSION_LESS "8.90.0")
        # old system
        include(${VTK_USE_FILE})
    endif()
endif()

# Gather source files
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Create executable
if(HEADERS)
    add_executable(intensityFieldContourView ${SOURCES} ${HEADERS})
else()
    add_executable(intensityFieldContourView ${SOURCES})
endif()

# Link with GInt library
target_link_libraries(intensityFieldContourView PUBLIC GInt)

# Link OpenMP if available
if(OpenMP_FOUND)
    target_link_libraries(intensityFieldContourView PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link VTK if enabled
if(VTK_SUPPORT_ENABLED)
    if(TARGET VTK::VTK)
        target_link_libraries(intensityFieldContourView PUBLIC VTK::VTK)
    elseif(VTK_LIBRARIES)
        target_link_libraries(intensityFieldContourView PUBLIC ${VTK_LIBRARIES})
    endif()
    
    message(STATUS "VTK_VERSION: ${VTK_VERSION}")
    if(NOT VTK_VERSION VERSION_LESS "8.90.0")
        message("auto init vtk modules")
        vtk_module_autoinit(TARGETS intensityFieldContourView GInt MODULES ${VTK_LIBRARIES})
    endif()
endif()
