CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_STANDARD 11)

ENABLE_LANGUAGE(CXX)

if(VTK_SUPPORT_ENABLED)
    find_package(VTK REQUIRED)# NO_MODULE)
    find_package(VTK COMPONENTS
        #incremental
        #txaio
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
        #vtkDomainsChemistryOpenGL2
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
        #vtkFiltersParallelDIY2
        #vtkFiltersParallelFlowPaths
        #vtkFiltersParallelGeometry
        vtkFiltersParallelImaging
        #vtkFiltersParallelMPI
        #vtkFiltersParallelStatistics
        #vtkFiltersParallelVerdict
        vtkFiltersPoints
        vtkFiltersProgrammable
        #vtkFiltersPython
        #vtkFiltersReebGraph
        vtkFiltersSMP
        vtkFiltersSelection
        vtkFiltersSources
        vtkFiltersStatistics
        vtkFiltersTexture
        vtkFiltersTopology
        vtkFiltersVerdict
        #vtkGUISupportQt
        #vtkGUISupportQtOpenGL
        #vtkGUISupportQtSQL
        vtkGeovisCore
        vtkIOAMR
        vtkIOCore
        vtkIOEnSight
        vtkIOExodus
        vtkIOExport
        #vtkIOExportOpenGL2
        #vtkIOFFMPEG
        #vtkIOGDAL
        #vtkIOGeoJSON
        vtkIOGeometry
        vtkIOImage
        vtkIOImport
        vtkIOInfovis
        vtkIOLSDyna
        vtkIOLegacy
        vtkIOMINC
        #vtkIOMPIImage
        #vtkIOMPIParallel
        vtkIOMovie
        #vtkIOMySQL
        #vtkIONetCDF
        #vtkIOODBC
        vtkIOPLY
        vtkIOParallel
        #vtkIOParallelExodus
        #vtkIOParallelLSDyna
        #vtkIOParallelNetCDF
        vtkIOParallelXML
        #vtkIOParallelXdmf3
        #vtkIOPostgreSQL
        vtkIOSQL
        vtkIOSegY
        #vtkIOTRUCHAS
        vtkIOTecplotTable
        #vtkIOVPIC
        vtkIOVideo
        vtkIOXML
        vtkIOXMLParser
        #vtkIOXdmf2
        #vtkIOXdmf3
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
        #vtkInfovisBoost
        #vtkInfovisBoostGraphAlgorithms
        vtkInfovisCore
        vtkInfovisLayout
        vtkInteractionImage
        vtkInteractionStyle
        vtkInteractionWidgets
        #vtkLocalExample
        #vtkMetaIO
        vtkParallelCore
        #vtkParallelMPI
        #vtkParallelMPI4Py
        #vtkParseOGLExt
        #vtkPython
        #vtkPythonInterpreter
        vtkRenderingAnnotation
        vtkRenderingContext2D
        vtkRenderingContextOpenGL2
        vtkRenderingCore
        vtkRenderingFreeType
        #vtkRenderingFreeTypeFontConfig
        vtkRenderingGL2PSOpenGL2
        vtkRenderingImage
        #vtkRenderingLICOpenGL2
        vtkRenderingLOD
        vtkRenderingLabel
        #vtkRenderingMatplotlib
        vtkRenderingOpenGL2
        #vtkRenderingParallel
        #vtkRenderingParallelLIC
        #vtkRenderingQt
        vtkRenderingSceneGraph
        #vtkRenderingTk
        vtkRenderingVolume
        #vtkRenderingVolumeAMR
        vtkRenderingVolumeOpenGL2
        #vtkTclTk
        vtkTestingCore
        #vtkTestingGenericBridge
        #vtkTestingIOSQL
        vtkTestingRendering
        #vtkUtilitiesBenchmarks
        #vtkUtilitiesEncodeString
        #vtkUtilitiesHashSource
        #vtkVPIC
        vtkViewsContext2D
        vtkViewsCore
        #vtkViewsGeovis
        vtkViewsInfovis
        #vtkViewsQt
        #vtkWrappingJava
        #vtkWrappingPythonCore
        #vtkWrappingTcl
        vtkWrappingTools
        #vtkalglib
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
        #vtklibproj4
        vtklibxml2
        vtklz4
        #vtkm
        #vtkmpi4py
        vtknetcdf
        #vtknetcdfcpp
        #vtkoggtheora
        vtkpng
        vtksqlite
        vtksys
        vtktiff
        vtkverdict
        #vtkxdmf2
        #vtkxdmf3
        #vtkzfp
        vtkzlib
        #wslink

        vtkFiltersSources
        vtkIOLegacy
        vtkFiltersExtraction
        vtkFiltersGeometry
        vtkFiltersModeling
        vtkRenderingFreeType

      vtkFiltersSources
      vtkIOLegacy
      vtkFiltersExtraction
      vtkFiltersGeometry
      vtkFiltersModeling
      vtkRenderingFreeType
      vtkCommonColor
      vtkCommonCore
      vtkFiltersSources
      vtkIOImage
      vtkInteractionStyle
      vtkRenderingContextOpenGL2
      vtkRenderingCore
      vtkRenderingFreeType
      vtkRenderingGL2PSOpenGL2
      vtkRenderingOpenGL2 QUIET)
    if (NOT VTK_FOUND)
      message("Skipping ImageWriter: ${VTK_NOT_FOUND_MESSAGE}")
      return ()
    endif()
    message (STATUS "VTK_VERSION: ${VTK_VERSION}")
    if (VTK_VERSION VERSION_LESS "8.90.0")
      # old system
      include(${VTK_USE_FILE})
  endif()
  #include(${VTK_USE_FILE})
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

FILE(GLOB SOURCES ${GInt_SOURCE_DIR}/intensityFieldContourView/*.cxx)

ADD_EXECUTABLE(intensityFieldContourView ${HEADERS} ${SOURCES})
  
TARGET_LINK_LIBRARIES(intensityFieldContourView LINK_PUBLIC GInt ${VTK_LIBRARIES})


message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (NOT VTK_VERSION VERSION_LESS "8.90.0")
    message("auto init vtk modules")
    vtk_module_autoinit(TARGETS intensityFieldContourView GInt MODULES ${VTK_LIBRARIES})
endif ()

if(APPLE)
  set_target_properties(intensityFieldContourView PROPERTIES XCODE_ATTRIBUTE_CC "${APPLE_CLANG_COMPILER}")
  set_target_properties(intensityFieldContourView PROPERTIES XCODE_ATTRIBUTE_CXX "${APPLE_CLANG_COMPILER}")
endif()
                             
 
