# msceer_dmt: Python bindings for the discrete Morse theory kernels.
#
# Header-only C++ under include/ (dmt_core + gi_dmt_*.h). Links GInt for the
# static offset tables of the topological grid; deliberately does not link
# msc_2d_lib, whose API is consumed elsewhere and stays untouched.

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx")
file(GLOB HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

pybind11_add_module(msceer_dmt MODULE ${SOURCES} ${HEADERS})

set_target_properties(msceer_dmt PROPERTIES
    OUTPUT_NAME "msceer_dmt"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

target_compile_features(msceer_dmt PRIVATE cxx_std_17)
if(MSVC)
    target_compile_options(msceer_dmt PRIVATE /bigobj)
endif()

target_link_libraries(msceer_dmt
    PRIVATE GInt
    PRIVATE pybind11::module
)
if(OpenMP_FOUND)
    target_link_libraries(msceer_dmt PRIVATE OpenMP::OpenMP_CXX)
endif()

# device='gpu' on Gradient3D.pls / pls_br / dmf when the CUDA path is built.
if(TARGET gpu_dgrad)
    target_link_libraries(msceer_dmt PRIVATE gpu_dgrad)
    target_compile_definitions(msceer_dmt PRIVATE DMT_PY_GPU)
endif()
