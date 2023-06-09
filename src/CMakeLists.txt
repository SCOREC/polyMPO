set(HEADERS
  pmt_mesh.hpp
  pmt_materialPoints.hpp
  pmt_wachspressBasis.hpp
  pmt_MPMesh.hpp
  pmt_utils.hpp
  pmt_assembly.hpp
)

set(SOURCES
  pmt_materialPoints.cpp
  pmt_mesh.cpp
  pmt_utils.cpp
)

add_library(polyMpmTest-core ${SOURCES})
set_property(TARGET polyMpmTest-core PROPERTY CXX_STANDARD "17")
set_property(TARGET polyMpmTest-core PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET polyMpmTest-core PROPERTY CXX_EXTENSIONS OFF)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(polyMpmTest-core PUBLIC ${MATH_LIBRARY})
endif()

target_include_directories(polyMpmTest-core INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(polyMpmTest-core PUBLIC pumipic::pumipic)
if(${NetCDF_FOUND})
  target_compile_definitions(polyMpmTest-core PRIVATE POLYMPO_HAS_NETCDF)
  target_link_libraries(polyMpmTest-core PRIVATE PkgConfig::NetCDF)
endif()

polyMpmTest_export_lib(polyMpmTest-core "${HEADERS}")

add_library(polyMpmTest INTERFACE)
target_link_libraries(polyMpmTest INTERFACE ${polyMpmTest_EXPORTED_TARGETS})
export_target(polyMpmTest)

