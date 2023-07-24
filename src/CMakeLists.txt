set(HEADERS
  pmt_mesh.hpp
  pmt_materialPoints.hpp
  pmt_wachspressBasis.hpp
  pmt_MPMesh.hpp
  pmt_utils.hpp
  pmt_assembly.hpp
  pmo_c.h
  pmo_createTestMPMesh.hpp
)

set(SOURCES
  pmt_materialPoints.cpp
  pmt_mesh.cpp
  pmt_utils.cpp
  pmt_MPMesh.cpp
  pmo_c.cpp
  pmo_createTestMPMesh.cpp
)

add_library(polyMPO-core ${SOURCES})
set_property(TARGET polyMPO-core PROPERTY CXX_STANDARD "17")
set_property(TARGET polyMPO-core PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET polyMPO-core PROPERTY CXX_EXTENSIONS OFF)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(polyMPO-core PUBLIC ${MATH_LIBRARY})
endif()

target_include_directories(polyMPO-core
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/test>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(polyMPO-core PUBLIC pumipic::pumipic)
if(${NetCDF_FOUND})
  target_compile_definitions(polyMPO-core PRIVATE POLYMPO_HAS_NETCDF)
  target_link_libraries(polyMPO-core PRIVATE PkgConfig::NetCDF)
endif()

polyMPO_export_lib(polyMPO-core "${HEADERS}")

add_library(polyMPO INTERFACE)
target_link_libraries(polyMPO INTERFACE ${polyMPO_EXPORTED_TARGETS})
export_target(polyMPO)

### Fortran inteface ###
add_library(polympo_fortran pmo_fortran.f90)
target_link_libraries(polympo_fortran PUBLIC MPI::MPI_Fortran)
if(${NetCDF_Fortran_FOUND})
  target_link_libraries(polympo_fortran PRIVATE PkgConfig::NetCDF_Fortran)
  # why doesn't target_link_libraries do the following?
  target_include_directories(polympo_fortran PRIVATE ${NetCDF_Fortran_INCLUDEDIR})
endif()
target_include_directories(polympo_fortran INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
set(HEADERS)
polyMPO_export_lib(polympo_fortran "${HEADERS}")

