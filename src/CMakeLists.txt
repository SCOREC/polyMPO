set(HEADERS
  pmpo_mesh.hpp
  pmpo_materialPoints.hpp
  pmpo_wachspressBasis.hpp
  pmpo_MPMesh.hpp
  pmpo_utils.hpp
  pmpo_assembly.hpp
  pmpo_c.h
  pmpo_createTestMPMesh.hpp
)

set(SOURCES
  pmpo_materialPoints.cpp
  pmpo_mesh.cpp
  pmpo_utils.cpp
  pmpo_MPMesh.cpp
  pmpo_c.cpp
  pmpo_createTestMPMesh.cpp
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
add_library(polympo_fortran pmpo_fortran.f90)
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

