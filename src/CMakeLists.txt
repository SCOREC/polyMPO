set(HEADERS
  pmpo_mesh.hpp
  pmpo_materialPoints.hpp
  pmpo_wachspressBasis.hpp
  pmpo_MPMesh.hpp
  pmpo_utils.hpp
  pmpo_MPMesh_assembly.hpp
  pmpo_c.h
  pmpo_createTestMPMesh.hpp
  pmpo_const_relation.hpp
)

set(SOURCES
  pmpo_materialPoints.cpp
  pmpo_mesh.cpp
  pmpo_utils.cpp
  pmpo_MPMesh.cpp
  pmpo_c.cpp
  pmpo_createTestMPMesh.cpp
  ../test/testMPAppIDs.cpp
)

add_library(polyMPO-core ${SOURCES})

# Default ON only for Cray GPU builds (Kokkos CUDA + Cray MPICH GTL available).
# Can still be overridden with -DpolyMPO_ENABLE_GPU_AWARE_MPI=ON/OFF.
set(polyMPO_GPU_AWARE_MPI_DEFAULT OFF)
if(Kokkos_ENABLE_CUDA
   AND DEFINED ENV{PE_MPICH_GTL_DIR_nvidia80}
   AND DEFINED ENV{PE_MPICH_GTL_LIBS_nvidia80})
  set(polyMPO_GPU_AWARE_MPI_DEFAULT ON)
endif()

option(polyMPO_ENABLE_GPU_AWARE_MPI
  "Use GPU-aware (CUDA-aware) MPI for halo exchange"
  ${polyMPO_GPU_AWARE_MPI_DEFAULT})
message(STATUS "polyMPO_ENABLE_GPU_AWARE_MPI: ${polyMPO_ENABLE_GPU_AWARE_MPI}")

if(polyMPO_ENABLE_GPU_AWARE_MPI)
  if(NOT Kokkos_ENABLE_CUDA)
    message(FATAL_ERROR
      "polyMPO_ENABLE_GPU_AWARE_MPI=ON requires Kokkos built with CUDA.")
  endif()

  if(DEFINED ENV{PE_MPICH_GTL_DIR_nvidia80}
     AND DEFINED ENV{PE_MPICH_GTL_LIBS_nvidia80})
    target_link_options(polyMPO-core PUBLIC
      "$ENV{PE_MPICH_GTL_DIR_nvidia80}"
      "$ENV{PE_MPICH_GTL_LIBS_nvidia80}"
    )
    message(STATUS "polyMPO: enabling Cray MPICH CUDA GTL linkage")
  else()
    message(FATAL_ERROR
      "Cray MPICH CUDA GTL environment variables are unavailable. "
      "Load cudatoolkit and craype-accel-nvidia80, or set "
      "CRAY_ACCEL_TARGET=nvidia80 before configuring.")
  endif()

  target_compile_definitions(polyMPO-core PUBLIC GPU_AWARE_MPI)
endif()

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

polyMPO_export_lib(polyMPO-core "${HEADERS}")

### Fortran inteface ###
add_library(polympo_fortran pmpo_fortran.f90)
target_link_libraries(polympo_fortran PUBLIC MPI::MPI_Fortran)
target_include_directories(polympo_fortran INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)
set(HEADERS)
polyMPO_export_lib(polympo_fortran "${HEADERS}")

add_library(polyMPO INTERFACE)
target_link_libraries(polyMPO INTERFACE ${polyMPO_EXPORTED_TARGETS})
bob_export_target(polyMPO)

bob_end_subdir()
