cmake_minimum_required(VERSION 3.19)
project(polyMpmTest VERSION 0.0.1 LANGUAGES CXX Fortran)

include(CTest)
enable_testing()

function(export_target tgt_name)
  install(TARGETS ${tgt_name} EXPORT ${tgt_name}-target
      RUNTIME DESTINATION bin
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib)
  install(EXPORT ${tgt_name}-target NAMESPACE ${PROJECT_NAME}::
          DESTINATION lib/cmake/${PROJECT_NAME})
  set(${PROJECT_NAME}_EXPORTED_TARGETS
      ${${PROJECT_NAME}_EXPORTED_TARGETS} ${tgt_name} PARENT_SCOPE)
endfunction(export_target)

macro(polyMpmTest_export_lib target headers)
  export_target(${target})
  install(FILES ${headers} DESTINATION include)
endmacro(polyMpmTest_export_lib)

#pumipic dependency
set(CMAKE_CXX_EXTENSIONS Off) #kokkos does not allow extensions
message(STATUS "Found PUMIpic: ${pumipic_VERSION}")
find_package(pumipic REQUIRED)
#Clear the omegah compilation flags that it passes to cuda.  Using the
# kokkos target, and nvcc_wrapper, provide sufficient flags.
set_property(TARGET Omega_h::omega_h PROPERTY INTERFACE_COMPILE_OPTIONS "")

#kokkos - pumipic should provide this
find_package(Kokkos REQUIRED)

option(POLYMPO_USE_NETCDF "Enable the use of NetCDF" ON)
if(POLYMPO_USE_NETCDF)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(NetCDF REQUIRED IMPORTED_TARGET netcdf>=4.7.3)
  pkg_check_modules(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran>=4.5.2)
endif()

# Set a default build type if none was specified
# from: https://www.kitware.com/cmake-and-the-default-build-type/
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#install fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

add_subdirectory (src)
add_subdirectory (test)
