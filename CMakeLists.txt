cmake_minimum_required(VERSION 3.19)
project(polyMpmTest VERSION 0.0.1 LANGUAGES CXX)

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

#kokkos does not allow extensions
set(CMAKE_CXX_EXTENSIONS Off)
find_package(Kokkos REQUIRED)

option(POLYMPO_USE_NETCDF "Enable the use of NetCDF" OFF)
if(POLYMPO_USE_NETCDF)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(NetCDF REQUIRED IMPORTED_TARGET netcdf>=4.7.3)
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

add_subdirectory (src)
add_subdirectory (test)
