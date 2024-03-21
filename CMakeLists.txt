cmake_minimum_required(VERSION 3.19)
project(polyMPO VERSION 0.0.2 LANGUAGES CXX Fortran)

include(CTest)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
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

macro(polyMPO_export_lib target headers)
  export_target(${target})
  install(FILES ${headers} DESTINATION include)
endmacro(polyMPO_export_lib)

#pumipic dependency
set(CMAKE_CXX_EXTENSIONS Off) #kokkos does not allow extensions
message(STATUS "Found PUMIpic: ${pumipic_VERSION}")
find_package(pumipic REQUIRED)
#Clear the omegah compilation flags that it passes to cuda.  Using the
# kokkos target, and nvcc_wrapper, provide sufficient flags.
set_property(TARGET Omega_h::omega_h PROPERTY INTERFACE_COMPILE_OPTIONS "")

#kokkos - pumipic should provide this
find_package(Kokkos REQUIRED)

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

option(IS_TESTING "Build for CTest" OFF)
if(IS_TESTING)
  add_subdirectory (test)
endif()

## export the library
set_target_properties(polyMPO PROPERTIES
  PUBLIC_HEADER "${POLYMPO_HEADERS}")
target_include_directories(polyMPO
  INTERFACE
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
install(
  TARGETS polyMPO polyMPO-core
  EXPORT polyMPO-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/polympo-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyMPO
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/polympo-config-version.cmake"
  COMPATIBILITY AnyNewerVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/polympo-config.cmake"
  "${PROJECT_BINARY_DIR}/polympo-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyMPO)

install(
  EXPORT polyMPO-targets
  NAMESPACE polyMPO::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/polyMPO)