set(HEADERS
  mesh.hpp
  materialPoints.hpp
  wachpressBasis.hpp
  MPM.hpp
  utils.hpp
  assembly.hpp
)

set(SOURCES
  mesh.cpp
)

add_library (polyMpmTest-core ${SOURCES})

target_include_directories(polyMpmTest-core INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:include>)
target_link_libraries(polyMpmTest-core PUBLIC Kokkos::kokkos)
target_link_libraries(polyMpmTest-core PRIVATE PkgConfig::NetCDF)

polyMpmTest_export_lib(polyMpmTest-core "${HEADERS}")

add_library (polyMpmTest INTERFACE)
target_link_libraries(polyMpmTest INTERFACE ${polyMpmTest_EXPORTED_TARGETS})
export_target(polyMpmTest)

