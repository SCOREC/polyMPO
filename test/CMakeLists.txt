add_library(pmtUtils testUtils.cpp)
target_link_libraries(pmtUtils PRIVATE polyMpmTest)
set_property(TARGET pmtUtils PROPERTY CXX_STANDARD "17")
set_property(TARGET pmtUtils PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET pmtUtils PROPERTY CXX_EXTENSIONS OFF)
 
target_include_directories(pmtUtils PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

function(pmt_add_exe exename srcname)
  if(BUILD_TESTING)
    add_executable(${exename} ${srcname})
  else()
    add_executable(${exename} EXCLUDE_FROM_ALL ${srcname})
  endif()
  target_link_libraries (${exename} polyMpmTest pmtUtils)
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD "17")
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD_REQUIRED ON)
  set_property(TARGET ${exename} PROPERTY CXX_EXTENSIONS OFF)
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

pmt_add_exe(testWachspress testWachspress.cpp)
pmt_add_exe(unitTest unitTest.cpp)
pmt_add_exe(mpUnitTest mpUnitTest.cpp)
pmt_add_exe(timeAssmblyWachspress testTiming.cpp)
pmt_add_exe(testReconstruction testReconstruction.cpp)
pmt_add_exe(testTracking testTracking.cpp)
pmt_add_exe(testRebuild testRebuild.cpp)

add_test(NAME unit_test COMMAND ./unitTest)
add_test(NAME mp_unit_test COMMAND ./mpUnitTest)
add_test(NAME test_wachspress COMMAND ./testWachspress)
add_test(NAME test_timing COMMAND ./timeAssmblyWachspress 5)
add_test(NAME test_reconstruction COMMAND ./testReconstruction)
add_test(NAME test_tracking COMMAND ./testTracking)
add_test(NAME test_rebuild COMMAND ./testRebuild)

function(pmt_add_fortran_exe exename srcname)
  pmt_add_exe(${exename} ${srcname})
  target_link_libraries(${exename} pmtUtils polympo_fortran polyMpmTest MPI::MPI_Fortran)
  if(${NetCDF_Fortran_FOUND})
    target_link_libraries(${exename} PkgConfig::NetCDF_Fortran)
    target_include_directories(${exename} PRIVATE ${NetCDF_Fortran_INCLUDEDIR})
  endif()
  set_property(TARGET ${exename} PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(${scrname} PROPERTIES COMPILE_FLAGS -cpp)
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

pmt_add_fortran_exe(testFortranInit testFortranInit.f90)
add_test(NAME testFortranInit COMMAND ./testFortranInit)

pmt_add_fortran_exe(testFortran testFortran.f90)
add_test(NAME testFortran COMMAND ./testFortran)

set(TEST_NC_FILE "${CMAKE_SOURCE_DIR}/test/sample_mpas_meshes/planar_nonuniform_cvt_for_square_737elms.nc")
pmt_add_fortran_exe(testFortranReadMPAS testFortranReadMPAS.f90)
add_test(NAME testFortranReadMPAS 
         COMMAND ./testFortranReadMPAS ${TEST_NC_FILE})

