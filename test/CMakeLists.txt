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

add_test(NAME unit_test COMMAND ./unitTest)
add_test(NAME mp_unit_test COMMAND ./mpUnitTest)
add_test(NAME test_wachspress COMMAND ./testWachspress)
add_test(NAME test_timing COMMAND ./timeAssmblyWachspress 5)
add_test(NAME test_reconstruction COMMAND ./testReconstruction)

pmt_add_exe(testFortran testFortran.f90)
target_link_libraries(testFortran pmtUtils polympo_fortran polyMpmTest MPI::MPI_Fortran)
set_property(TARGET testFortran PROPERTY LINKER_LANGUAGE Fortran)
set_source_files_properties(testFortran.f90 PROPERTIES COMPILE_FLAGS -cpp)
install(TARGETS testFortran DESTINATION bin)

add_test(NAME testFortran COMMAND ./testFortran)

