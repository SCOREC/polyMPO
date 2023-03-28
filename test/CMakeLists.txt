add_library(pmtUtils testUtils.cpp)
target_link_libraries(pmtUtils PRIVATE polyMpmTest)
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
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

pmt_add_exe(testWachspress testWachspress.cpp)
pmt_add_exe(unitTest unitTest.cpp)
pmt_add_exe(timeAssmblyWachspress testTiming.cpp)
pmt_add_exe(testT2L testT2L.cpp)

add_test(NAME unit_test COMMAND ./unitTest)
add_test(NAME test_wachspress COMMAND ./testWachspress)
add_test(NAME test_T2L COMMAND ./testT2L)
add_test(NAME test_timing COMMAND ./timeAssmblyWachspress 5)
