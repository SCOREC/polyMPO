add_library(readMPAS readMPAS.f90)
target_link_libraries(readMPAS pmpoUtils polympo_fortran polyMPO MPI::MPI_Fortran)
target_link_libraries(readMPAS PkgConfig::NetCDF_Fortran)
target_include_directories(readMPAS PRIVATE ${NetCDF_Fortran_INCLUDEDIR})

add_library(pmpoUtils testUtils.cpp)
target_link_libraries(pmpoUtils PRIVATE polyMPO)
target_link_libraries(pmpoUtils PRIVATE readMPAS)
set_property(TARGET pmpoUtils PROPERTY CXX_STANDARD "17")
set_property(TARGET pmpoUtils PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET pmpoUtils PROPERTY CXX_EXTENSIONS OFF)
 
target_include_directories(pmpoUtils PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

function(pmpo_add_exe exename srcname)
  if(BUILD_TESTING)
    add_executable(${exename} ${srcname})
  else()
    add_executable(${exename} EXCLUDE_FROM_ALL ${srcname})
  endif()
  target_link_libraries (${exename} polyMPO pmpoUtils)
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD "17")
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD_REQUIRED ON)
  set_property(TARGET ${exename} PROPERTY CXX_EXTENSIONS OFF)
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

pmpo_add_exe(testWachspress testWachspress.cpp)
pmpo_add_exe(unitTest unitTest.cpp)
pmpo_add_exe(mpUnitTest mpUnitTest.cpp)
pmpo_add_exe(testReconstruction testReconstruction.cpp)
pmpo_add_exe(testTracking testTracking.cpp)
pmpo_add_exe(testRebuild testRebuild.cpp)
pmpo_add_exe(timeAssmblyWachspress testTiming.cpp)

add_test(NAME unit_test COMMAND mpirun --bind-to core -np 1 ./unitTest)
add_test(NAME mp_unit_test COMMAND mpirun --bind-to core -np 1 ./mpUnitTest)
add_test(NAME test_wachspress COMMAND mpirun --bind-to core -np 1 ./testWachspress)
add_test(NAME test_reconstruction COMMAND mpirun --bind-to core -np 1 ./testReconstruction)
add_test(NAME test_tracking COMMAND mpirun --bind-to core -np 1 ./testTracking)
add_test(NAME test_rebuild COMMAND mpirun --bind-to core -np 1 ./testRebuild)

function(pmpo_add_fortran_exe exename srcname)
  pmpo_add_exe(${exename} ${srcname})
  target_link_libraries(${exename} pmpoUtils polympo_fortran polyMPO MPI::MPI_Fortran)
  if(${NetCDF_Fortran_FOUND})
    target_link_libraries(${exename} PkgConfig::NetCDF_Fortran)
    target_include_directories(${exename} PRIVATE ${NetCDF_Fortran_INCLUDEDIR})
  endif()
  set_property(TARGET ${exename} PROPERTY LINKER_LANGUAGE Fortran)
  set_source_files_properties(${scrname} PROPERTIES COMPILE_FLAGS -cpp)
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

pmpo_add_fortran_exe(testFortranInit testFortranInit.f90)
add_test(NAME testFortranInit COMMAND mpirun --bind-to core -np 1 ./testFortranInit)

pmpo_add_fortran_exe(testFortran testFortran.f90)
add_test(NAME testFortran COMMAND mpirun --bind-to core -np 1 ./testFortran)

set(TEST_NC_FILE_PLANAR "${CMAKE_SOURCE_DIR}/test/sample_mpas_meshes/planar_nonuniform_cvt_for_square_737elms.nc")
set(TEST_NC_FILE_SPHERICAL "${CMAKE_SOURCE_DIR}/test/sample_mpas_meshes/ocean.nc")
if(EXISTS ${TEST_NC_FILE_SPHERICAL})
    set(TEST_NC_FILE ${TEST_NC_FILE_SPHERICAL})
else()
    set(TEST_NC_FILE ${TEST_NC_FILE_PLANAR})
endif()
pmpo_add_fortran_exe(testFortranReadMPAS testFortranReadMPAS.f90)
add_test(NAME testFortranReadMPAS 
         COMMAND mpirun --bind-to core -np 1 ./testFortranReadMPAS ${TEST_NC_FILE})

add_test(NAME test_timing
         COMMAND mpirun --bind-to core -np 1 ./timeAssmblyWachspress 1 ${TEST_NC_FILE})
