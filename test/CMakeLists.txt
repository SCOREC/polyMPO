option(POLYMPO_USE_NETCDF "Enable the use of NetCDF" ON)
if(POLYMPO_USE_NETCDF)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(NetCDF REQUIRED IMPORTED_TARGET netcdf>=4.7.3 netcdf_parallel>=4.7.3)
  pkg_search_module(NetCDF_Fortran REQUIRED IMPORTED_TARGET netcdf-fortran>=4.5.2 netcdf-fortran_parallel>=4.5.2)
endif()


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

#add exe
function(pmpo_add_exe exename srcname)
  if(IS_TESTING)
    add_executable(${exename} ${srcname})
  else()
    add_executable(${exename} EXCLUDE_FROM_ALL ${srcname})
  endif()
  target_link_libraries (${exename} polyMPO pmpoUtils)
  if(POLYMPO_USE_NETCDF)
    target_link_libraries (${exename} PkgConfig::NetCDF)
    target_compile_definitions(${exename} PRIVATE POLYMPO_HAS_NETCDF)
    target_include_directories(${exename} PRIVATE ${NetCDF_INCLUDEDIR})
  endif()
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD "17")
  set_property(TARGET ${exename} PROPERTY CXX_STANDARD_REQUIRED ON)
  set_property(TARGET ${exename} PROPERTY CXX_EXTENSIONS OFF)
  install(TARGETS ${exename} DESTINATION bin)
endfunction()

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

pmpo_add_exe(testWachspress testWachspress.cpp)
pmpo_add_exe(unitTest unitTest.cpp)
pmpo_add_exe(mpUnitTest mpUnitTest.cpp)
pmpo_add_exe(testReconstruction testReconstruction.cpp)
pmpo_add_exe(testTracking testTracking.cpp)
pmpo_add_exe(testRebuild testRebuild.cpp)
pmpo_add_exe(timeAssmblyWachspress testTiming.cpp)
pmpo_add_fortran_exe(testFortranInit testFortranInit.f90)
pmpo_add_fortran_exe(testFortran testFortran.f90)
pmpo_add_fortran_exe(testFortranReadMPAS testFortranReadMPAS.f90)
pmpo_add_fortran_exe(testFortranMPMeshModule testFortranMPMeshModule.f90)
pmpo_add_fortran_exe(testFortranMPAppIDs testFortranMPAppIDs.f90)
pmpo_add_fortran_exe(testFortranMPAdvection testFortranMPAdvection.f90)
pmpo_add_fortran_exe(testFortranCreateRebuildMPs testFortranCreateRebuildMPs.f90)
pmpo_add_fortran_exe(testFortranMPReconstruction testFortranMPReconstruction.f90)

#add test
find_program(MPIRUN_EXECUTABLE NAMES mpirun) 

function(pmpo_add_test testname command)
  if(MPIRUN_EXECUTABLE)
    add_test(NAME ${testname} 
             COMMAND mpirun --bind-to core -np 1 ${command} ${ARGN})
  else()
    add_test(NAME ${testname} COMMAND ${command} ${ARGN})
  endif()
endfunction()

pmpo_add_test(unit_test ./unitTest)
pmpo_add_test(mp_unit_test ./mpUnitTest)
pmpo_add_test(test_reconstruction ./testReconstruction)
pmpo_add_test(test_tracking ./testTracking)
pmpo_add_test(test_rebuild ./testRebuild)
pmpo_add_test(testFortranInit ./testFortranInit)
pmpo_add_test(testFortran ./testFortran)
pmpo_add_test(testFortranMPMeshModule ./testFortranMPMeshModule)
pmpo_add_test(testFortranMPAppIDs ./testFortranMPAppIDs)

#set NC file for test
set(TEST_NC_FILE_PLANAR "${CMAKE_SOURCE_DIR}/test/sample_mpas_meshes/planar_nonuniform_cvt_for_square_673elms.nc")
set(TEST_NC_FILE_SPHERICAL
"${CMAKE_SOURCE_DIR}/test/sample_mpas_meshes/spherical_cvt_642elms.nc")
if(EXISTS ${TEST_NC_FILE_SPHERICAL})
    message(STATUS "Found Spherical NC File: ${TEST_NC_FILE_SPHERICAL}") 
    set(TEST_NC_FILE ${TEST_NC_FILE_SPHERICAL})
    # this test only works for spherical NC files!
    pmpo_add_test(testFortranMPAdvection ./testFortranMPAdvection 5 5 ${TEST_NC_FILE}) 
else()
    set(TEST_NC_FILE ${TEST_NC_FILE_PLANAR})
endif()
pmpo_add_test(testFortranReadMPAS ./testFortranReadMPAS ${TEST_NC_FILE})
pmpo_add_test(testFortranCreateRebuildMPs ./testFortranCreateRebuildMPs ${TEST_NC_FILE})
pmpo_add_test(testFortranMPReconstruction ./testFortranMPReconstruction ${TEST_NC_FILE})

pmpo_add_test(test_timing ./timeAssmblyWachspress 1 ${TEST_NC_FILE})
pmpo_add_test(test_wachspress ./testWachspress ${TEST_NC_FILE})

bob_end_subdir()
