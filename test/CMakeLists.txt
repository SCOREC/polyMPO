add_executable (testFoo testFoo.cpp testUtils.cpp)
add_executable (unitTest unitTest.cpp testUtils.cpp)
add_executable (timeAssmblyWachspress testTiming.cpp testUtils.cpp)

target_link_libraries (testFoo polyMpmTest m)
target_link_libraries (unitTest polyMpmTest m)
target_link_libraries (timeAssmblyWachspress polyMpmTest m)

install(TARGETS testFoo DESTINATION bin)
install(TARGETS timeAssmblyWachspress DESTINATION bin)

#enable_testing()

add_test(
    unit_test unitTest
    )
add_test(
    test_watchpress testFoo
    )
