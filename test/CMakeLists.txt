add_executable (testFoo testFoo.cpp testUtils.cpp)
add_executable (unitTest unitTest.cpp testUtils.cpp)

target_link_libraries (testFoo polyMpmTest m)
target_link_libraries (unitTest polyMpmTest m)
install(TARGETS testFoo DESTINATION bin)

enable_testing()

add_test(
    unit_test unitTest
    #COMMAND $TARGET_FILE:unitTest
    )
