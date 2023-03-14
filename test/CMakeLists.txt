add_executable (testFoo testFoo.cpp testUtils.cpp)
target_link_libraries (testFoo polyMpmTest m)
install(TARGETS testFoo DESTINATION bin)

