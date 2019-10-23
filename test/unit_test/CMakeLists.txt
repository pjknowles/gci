get_dependency(googletest)
include(GoogleTest)

#This is outdated
#add_executable(test-gci
#        testVibSpace.cpp
#        testHProduct.cpp
#        testMixedOperator.cpp
#        testHProductSet.cpp
#        )
#target_link_libraries(test-gci PUBLIC gci gmock_main)
#gtest_discover_tests(test-gci)

add_executable(test-gci-parallel
        main.cpp
        testArray.cpp
        testDavidson.cpp
        )
target_link_libraries(test-gci-parallel PUBLIC gci gmock_main)
add_test(NAME test-gci-parallel
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        ${CMAKE_CURRENT_BINARY_DIR}/test-gci-parallel ${MPIEXEC_POSTFLAGS}
        )

