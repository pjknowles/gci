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

add_executable(test-gci-array
        main.cpp
        testArray.cpp
        )
target_link_libraries(test-gci-array PUBLIC gci gmock_main)
add_test(NAME test-gci-array
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        ${CMAKE_CURRENT_BINARY_DIR}/test-gci-array ${MPIEXEC_POSTFLAGS}
        )

