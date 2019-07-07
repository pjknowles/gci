FetchContent_GetProperties(googletest_dep)
if (NOT googletest_dep_POPULATED)
    FetchContent_Populate(googletest_dep)
    add_subdirectory(${googletest_dep_SOURCE_DIR} ${googletest_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()

#add_executable(test_gci)
#target_link_libraries(test_gci PUBLIC gci-lib gmock_main)
