########## Dependencies ##########


########## Target ##########



########## Main executable ##########

add_executable(gci.exe gci.cpp)
target_link_libraries(gci.exe PUBLIC gci)
set_target_properties(gci.exe PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
