cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)

file (GLOB SRCS *.cpp *.h)
add_executable(gci.exe ${SRCS})

get_filename_component(COMPILER_NAME  ${CMAKE_CXX_COMPILER} NAME)
if (${COMPILER_NAME} STREQUAL "mpic++")
find_package(ppidd)
if (GA_FOUND)
include_directories(${PPIDD_INCLUDE_DIRS})
target_link_libraries(gci.exe ${PPIDD_LIBRARIES} ${GA_LIBRARIES})
add_definitions(-DGCI_PARALLEL)
endif (GA_FOUND)
endif (${COMPILER_NAME} STREQUAL "mpic++")


configure_file(FCIDUMP FCIDUMP COPYONLY)
