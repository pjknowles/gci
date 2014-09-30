cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
subdirs(FCIdump)

file (GLOB SRCS *.cpp *.h)
add_executable(gci.exe ${SRCS})

include_directories(FCIdump)
link_directories(${CMAKE_BINARY_DIR}/FCIdump)
target_link_libraries(gci.exe FCIdump)

get_filename_component(COMPILER_NAME  ${CMAKE_CXX_COMPILER} NAME)
if (${COMPILER_NAME} STREQUAL "mpic++")
find_package(ppidd)
if (GA_FOUND)
include_directories(${PPIDD_INCLUDE_DIRS})
target_link_libraries(gci.exe ${PPIDD_LIBRARIES} ${GA_LIBRARIES})
add_definitions(-DGCI_PARALLEL)
endif (GA_FOUND)
endif (${COMPILER_NAME} STREQUAL "mpic++")

find_package(MKL)
if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
target_link_libraries(gci.exe ${MKL_LIBRARY})
else (MKL_FOUND)
target_link_libraries(gci.exe -llapack)
endif (MKL_FOUND)

configure_file(gci.fcidump gci.fcidump COPYONLY)
