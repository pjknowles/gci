cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
subdirs(FCIdump symmetry_matrix IterativeSolver)
enable_language(Fortran)

file (GLOB SRCS src/*.cpp src/*.h src/*.F90 extra/*.cpp extra/*.h)
add_definitions(-DNOMAIN)
add_executable(gci.exe ${SRCS})

include_directories(src extra FCIdump IterativeSolver IterativeSolver/Eigen symmetry_matrix symmetry_matrix/memory)
target_link_libraries(gci.exe FCIdump memory IterativeSolver symmetry_matrix)
link_directories(${CMAKE_BINARY_DIR}/FCIdump ${CMAKE_BINARY_DIR}/symmetry_matrix ${CMAKE_BINARY_DIR}/IterativeSolver ${CMAKE_BINARY_DIR}/symmetry_matrix/memory)

get_filename_component(COMPILER_NAME  ${CMAKE_CXX_COMPILER} NAME)
if (${COMPILER_NAME} STREQUAL "mpic++" OR ${COMPILER_NAME} STREQUAL "mpicxx")
find_package(PPIDD)
if (GA_FOUND)
    message("Building in parallel, PPIDD includes at ${PPIDD_INCLUDE_DIRS}")
include_directories(${PPIDD_INCLUDE_DIRS})
target_link_libraries(gci.exe ${PPIDD_LIBRARIES} ${GA_LIBRARIES})
add_definitions(-DGCI_PARALLEL)
else(GA_FOUND)
    message("Building in serial")
endif (GA_FOUND)
endif (${COMPILER_NAME} STREQUAL "mpic++" OR ${COMPILER_NAME} STREQUAL "mpicxx")

find_package(MKL)
if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
message("MKL link libraries are ${MKL_LIBRARIES}")
target_link_libraries(gci.exe ${MKL_LIBRARIES})
else (MKL_FOUND)
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(gci.exe ${BLAS_LIBRARY} memory -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl)
message("Blas link library is ${BLAS_LIBRARIES}")
else (BLAS_FOUND)
target_link_libraries(gci.exe -llapack memory -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl)
endif (BLAS_FOUND)
endif (MKL_FOUND)

configure_file(gci.fcidump gci.fcidump COPYONLY)

find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(${PROJECT_NAME}-doc ALL
    DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    )
file(GLOB headers *.h *.fh)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
DEPENDS ${SRC_LIST} ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
