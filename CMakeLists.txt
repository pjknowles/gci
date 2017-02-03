cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
enable_language(Fortran)
set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fdefault-integer-8")
set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fdefault-integer-8")

file (GLOB SRCS src/*.cpp src/*.h src/*.F90 extra/*.cpp extra/*.h
    symmetry_matrix/memory/memory.F90 symmetry_matrix/memory/bytestreamC.cpp symmetry_matrix/memory/bytestreamF.F90
    symmetry_matrix/SMat.cpp symmetry_matrix/SMatfunction.cpp symmetry_matrix/Operator.cpp symmetry_matrix/SMatMat.cpp
    FCIdump/FCIdump.cpp
    IterativeSolver/IterativeSolver/ISDavidson.cpp
    IterativeSolver/IterativeSolver/ISDiis.cpp
    IterativeSolver/IterativeSolver/ISRSPT.cpp
    IterativeSolver/IterativeSolver/IterativeSolver.cpp
    IterativeSolver/IterativeSolver/PagedParameterVector.cpp
    IterativeSolver/IterativeSolver/Storage.cpp
    IterativeSolver/IterativeSolver/SimpleParameterVector.cpp
    )
add_definitions(-DNOMAIN)
add_executable(gci.exe ${SRCS})

include_directories(src extra FCIdump IterativeSolver/IterativeSolver IterativeSolver IterativeSolver/Eigen symmetry_matrix symmetry_matrix/memory)

get_filename_component(COMPILER_NAME  ${CMAKE_CXX_COMPILER} NAME)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(gci.exe ${MPI_CXX_LIBRARIES})

find_package(MKL)
if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
message("MKL link libraries are ${MKL_LIBRARIES}")
target_link_libraries(gci.exe ${MKL_LIBRARIES})
else (MKL_FOUND)
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
include_directories(${BLAS_INCLUDE_DIRS})
target_link_libraries(gci.exe ${BLAS_LIBRARY} -lm -ldl)
message("Blas link library is ${BLAS_LIBRARIES}")
else (BLAS_FOUND)
target_link_libraries(gci.exe -llapack -lm -ldl)
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
