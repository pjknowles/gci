cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
enable_language(Fortran)
set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-fdefault-integer-8" _gnu8)
if(_gnu8)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
endif()
check_fortran_compiler_flag("-i8" _i8)
if(_i8)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
endif()
check_fortran_compiler_flag("-fast" _fast)
if(_fast)
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fast")
endif()

file (GLOB SRCS src/*.cpp src/*.h src/*.F90 extra/*.cpp extra/*.h
    submodules/Profiler/Profiler.h
    submodules/Profiler/Profiler.cpp
    submodules/memory/memory.F90 submodules/memory/bytestreamC.cpp submodules/memory/bytestreamF.F90
    submodules/symmetry_matrix/SMat.cpp submodules/symmetry_matrix/SMatfunction.cpp submodules/symmetry_matrix/Operator.cpp submodules/symmetry_matrix/SMatMat.cpp
    submodules/FCIdump/FCIdump.cpp
    submodules/IterativeSolver/IterativeSolver/ISDavidson.cpp
    submodules/IterativeSolver/IterativeSolver/ISDiis.cpp
    submodules/IterativeSolver/IterativeSolver/ISRSPT.cpp
    submodules/IterativeSolver/IterativeSolver/IterativeSolver.cpp
    submodules/IterativeSolver/IterativeSolver/PagedParameterVector.cpp
    submodules/IterativeSolver/IterativeSolver/Storage.cpp
    submodules/IterativeSolver/IterativeSolver/SimpleParameterVector.cpp
    )
add_definitions(-DNOMAIN)
add_executable(gci.exe ${SRCS})

#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Eigen3 3.3 REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})
include_directories(src extra submodules/Profiler submodules/FCIdump submodules/IterativeSolver submodules/symmetry_matrix submodules/memory)

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
