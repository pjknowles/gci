cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
enable_language(CXX)
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-O3 -march=native" _native)
if (_native)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
endif()
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
  set(NPROC 1)
endif()
if(NPROC EQUAL 8)
#    hack to avoid overload on workstation / CI
  set(NPROC 4)
endif()
message("Number of processors ${NPROC}")
set(CTEST_BUILD_FLAGS -j${NPROC})

file (GLOB SRCS src/*.cpp
    submodules/Profiler/Profiler.cpp
    submodules/memory/bytestreamC.cpp
    submodules/symmetry_matrix/SMat.cpp submodules/symmetry_matrix/SMatfunction.cpp submodules/symmetry_matrix/Operator.cpp submodules/symmetry_matrix/SMatMat.cpp
    submodules/FCIdump/FCIdump.cpp
    submodules/molpro-plugin/PluginGuest.cpp
    )
add_definitions(-DNOMAIN)
add_executable(gci.exe ${SRCS})
include_directories(src submodules/Profiler submodules/FCIdump submodules/IterativeSolver submodules/symmetry_matrix submodules/memory submodules/molpro-plugin)

find_package (Eigen3 3.3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
target_link_libraries(gci.exe ${MPI_CXX_LIBRARIES})
add_definitions(-DHAVE_MPI_H)

find_package(MKL)
if (MKL_FOUND)
include_directories(${MKL_INCLUDE_DIRS})
add_definitions(-DEIGEN_USE_MKL_ALL)
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
    DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/html/index.xhtml
    )
file(GLOB headers src/*.h)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.xhtml
COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
DEPENDS ${SRC_LIST} ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

enable_testing()
file (GLOB TESTS test/*.fcidump)
foreach(fil ${TESTS})
  get_filename_component(nam ${fil} NAME_WE)
  add_test(${nam} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} --allow-run-as-root ${CMAKE_CURRENT_BINARY_DIR}/gci.exe ${fil})
endforeach()
