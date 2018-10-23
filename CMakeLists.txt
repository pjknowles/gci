cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
enable_language(CXX)

###############################
### Options
###############################
add_definitions(-DNOMAIN) # Do we want it global? Thish should really be a variable passed to a specific target

###############################
### Global Compiler Flags
###############################
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-O3 -march=native" _native)
if (_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")
endif ()
include(ProcessorCount)
ProcessorCount(NPROC)
if (NPROC EQUAL 0)
    set(NPROC 1)
endif ()
if (NPROC EQUAL 8)
    #    hack to avoid overload on workstation / CI
    set(NPROC 4)
endif ()
message("Number of processors ${NPROC}")
set(CTEST_BUILD_FLAGS -j${NPROC})


###############################
### Dependencies
###############################

find_package(Eigen3 3.3 REQUIRED)

find_package(MPI REQUIRED)
add_definitions(-DHAVE_MPI_H)                                                       # Should add this to target
# OpenMPI root guard: https://github.com/open-mpi/ompi/issues/4451
if ("$ENV{USER}" STREQUAL "root" OR "$ENV{HOME}" STREQUAL "/root")
    set(MPIOPTIONS "--allow-run-as-root") # used in tests
endif ()

find_package(MKL)
if (MKL_FOUND)
    add_definitions(-DEIGEN_USE_MKL_ALL)
    message("MKL link libraries are ${MKL_LIBRARIES}")
else (MKL_FOUND)
    find_package(BLAS REQUIRED)
    if (BLAS_FOUND)
        message("Blas link library is ${BLAS_LIBRARIES}")
    endif (BLAS_FOUND)
endif (MKL_FOUND)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.8.1"
)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(${PROJECT_NAME}-doc ALL
            DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.xhtml
            )
    file(GLOB headers src/*.h)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.xhtml
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SRC_LIST} ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)

###############################
### Targets
###############################

add_subdirectory(src)


###############################
### Tests
###############################

add_subdirectory(test)
