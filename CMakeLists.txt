cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(gci)
enable_language(CXX)
#TODO Need to split dependency libraries into CXX and Fortran  components and only compile the one that's needed
enable_language(Fortran)

########## Options ##########

# Allow messages from scripts to be shown
set(FETCHCONTENT_QUIET OFF)


########## Global Compiler Flags ##########

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
#    set(NPROC 4)
endif ()
message("Number of processors ${NPROC}")
set(CTEST_BUILD_FLAGS -j${NPROC})


########## Dependencies ##########

include(FetchContent)
set(SYMMETRY_MATRIX_DOUBLE ON CACHE BOOL "Whether to build and test symmetry_matrix with double data type" FORCE)
set(SYMMETRY_MATRIX_COMPLEXDOUBLE OFF CACHE BOOL "Whether to build and test symmetry_matrix with complex-double data type" FORCE)
set(SYMMETRY_MATRIX_FLOAT OFF CACHE BOOL "Whether to build and test symmetry_matrix with float data type" FORCE)
set(SYMMETRY_MATRIX_FORTRAN OFF CACHE BOOL "Whether to build symmetry_matrix fortran binding" FORCE)
add_subdirectory(dependencies)

if (NOT TARGET Eigen3::Eigen)
    FetchContent_GetProperties(eigen3_dep)
    if (NOT eigen3_dep_POPULATED)
        FetchContent_Populate(eigen3_dep)
        add_library(ExternalEigen3 INTERFACE)
        target_include_directories(ExternalEigen3 INTERFACE ${CMAKE_SOURCE_DIR}/dependencies/Eigen3)
        add_library(Eigen3::Eigen ALIAS ExternalEigen3)
        message("Using external Eigen3 library")
    endif ()
endif ()

find_package(MPI REQUIRED CXX Fortran)
#TODO Consider reflecting this mandatory provision of MPI in the code, ie remove #ifdef HAVE_MPI_H
add_definitions(-DHAVE_MPI_H)

########## Targets ##########

add_subdirectory(lib)
add_subdirectory(src)
add_subdirectory(test)


########## Documentation ##########

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

