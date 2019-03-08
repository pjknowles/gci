cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)

project(gci)
message("Build ${CMAKE_PROJECT_NAME} with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

########## Options ##########

# Allow messages from scripts to be shown
set(FETCHCONTENT_QUIET OFF)


########## Dependencies ##########

include(FetchContent)
set(SYMMETRY_MATRIX_DOUBLE ON CACHE BOOL "Whether to build and test symmetry_matrix with double data type" FORCE)
set(SYMMETRY_MATRIX_COMPLEXDOUBLE OFF CACHE BOOL "Whether to build and test symmetry_matrix with complex-double data type" FORCE)
set(SYMMETRY_MATRIX_FLOAT OFF CACHE BOOL "Whether to build and test symmetry_matrix with float data type" FORCE)
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

find_package(MPI REQUIRED CXX )
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
    file(GLOB headers src/*.h lib/*.h)
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.xhtml
            COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            DEPENDS ${SRC_LIST} ${headers} ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM
            )
endif (DOXYGEN_FOUND)

