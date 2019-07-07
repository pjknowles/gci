cmake_minimum_required(VERSION 3.11)
set(CMAKE_CXX_STANDARD 14)

project(gci VERSION 1.0)
option(FORTRAN "Whether to build with Fortran support" OFF)
message("Build ${CMAKE_PROJECT_NAME} with CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

########## Options ##########

# Allow messages from scripts to be shown
#set(FETCHCONTENT_QUIET OFF)


########## Dependencies ##########

add_subdirectory(dependencies)

find_package(MPI REQUIRED CXX )
#TODO Consider reflecting this mandatory provision of MPI in the code, ie remove #ifdef HAVE_MPI_H
add_definitions(-DHAVE_MPI_H)

########## Targets ##########

add_subdirectory(lib)
add_subdirectory(src)
enable_testing()
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

