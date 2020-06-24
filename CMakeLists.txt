cmake_minimum_required(VERSION 3.15)
project(gci CXX)

include(FetchContent)
FetchContent_Declare(
        LibraryManager
        GIT_REPOSITORY https://gitlab+deploy-token-178345:RQZa2rmdwnV7Hdwg8XEq@gitlab.com/molpro/LibraryManager.git
        GIT_TAG master
)
FetchContent_MakeAvailable(LibraryManager)
set(MPI ON CACHE BOOL "Whether to build with MPI" FORCE)
option(FORTRAN "Whether to build fortran sources" OFF)
LibraryManager_Project(MPI_OPTION FORTRAN_OPTION)

########## Options ##########
option(GCI_TEST "Whether to force tests to be built" OFF)

########## Dependencies ##########

add_subdirectory(dependencies)

########## Targets ##########

add_subdirectory(src)

if (${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME} OR GCI_TEST)
    enable_testing()
    add_subdirectory(test)
endif ()
########## Documentation ##########

# Documentation has to come last to make sure all dependencies have been added
include(DependencyManagerDocs)
DependencyManagerDocs_Add(${PROJECT_NAME}
        # FILES README.md
        DOC_URL "https://molpro.gitlab.io/gci/${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
        # projects whose documentation has to be built before the current project.
        DEPENDS FCIdump LinearAlgebra memory molpro-plugin Profiler symmetry_matrix utilities
        )

