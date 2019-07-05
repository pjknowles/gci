include(FetchContent)
function(get_dependency MAKESUBDIR)
    foreach (dep IN LISTS ARGN)
        FetchContent_GetProperties(${dep})
        if (NOT ${dep}_POPULATED)
            file(LOCK ${CMAKE_SOURCE_DIR}/dependencies/${dep}_lockfile GUARD FILE TIMEOUT 1000)
            FetchContent_Populate(${dep})
            if (MAKESUBDIR)
                add_subdirectory(${${dep}_SOURCE_DIR} ${${dep}_BINARY_DIR} EXCLUDE_FROM_ALL)
            endif ()
            file(LOCK ${CMAKE_SOURCE_DIR}/dependencies/${dep}_lockfile RELEASE)
        endif ()
        foreach (s SOURCE_DIR BINARY_DIR POPULATED)
            set(${dep}_${s} "${${dep}_${s}}" PARENT_SCOPE)
        endforeach ()
    endforeach ()
endfunction()


if (NOT TARGET Eigen3::Eigen)
    get_dependency(OFF eigen3_dep)
    add_library(ExternalEigen3 INTERFACE)
    target_include_directories(ExternalEigen3 INTERFACE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/dependencies/Eigen3>
            $<INSTALL_INTERFACE:include>
            )
    add_library(Eigen3::Eigen ALIAS ExternalEigen3)
    message("Using external Eigen3 library")
    list(APPEND DEPENDENCIES ExternalEigen3)
else ()
    list(APPEND DEPENDENCIES Eigen3::Eigen)
endif ()

set(SYMMETRY_MATRIX_DOUBLE ON CACHE BOOL "Whether to build and test symmetry_matrix with double data type" FORCE)
set(SYMMETRY_MATRIX_COMPLEXDOUBLE OFF CACHE BOOL "Whether to build and test symmetry_matrix with complex-double data type" FORCE)
set(SYMMETRY_MATRIX_FLOAT OFF CACHE BOOL "Whether to build and test symmetry_matrix with float data type" FORCE)
if (NOT CMAKE_Fortran_COMPILER)
    message("Building GCI without Fortran support")
    set(FORTRAN OFF)
else ()
    message("Building GCI with Fortran support")
    set(FORTRAN ON)
endif ()
foreach (dep fcidump_dep;iterativesolver_dep;memory_dep;molpro_plugin_dep;profiler_dep;symmetry_matrix_dep)
    get_dependency(ON ${dep})
endforeach ()
list(APPEND DEPENDENCIES Profiler FCIdump molpro_plugin)
#list(APPEND DEPENDENCIES Profiler FCIdump IterativeSolver memory molpro_plugin)


add_library(gci
        gci.h
        gciDeterminant.cpp gciDeterminant.h
        gciExcitationSet.cpp gciExcitationSet.h
        gciFile.cpp gciFile.h
        gciMolpro.cpp gciMolpro.h
        gciOptions.cpp gciOptions.h
        gciOrbitals.h
        gciOrbitalSpace.cpp gciOrbitalSpace.h
        gciPrintable.h
        gciRun.cpp gciRun.h
        gciState.cpp gciState.h
        gciString.cpp gciString.h
        gciStringSet.cpp gciStringSet.h
        gciSymmetrySpace.cpp gciSymmetrySpace.h
        gciTransitionDensity.cpp gciTransitionDensity.h
        gciWavefunction.cpp gciWavefunction.h
        sharedCounter.cpp sharedCounter.h
        )

if (Molpro_SOURCE_DIR)
    target_compile_definitions(gci PRIVATE MOLPRO)
endif ()

find_package(BLAS)
if (BLAS_FOUND)
    message("Found BLAS")
    target_compile_definitions(gci PRIVATE -DBLAS)
endif ()
STRING(FIND "${BLAS_LIBRARIES}" "mkl" POS)
if (POS GREATER -1)
    message("Found MKL")
    target_compile_definitions(gci PRIVATE -DEIGEN_USE_MKL_ALL)
endif ()

target_link_libraries(gci PUBLIC Eigen3::Eigen Profiler FCIdump IterativeSolver symmetry_matrix memory molpro_plugin ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})

target_include_directories(gci PRIVATE ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(gci PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_compile_options(gci PRIVATE ${MPI_CXX_COMPILE_FLAGS})

if (MPI_CXX_FOUND)
    target_link_libraries(gci PUBLIC MPI::MPI_CXX)
endif ()

include(libconfig.cmake)

