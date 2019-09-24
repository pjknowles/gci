unset(DEPENDENCIES)
set(SYMMETRY_MATRIX_DOUBLE ON CACHE BOOL "Whether to build and test symmetry_matrix with double data type" FORCE)
set(SYMMETRY_MATRIX_COMPLEXDOUBLE OFF CACHE BOOL "Whether to build and test symmetry_matrix with complex-double data type" FORCE)
set(SYMMETRY_MATRIX_FLOAT OFF CACHE BOOL "Whether to build and test symmetry_matrix with float data type" FORCE)
if (NOT CMAKE_Fortran_COMPILER)
    message(STATUS "Building GCI without Fortran support")
    set(FORTRAN OFF)
else ()
    message(STATUS "Building GCI with Fortran support")
    set(FORTRAN ON)
endif ()
foreach (dep FCIdump IterativeSolver memory molpro_plugin Profiler symmetry_matrix)
    string(TOLOWER ${dep} depl)
    get_dependency(${depl})
    list(APPEND DEPENDENCIES ${dep})
endforeach ()
if (NOT TARGET ga::ga)
    get_dependency(ga)
    list(APPEND DEPENDENCIES ga)
endif ()
if (Molpro_SOURCE_DIR)
    get_dependency(ppidd)
    list(APPEND DEPENDENCIES ppidd)
endif ()


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
        gciHProduct.cpp gciHProduct.h
        gciHProductSet.cpp gciHProductSet.h
        gciMixedOperator.cpp gciMixedOperator.h
        gciMixedWavefunction.cpp gciMixedWavefunction.h
        gciVibSpace.cpp gciVibSpace.h
        gciDavidson.h gciDavidson.cpp
        gciMixedOperatorSecondQuant.h gciMixedOperatorSecondQuant.cpp
        gciVibOperator.cpp gciVibOperator.h
        gciVibExcitation.h)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "gci.h;gciRun.h;${CMAKE_CURRENT_BINARY_DIR}/gci-config.h")

if (Molpro_SOURCE_DIR)
    target_compile_definitions(gci PRIVATE MOLPRO)
endif ()

find_package(BLAS)
if (BLAS_FOUND)
    messagev("Found BLAS")
    target_compile_definitions(gci PRIVATE -DBLAS)
endif ()
STRING(FIND "${BLAS_LIBRARIES}" "mkl" POS)
if (POS GREATER -1)
    messagev("Found MKL")
    target_compile_definitions(gci PRIVATE -DEIGEN_USE_MKL_ALL)
endif ()

target_link_libraries(gci PUBLIC
        ${BLAS_LINKER_FLAGS} ${BLAS_LIBRARIES})

#target_include_directories(gci PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_link_libraries(gci PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
#target_compile_options(gci PRIVATE ${MPI_CXX_COMPILE_FLAGS})

if (MPI_CXX_FOUND)
    target_link_libraries(gci PUBLIC MPI::MPI_CXX)
endif ()

if (TARGET ga::ga AND FALSE)
    target_link_libraries(gci PUBLIC ga::ga)
endif()

configure_library(${PROJECT_NAME} ${DEPENDENCIES})
