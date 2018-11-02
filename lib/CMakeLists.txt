########## Dependencies ##########

FetchContent_GetProperties(fcidump_dep)
if (NOT fcidump_dep_POPULATED)
    FetchContent_Populate(fcidump_dep)
    add_subdirectory(${fcidump_dep_SOURCE_DIR} ${fcidump_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
FetchContent_GetProperties(iterativesolver_dep)
if (NOT iterativesolver_dep_POPULATED)
    FetchContent_Populate(iterativesolver_dep)
    add_subdirectory(${iterativesolver_dep_SOURCE_DIR} ${iterativesolver_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
FetchContent_GetProperties(memory_dep)
if (NOT memory_dep_POPULATED)
    FetchContent_Populate(memory_dep)
    add_subdirectory(${memory_dep_SOURCE_DIR} ${memory_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
FetchContent_GetProperties(molpro_plugin_dep)
if (NOT molpro_plugin_dep_POPULATED)
    FetchContent_Populate(molpro_plugin_dep)
    add_subdirectory(${molpro_plugin_dep_SOURCE_DIR} ${molpro_plugin_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
FetchContent_GetProperties(profiler_dep)
if (NOT profiler_dep_POPULATED)
    FetchContent_Populate(profiler_dep)
    add_subdirectory(${profiler_dep_SOURCE_DIR} ${profiler_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()
FetchContent_GetProperties(symmetry_matrix_dep)
if (NOT symmetry_matrix_dep_POPULATED)
    FetchContent_Populate(symmetry_matrix_dep)
    add_subdirectory(${symmetry_matrix_dep_SOURCE_DIR} ${symmetry_matrix_dep_BINARY_DIR} EXCLUDE_FROM_ALL)
endif ()


########## Target ##########

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
        gciRHF.cpp gciRHF.h
        gciBBO_RHF.cpp gciBBO_RHF.h
        gciHProduct.cpp gciHProduct.h
        gciHProductSet.cpp gciHProductSet.h
        gciMixedOperator.cpp gciMixedOperator.h
        gciMixedWavefunction.cpp gciMixedWavefunction.h
        gciOperatorBBO.cpp gciOperatorBBO.h
        gciVibSpace.cpp gciVibSpace.h
        )
target_include_directories(gci INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(gci PUBLIC Eigen3::Eigen Profiler FCIdump IterativeSolver symmetry_matrix memory molpro_plugin)

target_include_directories(gci PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(gci PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_compile_options(gci PUBLIC ${MPI_CXX_COMPILE_FLAGS})

if (MKL_FOUND)
    target_include_directories(gci PRIVATE ${MKL_INCLUDE_DIRS})
    target_link_libraries(gci PRIVATE ${MKL_LIBRARIES})
else ()
    if (BLAS_FOUND)
        target_include_directories(gci PRIVATE ${BLAS_INCLUDE_DIRS})
        target_link_libraries(gci PUBLIC ${BLAS_LIBRARY} -lm -ldl)
    else ()
        target_link_libraries(gci PUBLIC -llapack -lm -ldl)
    endif ()
endif ()

if (MPI_CXX_FOUND)
    target_link_libraries(gci PUBLIC MPI::MPI_CXX)
endif ()


