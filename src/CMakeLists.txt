###############################
### Library of the project
###############################

add_library(gci-lib
        gciDeterminant.cpp gciDeterminant.h
        gciExcitationSet.cpp gciExcitationSet.h
        gciFile.cpp gciFile.h
        gciMolpro.cpp gciMolpro.h
        gciOperator.cpp gciOperator.h
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
        ${CMAKE_SOURCE_DIR}/submodules/FCIdump/FCIdump.cpp
        ${CMAKE_SOURCE_DIR}/submodules/IterativeSolver/IterativeSolver.cpp
        ${CMAKE_SOURCE_DIR}/submodules/memory/bytestreamC.cpp
        ${CMAKE_SOURCE_DIR}/submodules/molpro-plugin/PluginGuest.cpp
        ${CMAKE_SOURCE_DIR}/submodules/Profiler/Profiler.cpp
        ${CMAKE_SOURCE_DIR}/submodules/symmetry_matrix/SMat.cpp
        ${CMAKE_SOURCE_DIR}/submodules/symmetry_matrix/SMatfunction.cpp
        ${CMAKE_SOURCE_DIR}/submodules/symmetry_matrix/Operator.cpp
        ${CMAKE_SOURCE_DIR}/submodules/symmetry_matrix/SMatMat.cpp
        )
target_include_directories(gci-lib
        PUBLIC
        ${CMAKE_SOURCE_DIR}/submodules/Profiler
        ${CMAKE_SOURCE_DIR}/submodules/FCIdump
        ${CMAKE_SOURCE_DIR}/submodules/IterativeSolver
        ${CMAKE_SOURCE_DIR}/submodules/symmetry_matrix
        ${CMAKE_SOURCE_DIR}/submodules/memory
        ${CMAKE_SOURCE_DIR}/submodules/molpro-plugin
        INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

### Dependencies

target_link_libraries(gci-lib PUBLIC Eigen3::Eigen)

target_include_directories(gci-lib PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(gci-lib PUBLIC ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
target_compile_options(gci-lib PUBLIC ${MPI_CXX_COMPILE_FLAGS})

if (MKL_FOUND)
    target_include_directories(gci-lib PRIVATE ${MKL_INCLUDE_DIRS})
    target_link_libraries(gci-lib PRIVATE ${MKL_LIBRARIES})
else (MKL_FOUND)
    if (BLAS_FOUND)
        target_include_directories(gci-lib PRIVATE ${BLAS_INCLUDE_DIRS})
        target_link_libraries(gci-lib PUBLIC ${BLAS_LIBRARY} -lm -ldl)
    else (BLAS_FOUND)
        target_link_libraries(gci-lib PUBLIC -llapack -lm -ldl)
    endif (BLAS_FOUND)
endif (MKL_FOUND)

###############################
### Main executable
###############################

add_executable(gci.exe gci.cpp gci.h)
target_link_libraries(gci.exe PUBLIC gci-lib)
