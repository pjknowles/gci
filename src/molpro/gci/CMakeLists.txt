LibraryManager_Append(${PROJECT_NAME}
        SOURCES
        gciDeterminant.cpp gciExcitationSet.cpp gciFile.cpp gciMolpro.cpp gciOptions.cpp gciOrbitalSpace.cpp gciRun.cpp
        gciState.cpp gciString.cpp gciStringSet.cpp gciSymmetrySpace.cpp gciTransitionDensity.cpp gciWavefunction.cpp
        gciHProduct.cpp gciHProductSet.cpp gciMixedOperator.cpp gciMixedWavefunction.cpp
        gciVibSpace.cpp gciDavidson.cpp gciPersistentOperator.cpp gciMixedOperatorSecondQuant.cpp gciVibOperator.cpp
        gciArray.cpp gciUtils.cpp
        PUBLIC_HEADER
        gci.h gciDeterminant.h gciExcitationSet.h gciFile.h gciMolpro.h gciOptions.h gciOrbitals.h gciOrbitalSpace.h
        gciPrintable.h gciRun.h gciState.h gciString.h gciStringSet.h gciSymmetrySpace.h gciTransitionDensity.h
        gciWavefunction.h gciHProduct.h gciHProductSet.h gciMixedOperator.h gciMixedWavefunction.h
        gciVibSpace.h gciDavidson.h gciPersistentOperator.h gciMixedOperatorSecondQuant.h gciVibOperator.h
        gciVibExcitation.h gciArray.h gciUtils.h)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

add_subdirectory(schedule)

LibraryManager_BLAS(${PROJECT_NAME} PUBLIC Intel10_64lp_seq Intel10_64lp)
target_compile_definitions(${PROJECT_NAME} PRIVATE -DBLAS)

target_link_libraries(${PROJECT_NAME} PUBLIC MPI::MPI_CXX ga::ga)

set(HDF5_PREFER_PARALLEL ON)
#set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED COMPONENTS C)
target_link_libraries(${PROJECT_NAME} PUBLIC ${HDF5_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
target_compile_definitions(${PROJECT_NAME} PUBLIC ${HDF5_DEFINITIONS})

set(SYMMETRY_MATRIX_DOUBLE ON CACHE BOOL "Whether to build and test symmetry_matrix with double data type")
set(SYMMETRY_MATRIX_COMPLEXDOUBLE OFF CACHE BOOL "Whether to build and test symmetry_matrix with complex-double data type")
set(SYMMETRY_MATRIX_FLOAT OFF CACHE BOOL "Whether to build and test symmetry_matrix with float data type")
DependencyManager_Populate(FCIdump)
DependencyManager_Populate(LinearAlgebra)
DependencyManager_Populate(memory)
DependencyManager_Populate(molpro-plugin)
DependencyManager_Populate(Profiler)
DependencyManager_Populate(symmetry_matrix)
DependencyManager_Populate(utilities)
target_link_libraries(${PROJECT_NAME} PUBLIC
        molpro::FCIdump
        molpro::LinearAlgebra
        molpro::memory
        molpro::molpro-plugin
        molpro::Profiler
        molpro::symmetry_matrix
        molpro::utilities
        )
