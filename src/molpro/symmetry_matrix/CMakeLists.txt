# Source files at this level define the public interface to the library
LibraryManager_Append(${PROJECT_NAME}
        SOURCES
        SMat.cpp
        SMat-implementation.h
        SMatMat-implementation.h
        Operator-implementation.h Operator-double.cpp
        PUBLIC_HEADER Operator.h SMat.h SMatMat.h
        )
if (SYMMETRY_MATRIX_DOUBLE)
    LibraryManager_Append(${PROJECT_NAME} SOURCES SMat-double.cpp SMatMat-double.cpp Operator-double.cpp)
endif ()
if (SYMMETRY_MATRIX_COMPLEXDOUBLE)
    message(STATUS "Building for complex double")
    LibraryManager_Append(${PROJECT_NAME} SOURCES SMat-complex-double.cpp SMatMat-complex-double.cpp Operator-complex-double.cpp)
endif ()
if (SYMMETRY_MATRIX_FLOAT)
    LibraryManager_Append(${PROJECT_NAME} SOURCES SMat-float.cpp SMatMat-float.cpp Operator-float.cpp)
endif ()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_14)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

include(DependencyManager)
#DependencyManager_Populate(memory)
#target_link_libraries(${PROJECT_NAME} PUBLIC molpro::memory)
DependencyManager_Populate(utilities)
target_link_libraries(${PROJECT_NAME} PUBLIC molpro::utilities)

LibraryManager_LAPACK(${PROJECT_NAME} Intel10_64lp_seq Intel10_64lp)
if (MKL_TYPE STREQUAL "lp64")
    message(STATUS "Eigen will use MKL")
    target_compile_definitions(${PROJECT_NAME} PUBLIC EIGEN_USE_MKL_ALL)
endif ()
if (FORTRAN AND NOT LAPACK_FOUND)
    find_package(LAPACK REQUIRED)
    target_link_libraries(${PROJECT_NAME} PUBLIC ${LAPACK_LIBRARIES})
endif()

if (NOT TARGET Eigen3::Eigen)
    message(VERBOSE "Eigen3::Eigen is not yet a target")
    find_package(Eigen3 3.3.0)
    if (EIGEN3_FOUND)
        message(STATUS "Eigen3 found on system")
    else ()
        if (NOT FORTRAN)
            message(FATAL_ERROR "Eigen3 cannot be configured from source with FORTRAN=OFF; install it in the system first")
        endif ()
        include(FetchContent)
        FetchContent_Declare(eigen3
                GIT_REPOSITORY https://gitlab.com/libeigen/eigen
                GIT_TAG 3.3.7
                )
        FetchContent_MakeAvailable(eigen3)
        add_library(Eigen3::Eigen ALIAS eigen)
        message(STATUS "Eigen3::Eigen not found on system, and downloaded")
    endif ()
else ()
    message(STATUS "SymmetryMatrix: Eigen3::Eigen target already defined")
endif ()

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_MPL2_ONLY)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_DONT_VECTORIZE)
endif ()



