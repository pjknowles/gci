FetchContent_Declare(
        eigen3_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/Eigen3"
        GIT_REPOSITORY https://github.com/eigenteam/eigen-git-mirror.git
        GIT_TAG "3.3.7"
        GIT_SHALLOW TRUE
)

FetchContent_Declare(
        googletest_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/googletest"
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.8.1"
        GIT_SHALLOW TRUE
)

file(STRINGS "FCIdump_SHA1" FCIdump_GIT_TAG)
FetchContent_Declare(
        fcidump_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/FCIdump"
#        GIT_REPOSITORY git@gitlab.com:molpro/fcidump.git
        GIT_REPOSITORY https://gitlab.com/molpro/fcidump.git
        GIT_TAG ${FCIdump_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh FCIdump ${CMAKE_BUILD_TYPE}
)

file(STRINGS "IterativeSolver_SHA1" IterativeSolver_GIT_TAG)
FetchContent_Declare(
        iterativesolver_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/IterativeSolver"
#        GIT_REPOSITORY https://gitlab+deploy-token-74714:VcwszaF11Jm3UJCadMFVVcwszaF11Jm3UJCadMFV@gitlab.com/molpro/iterativesolver.git
        GIT_REPOSITORY git@gitlab.com:molpro/iterativesolver.git
        GIT_TAG ${IterativeSolver_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh IterativeSolver ${CMAKE_BUILD_TYPE}
)

file(STRINGS "memory_SHA1" memory_GIT_TAG)
FetchContent_Declare(
        memory_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/memory"
        GIT_REPOSITORY https://gitlab.com/molpro/memory.git
        GIT_TAG ${memory_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh memory ${CMAKE_BUILD_TYPE}
)

file(STRINGS "molpro-plugin_SHA1" molpro_plugin_GIT_TAG)
FetchContent_Declare(
        molpro_plugin_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/molpro-plugin"
        GIT_REPOSITORY https://gitlab.com/molpro/molpro-plugin.git
        GIT_TAG ${molpro_plugin_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh molpro-plugin ${CMAKE_BUILD_TYPE}
)

file(STRINGS "Profiler_SHA1" Profiler_GIT_TAG)
FetchContent_Declare(
        profiler_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/Profiler"
        GIT_REPOSITORY https://gitlab.com/molpro/profiler.git
        GIT_TAG ${Profiler_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh Profiler ${CMAKE_BUILD_TYPE}
)

file(STRINGS "symmetry_matrix_SHA1" symmetry_matrix_GIT_TAG)
FetchContent_Declare(
        symmetry_matrix_dep
        SOURCE_DIR "${CMAKE_SOURCE_DIR}/dependencies/symmetry_matrix"
        GIT_REPOSITORY https://gitlab.com/molpro/symmetry_matrix.git
        GIT_TAG ${symmetry_matrix_GIT_TAG}
        UPDATE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/checkout.sh symmetry_matrix ${CMAKE_BUILD_TYPE}
)

