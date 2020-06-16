include(FetchContent)
FetchContent_Declare(
        dependency_manager
        GIT_REPOSITORY https://gitlab.com/dependencymanager/dependency-manager.git
        GIT_TAG 0.1.5
)
FetchContent_MakeAvailable(dependency_manager)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)

DependencyManager_Declare(FCIdump https://gitlab.com/molpro/fcidump.git)
DependencyManager_Declare(LinearAlgebra https://gitlab+deploy-token-186248:BVoGyrUG_o7zpBfCcsLD@gitlab.com/molpro/linearalgebra.git)
DependencyManager_Declare(memory git@github.com:molpro/memory)
DependencyManager_Declare(molpro-plugin https://gitlab.com/molpro/molpro-plugin.git)
DependencyManager_Declare(Profiler https://gitlab+deploy-token-74792:Bj6FL_cWqz_8y2PNDycx@gitlab.com/molpro/profiler.git)
DependencyManager_Declare(symmetry_matrix https://gitlab+deploy-token-74793:dLXMspxsDJozbVw7VBRF@gitlab.com/molpro/symmetry_matrix.git)
DependencyManager_Declare(utilities https://gitlab+deploy-token-181030:whWoM4qkenzjGjWsFTxj@gitlab.com/molpro/utilities.git)

find_program(ga-config ga-config)
set(GA_RUNTIME "MPI_RMA" CACHE STRING "overwriting default for GA runtime protocol")
find_program(ga-config ga-config)
if (ga-config)
    foreach (f cflags cppflags ldflags libs)
        execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
        string(STRIP "${ga-${f}}" ga-${f})
    endforeach ()
    add_library(ga::ga INTERFACE IMPORTED GLOBAL)
    target_link_libraries(ga::ga INTERFACE ${ga-libs})
    target_link_options(ga::ga INTERFACE ${ga-ldflags})
    target_compile_options(ga::ga INTERFACE ${ga-cflags} ${ga-cppflags})
else (ga-config)
    include(FetchContent)
    FetchContent_Declare(
            ga
            GIT_REPOSITORY https://github.com/GlobalArrays/ga
            GIT_TAG 2518e23433385bfa3726d507b8cd0d7ed038021b
    )
    FetchContent_MakeAvailable(ga)
endif (ga-config)

set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG "release-1.8.1"
)
