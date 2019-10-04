
declare_dependency(
        googletest
        https://github.com/google/googletest.git
)

declare_dependency(
        FCIdump
        https://gitlab.com/molpro/fcidump.git
)

declare_dependency(
        IterativeSolver
        https://gitlab+deploy-token-74790:BQ4_qi9fyN7sypf6H__m@gitlab.com/molpro/IterativeSolver.git
)

declare_dependency(
        memory
        git@github.com:molpro/memory
)

declare_dependency(
        molpro_plugin
        https://gitlab.com/molpro/molpro-plugin.git
)

declare_dependency(
        Profiler
        https://gitlab+deploy-token-74792:Bj6FL_cWqz_8y2PNDycx@gitlab.com/molpro/profiler.git
)

declare_dependency(
        symmetry_matrix
        https://gitlab+deploy-token-74793:dLXMspxsDJozbVw7VBRF@gitlab.com/molpro/symmetry_matrix.git
)

find_program(ga-config ga-config)
if (ga-config AND FALSE)
    foreach (f cflags cppflags ldflags libs)
        execute_process(COMMAND ${ga-config} --${f} OUTPUT_VARIABLE ga-${f})
        string(STRIP "${ga-${f}}" ga-${f})
    endforeach ()
    add_library(ga::ga INTERFACE IMPORTED GLOBAL)
    target_link_libraries(ga::ga INTERFACE ${ga-libs})
    target_link_options(ga::ga INTERFACE ${ga-ldflags})
    target_compile_options(ga::ga INTERFACE ${ga-cflags} ${ga-cppflags})
    message("Using local global array installation")
else ()
    set(GA_RUNTIME "MPI_PR" CACHE STRING "overwriting default for GA runtime protocol")
    declare_dependency(
            ga
            https://github.com/GlobalArrays/ga.git
    )
endif ()
