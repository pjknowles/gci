
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
set(GA_RUNTIME "MPI_PR" CACHE STRING "overwriting default for GA runtime protocol")
declare_dependency(
        ga
        https://github.com/GlobalArrays/ga.git
)
