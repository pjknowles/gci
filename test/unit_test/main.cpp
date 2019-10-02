#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mpi.h>
#include <gci.h>
#include <ga.h>
#include <SharedCounter.h>


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int result;
    gci::profiler = std::make_unique<Profiler>(Profiler("GCI"));
    MPI_Init(&argc, &argv);
    GA_Initialize();
    MPI_Comm_rank(MPI_COMM_COMPUTE, &gci::parallel_rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &gci::parallel_size);
    result = RUN_ALL_TESTS();
    gci::profiler.reset();
    gci::_nextval_counter[MPI_COMM_COMPUTE].reset(nullptr);
    GA_Terminate();
    MPI_Finalize();
    return result;
}

