#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <gci.h>
#include <SharedCounter.h>


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int result;
    gci::profiler = std::make_unique<Profiler>(Profiler("GCI"));
    MPI_Init(&argc, &argv);
    GA_Initialize();
    gci::mpi_comm_compute = GA_MPI_Comm();
    MPI_Comm_rank(gci::mpi_comm_compute, &gci::parallel_rank);
    MPI_Comm_size(gci::mpi_comm_compute, &gci::parallel_size);
    result = RUN_ALL_TESTS();
    gci::profiler.reset();
    gci::_nextval_counter[gci::mpi_comm_compute].reset(nullptr);
    GA_Terminate();
    MPI_Finalize();
    return result;
}

