#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <mpi.h>
#include <gci.h>
#include <sharedCounter.h>

using namespace gci;


int gci::parallel_size = 1;
int gci::parallel_rank = 0;
bool gci::molpro_plugin = false;
MPI_Comm gci::molpro_plugin_intercomm = MPI_COMM_NULL;
std::unique_ptr<sharedCounter> gci::_nextval_counter = nullptr;
int64_t gci::__my_first_task = 0;
int64_t gci::__task = 0;
int64_t gci::__task_granularity = 1;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int result;
    profiler = std::make_unique<Profiler>(Profiler("GCI"));
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_COMPUTE, &parallel_rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &parallel_size);
    result = RUN_ALL_TESTS();
    profiler.reset();
    _nextval_counter.reset();
    MPI_Finalize();
    return result;
}

