#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <macdecls.h>
#include <molpro/gci/gci.h>
#include <molpro/gci/schedule/SharedCounterGA.h>

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  GA_Initialize();
  int mem = 10000000;
  MA_init(C_DBL, mem, mem);
  molpro::gci::mpi_comm_compute = GA_MPI_Comm();
  MPI_Comm_rank(molpro::gci::mpi_comm_compute, &molpro::gci::parallel_rank);
  MPI_Comm_size(molpro::gci::mpi_comm_compute, &molpro::gci::parallel_size);
  molpro::gci::profiler = std::make_unique<molpro::Profiler>(molpro::Profiler("GCI"));
  molpro::gci::_nextval_counter[molpro::gci::mpi_comm_compute].reset(nullptr);
  if (!GA_Create_mutexes(1))
    GA_Error((char *)"Failed to create mutexes", 1);
  int result = RUN_ALL_TESTS();
  GA_Destroy_mutexes();
  molpro::gci::profiler.reset();
  GA_Terminate();
  MPI_Finalize();
  return result;
}
