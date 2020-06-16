#ifndef SHAREDCOUNTER_H
#define SHAREDCOUNTER_H
#ifdef MOLPRO
#include "molpro_config.h"
#endif
#if !defined(MOLPRO)

#include <mpi.h>

#ifdef MOLPRO
#include "ppidd.h"
//#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
//#define MPI_COMM_COMPUTE MPI_COMM_WORLD
#endif
#else
#define SHAREDCOUNTER_DUMMY
#define MPI_Comm int
#define MPI_Win int
//#define MPI_COMM_COMPUTE 0
#endif

#include <cstddef>
#include <vector>

namespace gci {
class SharedCounter {
public:
  explicit SharedCounter(const MPI_Comm &communicator);
  SharedCounter(const SharedCounter &) = delete;
  ~SharedCounter();
  int increment(int amount = 1);
  void reset();

private:
  MPI_Comm m_communicator;
  int m_hostrank;
  long int m_myval;
  int m_ga_handle;
  int m_ga_pgroup;
  int m_rank;
  int m_size;
};

} // namespace gci
#endif // SHAREDCOUNTER_H
