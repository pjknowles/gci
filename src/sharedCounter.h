#ifndef SHAREDCOUNTER_H
#define SHAREDCOUNTER_H
#ifdef MOLPRO
#include "common/molpro_config.h"
#endif
#if !defined(MOLPRO) || defined(HAVE_MPI_H)
#include "mpi.h"
#ifdef MOLPRO
#include "ppidd.h"
#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
#define MPI_COMM_COMPUTE MPI_COMM_WORLD
#endif
#else
#define SHAREDCOUNTER_DUMMY
#define MPI_Comm int
#define MPI_Win int
#define MPI_COMM_WORLD 0
#endif
#include <vector>
#include <cstddef>


class sharedCounter
{
public:
  sharedCounter(const MPI_Comm& communicator=MPI_COMM_COMPUTE);
  ~sharedCounter();
  int increment(int amount=1);
  void reset();
private:
  MPI_Win m_win;
  int  m_hostrank ;
  int  m_myval;
  std::vector<int> m_data;
  int m_rank, m_size;
};

#endif // SHAREDCOUNTER_H
