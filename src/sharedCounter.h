#ifndef SHAREDCOUNTER_H
#define SHAREDCOUNTER_H
#if !defined(MOLPRO) || defined(GA_MPI) || defined(MPI2)
#include "mpi.h"
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
  sharedCounter(MPI_Comm communicator=MPI_COMM_WORLD);
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
