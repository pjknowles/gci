#ifndef SHAREDCOUNTER_H
#define SHAREDCOUNTER_H
#include "mpi.h"
#include <vector>


class sharedCounter
{
public:
  sharedCounter(MPI_Comm communicator=MPI_COMM_WORLD);
  ~sharedCounter();
  int increment(int amount=1);
private:
  MPI_Win m_win;
  int  m_hostrank ;
  int  m_myval;
  std::vector<int> m_data;
  int m_rank, m_size;
};

#endif // SHAREDCOUNTER_H
