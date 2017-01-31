#include "sharedCounter.h"
#include <numeric>

sharedCounter::sharedCounter(MPI_Comm communicator)
{
  MPI_Comm_size(communicator,&m_size);
  m_hostrank=0;
  MPI_Comm_rank(communicator,&m_rank);

  if (m_rank == m_hostrank) {
      m_data.assign((size_t)m_size,0);
      MPI_Win_create(&m_data[0], m_size * sizeof(int), sizeof(int),
          MPI_INFO_NULL, communicator, &m_win);
    } else {
      MPI_Win_create(&m_data[0], 0, sizeof(int),
          MPI_INFO_NULL, communicator, &m_win);
    }
  m_myval = 0;
}
sharedCounter::~sharedCounter()
{
  MPI_Win_free(&m_win);
}

int sharedCounter::increment(int amount) {
  std::vector<int> vals(m_size);

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m_hostrank, 0, m_win);

  for (int i=0; i<m_size; i++) {
      if (i == m_rank) {
          MPI_Accumulate(&amount, 1, MPI_INT, 0, i, 1, MPI_INT, MPI_SUM,
                         m_win);
        } else {
          MPI_Get(&vals[i], 1, MPI_INT, 0, i, 1, MPI_INT, m_win);
        }
    }

  MPI_Win_unlock(0, m_win);
  m_myval += amount;

  vals[m_rank] = m_myval;
  return std::accumulate(vals.begin(),vals.end(),0);
}

