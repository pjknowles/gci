#include "gciMolpro.h"
#include "sharedCounter.h"
#include <numeric>

sharedCounter::sharedCounter(const MPI_Comm &communicator)
    : m_hostrank(0), m_rank(0), m_size(1) {
#ifndef SHAREDCOUNTER_DUMMY
  MPI_Comm_size(communicator, &m_size);
  MPI_Comm_rank(communicator, &m_rank);
#endif

  if (m_rank == m_hostrank) {
    m_data.assign((size_t) m_size, 0);
#ifndef SHAREDCOUNTER_DUMMY
    MPI_Win_create(&m_data[0], m_size * sizeof(int), sizeof(int),
                   MPI_INFO_NULL, communicator, &m_win);
#endif
  } else {
#ifndef SHAREDCOUNTER_DUMMY
    MPI_Win_create(&m_data[0], 0, sizeof(int),
                   MPI_INFO_NULL, communicator, &m_win);
#endif
  }
}
sharedCounter::~sharedCounter() {
#ifndef SHAREDCOUNTER_DUMMY
  MPI_Win_free(&m_win);
#endif
}

void sharedCounter::reset() {
  m_myval = 0;
  if (m_rank == m_hostrank)
    m_data.assign((size_t) m_size, 0);
}

#include <iostream>
int sharedCounter::increment(int amount) {
  std::vector<int> vals(static_cast<unsigned long>(m_size));

#ifndef SHAREDCOUNTER_DUMMY
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m_hostrank, 0, m_win);

  for (int i = 0; i < m_size; i++) {
    if (i == m_rank) {
      MPI_Accumulate(&amount, 1, MPI_INT, 0, i, 1, MPI_INT, MPI_SUM,
                     m_win);
    } else {
      MPI_Get(&vals[i], 1, MPI_INT, 0, i, 1, MPI_INT, m_win);
    }
  }

  MPI_Win_unlock(0, m_win);
#endif

  vals[m_rank] = m_myval;
  m_myval += amount;
//  std::cout << "vals[m_rank], m_myval"<<vals[m_rank]<<", "<<m_myval<<std::endl;
  return std::accumulate(vals.begin(), vals.end(), 0); // first returned value is zero
}
