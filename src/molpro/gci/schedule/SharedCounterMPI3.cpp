#include "SharedCounterMPI3.h"
#include <iostream>

namespace molpro::gci::schedule {
SharedCounterMPI3::SharedCounterMPI3(const MPI_Comm &communicator, int target_rank)
    : m_communicator(communicator), m_rank(0), m_target_rank(target_rank) {
  int size;
  MPI_Comm_size(m_communicator, &size);
  MPI_Comm_rank(m_communicator, &m_rank);
  if (target_rank < 0 || target_rank >= size)
    throw std::runtime_error("target rank not in the communicator");
  MPI_Aint n = m_rank == target_rank ? 1 : 0;
  auto sizeoftype = sizeof(unsigned long int);
  MPI_Win_allocate(n * sizeoftype, sizeoftype, MPI_INFO_NULL, m_communicator, &m_win_buffer, &m_win);
  SharedCounterMPI3::reset();
}

SharedCounterMPI3::~SharedCounterMPI3() {
  MPI_Win_free(&m_win);
  //    if (m_rank == m_target_rank)
  //        MPI_Free_mem(m_win_buffer);
}

void SharedCounterMPI3::reset() {
  MPI_Win_fence(0, m_win);
  m_curr_val = 0;
  if (m_rank == m_target_rank)
    *m_win_buffer = 0;
  MPI_Win_fence(0, m_win);
}

unsigned long int SharedCounterMPI3::increment(int amount) {
  MPI_Win_lock(MPI_LOCK_SHARED, m_target_rank, 0, m_win);
  unsigned long int a = amount;
  MPI_Fetch_and_op(&a, &m_curr_val, MPI_UNSIGNED_LONG, m_target_rank, 0, MPI_SUM, m_win);
  MPI_Win_unlock(m_target_rank, m_win);
  return m_curr_val;
}
} // namespace molpro::gci::schedule
