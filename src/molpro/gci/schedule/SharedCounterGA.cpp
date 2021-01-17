#include "SharedCounterGA.h"

#include "molpro/gci/gci.h"

namespace molpro::gci::schedule {
SharedCounterGA::SharedCounterGA(const MPI_Comm& communicator)
    : m_communicator(communicator), m_hostrank(0), m_ga_handle(0), m_rank(0), m_size(1) {
  int glob_size, glob_rank;
  MPI_Comm_size(m_communicator, &m_size);
  MPI_Comm_rank(m_communicator, &m_rank);
  MPI_Comm_size(mpi_comm_compute, &glob_size);
  MPI_Comm_rank(mpi_comm_compute, &glob_rank);
  auto glob_ranks = std::vector<int>(m_size, 0);
  MPI_Allgather(&glob_rank, 1, MPI_INT, glob_ranks.data(), 1, MPI_INT, m_communicator);
  // create new processor group from the communicator
  if (m_size > 1) {
    m_ga_pgroup = GA_Pgroup_create(glob_ranks.data(), m_size);
    m_ga_handle = NGA_Create_handle();
    NGA_Set_pgroup(m_ga_handle, m_ga_pgroup);
    int dims = 1, chunk = 1;
    NGA_Set_data(m_ga_handle, 1, &dims, C_LONG);
    NGA_Set_array_name(m_ga_handle, (char*)"SharedCounter");
    NGA_Set_chunk(m_ga_handle, &chunk);
    auto succ = GA_Allocate(m_ga_handle);
    if (!succ)
      GA_Error((char*)"Failed to allocate", 0);
    GA_Check_handle(m_ga_handle, (char*)"Failed in SharedCounter constructor");
  }
}

SharedCounterGA::~SharedCounterGA() {
  if (m_size > 1)
    GA_Destroy(m_ga_handle);
}

void SharedCounterGA::reset() {
  m_curr_val = 0;
  if (m_size > 1) {
    GA_Check_handle(m_ga_handle, (char*)"Failed in SharedCounter::reset()");
    GA_Zero(m_ga_handle);
  }
}

unsigned long int SharedCounterGA::increment(int amount) {
  m_curr_val += amount;
  if (m_size > 1) {
    int subscript = 0;
    m_curr_val = NGA_Read_inc(m_ga_handle, &subscript, (long int)amount);
  }
  return m_curr_val;
}
} // namespace molpro::gci::schedule
