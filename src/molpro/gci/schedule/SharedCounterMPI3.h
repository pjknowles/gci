#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H

#include <map>
#include <memory>
#include <mpi.h>

namespace molpro {
namespace gci {
/*!
 * @brief Shared counter using MPI3 remote memory access.
 *
 * The counter is stored on a single target processor.
 *
 */
class SharedCounterMPI3 {
public:
  explicit SharedCounterMPI3(const MPI_Comm &communicator, int target_rank = 0);
  SharedCounterMPI3(const SharedCounterMPI3 &) = delete;
  static std::shared_ptr<SharedCounterMPI3> instance(const MPI_Comm &communicator);

  static void clear() { m_counters.clear(); }

  ~SharedCounterMPI3();
  unsigned long int increment(unsigned long int amount = 1);
  //! Resets the counter. Collective.
  void reset();
  unsigned long int myval();

protected:
  MPI_Comm m_communicator;   //! Counter is accessible by processes in this communicator
  unsigned long int m_myval; //! counter incremented by this process only
  int m_rank;                //! rank of current process
  int m_target_rank;         //! rank of process storing the counter
  MPI_Win m_win;
  size_t *m_win_buffer;
  static std::map<MPI_Comm, std::shared_ptr<SharedCounterMPI3>> m_counters;
};

} // namespace gci
} // namespace molpro
#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H
