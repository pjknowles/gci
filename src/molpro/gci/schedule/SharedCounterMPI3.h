#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H

#include <map>
#include <memory>
#include <mpi.h>

#include "molpro/gci/schedule/SharedCounterBase.h"

namespace molpro::gci::schedule {
/*!
 * @brief Shared counter using MPI3 remote memory access.
 *
 * The counter is stored on a single target processor.
 *
 */
class SharedCounterMPI3 : public SharedCounterBase {
public:
  explicit SharedCounterMPI3(const MPI_Comm &communicator, int target_rank = 0);
  virtual ~SharedCounterMPI3();
  virtual unsigned long int increment(int amount = 1);
  virtual void reset();

protected:
  MPI_Comm m_communicator; //! Counter is accessible by processes in this communicator
  int m_rank;              //! rank of current process
  int m_target_rank;       //! rank of process storing the counter
  MPI_Win m_win;
  size_t *m_win_buffer;
};

} // namespace molpro::gci::schedule
#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERMPI3_H
