#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
#if !defined(MOLPRO_NO_LONGER) // TODO MOLPRO might not have global arrays

#include <mpi.h>

#include <ga-mpi.h>
#include <ga.h>
#include <molpro/mpi.h>

#else
#define SHAREDCOUNTER_DUMMY
//#define MPI_Comm int
//#define MPI_Win int
//#define MPI_COMM_COMPUTE 0
#endif

#include "molpro/gci/schedule/SharedCounterBase.h"
#include <map>
#include <memory>

namespace molpro::gci::schedule {
/*!
 * @brief Shared counter using Global Arrays
 */
class SharedCounterGA : public SharedCounterBase {
public:
  explicit SharedCounterGA(const MPI_Comm &communicator);
  virtual ~SharedCounterGA();
  virtual unsigned long int increment(int amount = 1);
  virtual void reset();

protected:
  MPI_Comm m_communicator; //! MPI communicator whose processes share this counter
  int m_hostrank;
  int m_ga_handle;
  int m_ga_pgroup;
  int m_rank; //! rank of process in communicator
  int m_size; //! size of communicator
};

} // namespace molpro::gci::schedule
#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
