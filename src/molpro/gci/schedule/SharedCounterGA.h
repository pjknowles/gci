#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
#ifdef MOLPRO
#include "molpro_config.h"
#endif
#if !defined(MOLPRO)

#include <mpi.h>

#include <ga-mpi.h>
#include <ga.h>

#ifdef MOLPRO
#include "ppidd.h"
//#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
//#define MPI_COMM_COMPUTE MPI_COMM_WORLD
#endif
#else
#define SHAREDCOUNTER_DUMMY
#define MPI_Comm int
#define MPI_Win int
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
  virtual int increment(int amount = 1);
  virtual void reset();

protected:
  MPI_Comm m_communicator;
  int m_hostrank;
  long int m_myval;
  int m_ga_handle;
  int m_ga_pgroup;
  int m_rank;
  int m_size;
};

} // namespace molpro::gci::schedule
#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERGA_H
