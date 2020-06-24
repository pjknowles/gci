#ifndef SHAREDCOUNTER_H
#define SHAREDCOUNTER_H
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

#include <map>
#include <memory>

namespace molpro {
namespace gci {
namespace schedule {
class SharedCounterGA {
public:
  explicit SharedCounterGA(const MPI_Comm &communicator);
  SharedCounterGA(const SharedCounterGA &) = delete;
  ~SharedCounterGA();
  static std::shared_ptr<SharedCounterGA> instance(const MPI_Comm &communicator);
  static void clear() { m_counters.clear(); }
  int increment(int amount = 1);
  void reset();

protected:
  MPI_Comm m_communicator;
  int m_hostrank;
  long int m_myval;
  int m_ga_handle;
  int m_ga_pgroup;
  int m_rank;
  int m_size;
  static std::map<MPI_Comm, std::shared_ptr<SharedCounterGA>> m_counters;
};

} // namespace schedule
} // namespace gci
} // namespace molpro
#endif // SHAREDCOUNTER_H
