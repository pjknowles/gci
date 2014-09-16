#ifndef GCI_H
#define GCI_H
#include "Profiler.h"

#ifndef MOLPRO
#define xout std::cout
#if defined(GA_MPI) || defined(GA_MPI2)
#define GCI_PARALLEL
#endif
#else
#include <cic/ItfCommon.h>
#include <cic/ItfFortranInt.h>
#endif

#ifdef GCI_PARALLEL
  extern "C" {
#include "ppidd_c.h"
  }
#endif


namespace gci {

extern Profiler profiler; // global profiler

}


#endif // GCI_H
