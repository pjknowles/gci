#ifndef GCI_H
#define GCI_H
#include "Profiler.h"

#ifndef MOLPRO
#define xout std::cout
#else
#include <cic/ItfCommon.h>
#include <cic/ItfFortranInt.h>
#endif



namespace gci {

extern Profiler profiler; // global profiler

}


#endif // GCI_H
