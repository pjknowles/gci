#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTER_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTER_H

#if GCI_SCHEDULE_SharedCounterGA
#define GCI_SCHEDULE_SHAREDCOUNTER_TYPE SharedCounterGA
#include "molpro/gci/schedule/SharedCounterGA.h"
#elif GCI_SCHEDULE_SharedCounterMPI3
#define GCI_SCHEDULE_SHAREDCOUNTER_TYPE SharedCounterMPI3
#include "molpro/gci/schedule/SharedCounterMPI3.h"
#else
#define GCI_SCHEDULE_SHAREDCOUNTER_TYPE SharedCounterMPI3
#include "molpro/gci/schedule/SharedCounterMPI3.h"
#endif

namespace molpro::gci::schedule {
using SharedCounter = GCI_SCHEDULE_SHAREDCOUNTER_TYPE;
}
#undef GCI_SCHEDULE_SHAREDCOUNTER_TYPE
#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTER_H
