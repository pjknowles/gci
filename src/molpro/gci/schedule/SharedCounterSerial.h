#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERSERIAL_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERSERIAL_H

#include <molpro/gci/schedule/SharedCounterBase.h>

namespace molpro{
namespace gci{
namespace schedule{

/*!
 * @brief A dummy shared counter that does not use MPI and works in serial only.
 */
class SharedCounterSerial : public SharedCounterBase {};
}
}
}

#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERSERIAL_H
