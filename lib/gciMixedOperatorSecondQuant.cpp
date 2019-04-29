#include "gciMixedOperatorSecondQuant.h"
#include "gciRun.h"

namespace gci {

MixedOperatorSecondQuant::MixedOperatorSecondQuant(const FCIdump &fcidump) :
        nMode(fcidump.parameter("NMODE", std::vector<int>{0})[0]), Hel(constructOperator(fcidump)) { }
} // namespace gci
