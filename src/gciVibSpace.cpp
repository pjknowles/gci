#include "gciVibSpace.h"

#include <stdexcept>

namespace gci {

VibSpace::VibSpace(int mode, int modal, int modeCouplingLvl) : nMode(mode), nModal(modal),
                                                                    modeCoupling(modeCouplingLvl) {
    if (modeCoupling > nMode) throw std::logic_error("Mode coupling level cannot be > number of modes");
}

bool VibSpace::operator==(const VibSpace &other) const {
    return (nMode == other.nMode) && (nModal == other.nModal) && (modeCoupling == other.modeCoupling);
}

}  // namespace gci
