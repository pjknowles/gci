#include "gciVibSpace.h"

#include <stdexcept>

namespace gci {

VibSpace::VibSpace(int mode, int modal, int excLevel) : nMode(mode), nModal(modal), excLvl(excLevel) {
  if (excLvl > nMode)
    throw std::logic_error("Mode coupling level cannot be > number of modes");
}

bool VibSpace::operator==(const VibSpace &other) const {
  return (nMode == other.nMode) && (nModal == other.nModal) && (excLvl == other.excLvl);
}

} // namespace gci
