#ifndef GCITRANSITIONDENSITY_H
#define GCITRANSITIONDENSITY_H
#include <vector>
#include <gciStringSet.h>
#include <gciOrbitalSpace.h>

namespace gci {

/*!
 * \brief Class to hold transition density matrix,
 * defined by an array of ExcitationSet objects
 */
class TransitionDensity : public std::vector<double>
{
public:
    TransitionDensity();
};
}

using namespace gci;

#endif // GCITRANSITIONDENSITY_H
