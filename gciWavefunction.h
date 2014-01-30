#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"
#include "gciState.h"
#include "gciString.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State
{
public:
    std::vector<String> alphaStrings; ///< The alpha-spin strings defining the CI basis
    std::vector<String> betaStrings; ///< The beta-spin strings defining the CI basis
    double* buffer; ///< buffer to hold coefficients describing the object

    void buildStrings(); ///< build alphaStrings and betaStrings
private:
};
}
using namespace gci;
#endif // GCIWAVEFUNCTION_H
