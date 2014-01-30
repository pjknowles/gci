#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"
#include "gciState.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State
{
public:
    double* buffer; ///< buffer to hold coefficients describing the object
private:
};
}
using namespace gci;
#endif // GCIWAVEFUNCTION_H
