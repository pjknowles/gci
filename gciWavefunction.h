#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"
#include "gciState.h"
#include "gciStringSet.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State
{
public:
    /**
 * @brief
 *
 * @param filename is the file containing the FCIDUMP. If present, load is called.
 */
    Wavefunction(std::string filename="");
    /**
     * @brief
     *
     * @param dump points to an FCIdump object. If present, load is called.
     */
    Wavefunction(FCIdump* dump);

    /*!
     * \brief Construct a Wavefunction object linked to a Hamiltonian
     * \param h The hamiltonian
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
    Wavefunction(Hamiltonian *h, int nelec, int symmetry, int ms2);

    StringSet alphaStrings[8]; ///< The alpha-spin strings defining the CI basis
    StringSet betaStrings[8]; ///< The beta-spin strings defining the CI basis
    double* buffer; ///< buffer to hold coefficients describing the object

    void buildStrings(); ///< build alphaStrings and betaStrings
    /*!
     \brief
    printable form of the String.
     \return std::string
    */

private:
};
}
using namespace gci;
#endif // GCIWAVEFUNCTION_H
