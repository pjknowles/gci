#ifndef GCIPARAMETERS_H
#define GCIPARAMETERS_H
#include "gci.h"
#include <string>
#include <iostream>

namespace gci {
/**
 * @brief
 * Class holds dimension and other parameters for gci basis
 *
 */
class Parameters
{
public:

/**
 * @brief
 *
 * @param filename is the file containing the FCIDUMP. If present, loadParameters is called.
 */
    Parameters(std::string filename="");

    ~Parameters();
    /*!
     \brief
    load basis size, number of electrons, spin from FCIDUMP file.
     \param filename is the file containing the FCIDUMP.
    */
    void load(std::string filename="FCIDUMP"); /**< something */
    unsigned int basisSize; /**< number of orbitals */
    unsigned int nelec; /**< number of electrons */
    int ms2; /**< twice the spin quantum number, ie multiplicity minus one */


};
}

using namespace gci;

#endif // GCIPARAMETERS_H
