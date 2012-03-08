#ifndef TREECIPARAMETERS_H
#define TREECIPARAMETERS_H
#include <string>
#include <iostream>

namespace TreeCI {
/**
 * @brief
 * Class holds dimension and other parameters for FCI basis
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
    unsigned int ms2; /**< twice the spin quantum number, ie multiplicity minus one */


};
}

using namespace TreeCI;

#endif // TREECIPARAMETERS_H
