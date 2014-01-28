#ifndef gciPARAMETERS_H
#define gciPARAMETERS_H
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
 * @param filename is the file containing the gciDUMP. If present, loadParameters is called.
 */
    Parameters(std::string filename="");

    ~Parameters();
    /*!
     \brief
    load basis size, number of electrons, spin from gciDUMP file.
     \param filename is the file containing the gciDUMP.
    */
    void load(std::string filename="gciDUMP"); /**< something */
    unsigned int basisSize; /**< number of orbitals */
    unsigned int nelec; /**< number of electrons */
    unsigned int ms2; /**< twice the spin quantum number, ie multiplicity minus one */


};
}

using namespace gci;

#endif // gciPARAMETERS_H
