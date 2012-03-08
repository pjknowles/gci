#ifndef TREECIPARAMETERS_H
#define TREECIPARAMETERS_H
#include <string>
#include <iostream>

/**
 * @brief
 * Class holds dimension and other parameters for FCI basis
 *
 */
class TreeCIParameters
{
public:

/**
 * @brief
 *
 * @param filename is the file containing the FCIDUMP. If present, loadParameters is called.
 */
    TreeCIParameters(std::string filename="");

    ~TreeCIParameters();
    /*!
     \brief
    load basis size, number of electrons, spin from FCIDUMP file.
     \param filename is the file containing the FCIDUMP.
    */
    void load(std::string filename="FCIDUMP"); /**< something */
    int basisSize; /**< number of orbitals */
    int nelec; /**< number of electrons */
    int ms2; /**< twice the spin quantum number, ie multiplicity minus one */


};

#endif // TREECIPARAMETERS_H
