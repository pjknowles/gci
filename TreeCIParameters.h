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
 */
    TreeCIParameters();
    ~TreeCIParameters();
/**
 * @brief
 *
 * @param file
 */
    TreeCIParameters(std::string file);
    /**
     * @brief
     *
     * @param file
     */
    void loadParameters(std::string file="FCIDUMP");
    int basisSize; /**< number of orbitals */
    int nelec; /**< number of electrons */
    int ms2; /**< twice the spin quantum number, ie multiplicity minus one */


};

#endif // TREECIPARAMETERS_H
