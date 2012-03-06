#ifndef TREECIHAMILTONIAN_H
#define TREECIHAMILTONIAN_H
#include <string>
using namespace std;
/**
 * @brief
 *
 */
class TreeCIHamiltonian
{
public:
/**
 * @brief
 *
 */
    TreeCIHamiltonian();
    void loadHamiltonian(string filename="FCIDUMP"); /**< load integrals from FCIDUMP */
    int basisSize; /**< TODO */
    bool spinUnrestricted; /**< TODO */
    double coreEnergy; /**< TODO */
    double *integrals_a, *integrals_b, *integrals_aa, *integrals_ab, *integrals_bb; /**< TODO */
};

#endif // TREECIHAMILTONIAN_H
