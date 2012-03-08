#ifndef TREECIHAMILTONIAN_H
#define TREECIHAMILTONIAN_H
#include <string>
#include <vector>
#include "TreeCIParameters.h"
using namespace std;
/**
 * @brief
 * Class holds hamiltonian operator for FCI or other calculation
 *
 */
class TreeCIHamiltonian : public TreeCIParameters
{
public:

/*!
 \brief

 \param filename : if present, call load
*/
    TreeCIHamiltonian(string filename="");
    ~TreeCIHamiltonian();
    void load(string filename="FCIDUMP"); /**< load integrals from FCIDUMP */
    void unload(); /**< destroy loaded integrals */
    bool loaded;  /**< whether the integrals are loaded */
    bool spinUnrestricted; /**< whether alpha and beta spin orbitals are different */
    double coreEnergy; /**< core energy */
    vector<double> *integrals_a;  /**< point to aa integrals */
    vector<double> *integrals_b; /**< point to bb integrals */
    vector<double> *integrals_aa; /**< point to aaaa integrals */
    vector<double> *integrals_ab; /**< point to aabb integrals */
    vector<double> *integrals_bb; /**< point to bbbb integrals */
private:
    int ijSize;
    int ijklSize;
    int verbosity;
};

#endif // TREECIHAMILTONIAN_H
