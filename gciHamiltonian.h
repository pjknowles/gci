#ifndef GCIHAMILTONIAN_H
#define GCIHAMILTONIAN_H
#include "gci.h"
#include "FCIdump.h"
#include <string>
#include <vector>
using namespace std;

namespace gci {
/**
 * @brief
 * Class holds hamiltonian operator for gci or other calculation
 *
 */
class Hamiltonian
{
public:

/*!
 \brief

 \param filename : if present, call load
*/
    Hamiltonian(string filename="");
    Hamiltonian(FCIdump* dump);
    ~Hamiltonian();
    void load(string filename="FCIDUMP"); /**< load integrals from FCIDUMP */
    void load(FCIdump* dump); /**< load integrals from FCIDUMP */
    void unload(); /**< destroy loaded integrals */
    bool loaded;  /**< whether the integrals are loaded */
    bool spinUnrestricted; /**< whether alpha and beta spin orbitals are different */
    double coreEnergy; /**< core energy */
    vector<double> *integrals_a;  /**< point to aa integrals */
    vector<double> *integrals_b; /**< point to bb integrals */
    vector<double> *integrals_aa; /**< point to aaaa integrals */
    vector<double> *integrals_ab; /**< point to aabb integrals */
    vector<double> *integrals_bb; /**< point to bbbb integrals */
    unsigned int basisSize;
private:
    unsigned int ijSize;
    unsigned int ijklSize;
    int verbosity;
    vector<unsigned int> orbital_symmetries;
    vector<unsigned int> symmetry_dimensions;
};
}

using namespace gci;

#endif // GCIHAMILTONIAN_H
