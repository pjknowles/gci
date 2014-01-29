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
 \brief construct Hamiltonian object

 \param filename : if present, call load
*/
    Hamiltonian(string filename="");
    /*!
     \brief construct Hamiltonian object

     \param dump : if present, call load
    */   Hamiltonian(FCIdump* dump);
    ~Hamiltonian();
    void load(string filename="FCIDUMP"); /**< \brief load integrals from FCIDUMP */
    void load(FCIdump* dump); /**< \brief load integrals from FCIDUMP */
    void unload(); /**< \brief destroy loaded integrals */
    /*!
     * \brief Construct a printable representation of the hamiltonian
     * \param verbosity how much information to include
     * \return printable representation of the hamiltonian
     */
    std::string printable(int verbosity=0);
    bool loaded;  /**< \brief whether the integrals are loaded */
    bool spinUnrestricted; /**< \brief whether alpha and beta spin orbitals are different */
    double coreEnergy; /**< \brief core energy */
    vector<double> *integrals_a;  /**< \brief point to aa integrals */
    vector<double> *integrals_b; /**< \brief point to bb integrals */
    vector<double> *integrals_aa; /**< \brief point to aaaa integrals */
    vector<double> *integrals_ab; /**< \brief point to aabb integrals */
    vector<double> *integrals_bb; /**< \brief point to bbbb integrals */
    unsigned int basisSize;///< \brief size of orbital basis set
private:
    unsigned int ijSize;
    unsigned int ijklSize;
    int verbosity;
    vector<int> orbital_symmetries;
    vector<unsigned int> symmetry_dimensions;
};
}

using namespace gci;

#endif // GCIHAMILTONIAN_H
