#ifndef GCI_GCIRHF_H
#define GCI_GCIRHF_H

#include "gciOperator.h"
#include <vector>

namespace gci{
namespace nm_RHF{

/*!
 * @brief Utility class for storing the MO coefficients and constructing density operator
 */
class Density {
public:
    dim_t dim; //! full dimensionality of MO's
    dim_t occ; //! Occupied orbitals
    std::vector<int> symmetries; //! symmetry of each orbital
    int symmetry; //! symmetry of the target state
    SMat Cmat; //! MO coefficients
    SMat Csplice; //! Occupied MO's
    Operator P;

    Density(dim_t &dim, dim_t &occ, std::vector<int> &symmetries, int symmetry);

    ~Density() = default;

    /*!
     * @brief Updates the density matrix. Should be called after Density::Cmat was modified.
     */
    void update();
};

/*!
 * @brief Expectation value of electronic Hamiltonian over HF wavefunction
 * @param P Electronic density matrix
 * @param Hel Electronic Hamiltonian
 * @param energy HF energy
 * @return Fock operator
 */
Operator electronicEnergy(Operator &P, Operator &Hel, double &energy);

} // namespace nm_RHF
} // namespace gci

#endif //GCI_GCIRHF_H
