#ifndef GCI_GCIRHF_H
#define GCI_GCIRHF_H

#include <Operator.h>
#include <vector>

namespace gci{
namespace nm_RHF{

/*!
 * @brief Utility class for storing the MO coefficients and constructing density operator
 */
class Density {
public:
    SymmetryMatrix::dim_t dim; //! full dimensionality of MO's
    SymmetryMatrix::dim_t occ; //! Occupied orbitals
    std::vector<int> symmetries; //! symmetry of each orbital
    int symmetry; //! symmetry of the target state
    SymmetryMatrix::SMat Cmat; //! MO coefficients
    SymmetryMatrix::SMat Csplice; //! Occupied MO's
    SymmetryMatrix::Operator P;

    Density(SymmetryMatrix::dim_t &dim, SymmetryMatrix::dim_t &occ, std::vector<int> &symmetries, int symmetry);

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
SymmetryMatrix::Operator electronicEnergy(const SymmetryMatrix::Operator &P, const SymmetryMatrix::Operator &Hel, double &energy);

} // namespace nm_RHF
} // namespace gci

#endif //GCI_GCIRHF_H
