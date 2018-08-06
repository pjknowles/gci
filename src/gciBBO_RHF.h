#ifndef GCI_GCIBBO_RHF_H
#define GCI_GCIBBO_RHF_H

#include <string>
#include <vector>

#include "SMat.h"
#include "gciOperator.h"
#include "gciOperatorBBO.h"

namespace gci {
namespace nm_BBO_RHF {

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
 * @brief Write the format for convergence print out during RHF iterations
 */
void writeFormat();

/*!
 * @brief Writes the convergence print out during RHF iterations
 */
void writeIter(int iIter, double eTot, double eTotPrev, std::vector<double> &energy,
               std::vector<double> &energyPrev);

/*!
 * @brief Solves the Fock equations, generating a new set of MO's and modals
 * @param molHam Molecular Hamiltonian
 * @param density Electronic density and MO's
 * @param U Vector of Modal coefficients
 */
void solveFock(gci::OperatorBBO &molHam, Density &density, std::vector<SMat> &U);

} // namespace nm_BBO_RHF
} //  namespace gci

#endif //GCI_GCIBBO_RHF_H
