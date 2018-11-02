#ifndef GCI_GCIBBO_RHF_H
#define GCI_GCIBBO_RHF_H

#include <string>
#include <vector>

#include "gciOperatorBBO.h"
#include "gciRHF.h"

namespace gci {
namespace nm_BBO_RHF {

using nm_RHF::Density;

/*!
 * @brief Write the format for convergence print out during RHF iterations
 */
void writeFormat();

/*!
 * @brief Writes the convergence print out during RHF iterations
 * @param iIter Number of SCF cycle
 * @param energy Energies: electronic, vibrational (per mode), interaction (per mode)
 * @param energyPrev Energies from a previous cycle
 * @param nMode Number of modes
 */
void writeIter(int iIter, std::valarray<double> &energy, std::valarray<double> &energyPrev, int nMode);

/*!
 * @brief Solves the Fock equations, generating a new set of MO's and modals
 * @param molHam Molecular Hamiltonian
 * @param density Electronic density and MO's
 * @param U Vector of Modal coefficients
 */
void solveFock(gci::OperatorBBO &molHam, Density &density, std::vector<SMat> &U);

/*!
 * @brief Rotates the unitary matrix, mixing the different orbitals
 * @param mat Unitary matrix
 * @param ang Rotation constan
 */
void rotate(SMat &mat, double ang);


} // namespace nm_BBO_RHF
} //  namespace gci

#endif //GCI_GCIBBO_RHF_H
