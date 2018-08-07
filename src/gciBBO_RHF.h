#ifndef GCI_GCIBBO_RHF_H
#define GCI_GCIBBO_RHF_H

#include <string>
#include <vector>

#include "SMat.h"
#include "gciOperator.h"
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
 */
void writeIter(int iIter, std::valarray<double> &energy, std::valarray<double> &energyPrev);

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
