#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "Operator.h"
#include <utility>
#include "FCIdump.h"
#include "Eigen/Dense"
#include "gciDeterminant.h"
#include "gciOrbitalSpace.h"

namespace gci {

  /*!
   * \brief Construct an Operator from an FCIdump. If the FCIdump
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a FCIdump.
   */
  SymmetryMatrix::Operator constructOperator(const FCIdump &dump);
  inline SymmetryMatrix::Operator constructOperator(FCIdump &&dump) {
    return constructOperator(dump);
  }

/*!
 * \brief Build a Fock operator from the density arising from a single Slater determinant
 * \param hamiltonian The hamiltonian
 * \param reference The Slater determinant
 * \param description Descriptive text
 */
SymmetryMatrix::Operator fockOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference, std::string description = "Fock");

void gsum(SymmetryMatrix::Operator& op);

/*!
 * \brief Build a same-spin operator from the density arising from a single Slater determinant
 * \param reference The Slater determinant
 * \param description Descriptive text
 */
SymmetryMatrix::Operator sameSpinOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference,
                          std::string description = "Same Spin Hamiltonian");

/*!
 * \brief Construct an operator templated on source, but with a special specification
 * \param special
 * - "Q" projector onto space containing satellite orbital
 * - "P" 1-Q
 * \param forceSpinUnrestricted whether to force conversion to a UHF object
 */
SymmetryMatrix::Operator *projector(const SymmetryMatrix::Operator& source, std::string special, bool forceSpinUnrestricted) ;

/*!
 * @brief Write an Operator to an FCIdump
 * @param op The operator
 * @param filename
 * @param orbital_symmetries
 */
void FCIDump(const SymmetryMatrix::Operator& op, const std::string filename, std::vector<int> orbital_symmetries=std::vector<int>(0));

/*!
   * \brief int1 Generate array of diagonal one-electron integrals
   * \param hamiltonian The hamiltonian
   * \param spin positive for alpha, negative for beta
   * \return one-dimensional array with h(i,i) at i-1
   */
Eigen::VectorXd int1(const SymmetryMatrix::Operator& hamiltonian, int spin);

/*!
   * \brief intJ Generate array of two-electron exchange integrals
   * \param hamiltonian The hamiltonian
   * \param spini positive for alpha, negative for beta, first index
   * \param spinj positive for alpha, negative for beta, second index
   * \return array with (ii|jj)
   */
Eigen::MatrixXd intJ(const SymmetryMatrix::Operator& hamiltonian, int spini, int spinj);
/*!
   * \brief intK Generate array of two-electron Coulomb integrals
   * \param hamiltonian The hamiltonian
   * \param spin positive for alpha, negative for beta
   * \return array with (ij|ji)
   */
Eigen::MatrixXd intK(const SymmetryMatrix::Operator& hamiltonian, int spin);

}

#endif // GCIOPERATOR_H
