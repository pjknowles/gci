#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "Operator.h"
#include <utility>
#include "FCIdump.h"
#include "Eigen/Dense"
#include "gciDeterminant.h"
#include "gciOrbitalSpace.h"

namespace gci {

class Operator : public SymmetryMatrix::Operator {
 public:
  /*!
   * \brief Construct a real hermitian fermionic operator
   * \param dimension Numbers of orbitals in each symmetry block.
   * \param orbital_symmetries Symmetries of orbitals, in the range 1 to 8
   * \param rank The rank (0, 1 or 2) of the operator.
   * \param uhf Whether the underlying 1-particle spaces are different for alpha and beta spin.
   * \param symmetry Overall symmetry of operator (0-7).
   * \param covariant Whether the operator transforms like spatial basis functions, rather than like a density matrix
   * \param hermitian Whether the representation of the operator is hermitian.
   * \param description A string describing the object.
   */
  explicit Operator(const dim_t &dimension,
                    int rank = 2,
                    bool uhf = false,
                    unsigned int symmetry = 0,
                    bool covariant = true,
                    bool hermitian = true,
                    std::string description = "")
      : SymmetryMatrix::Operator(dims_t{dimension, dimension, dimension, dimension},
                                 rank,
                                 uhf,
                                 std::vector<int>{hermitian ? 1 : 0, hermitian ? 1 : 0},
                                 std::vector<int>{-1, -1},
                                 symmetry,
                                 covariant,
                                 std::move(description)) {
  }

  Operator(const SymmetryMatrix::Operator &source)
      : SymmetryMatrix::Operator(source) {}

  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a bytestream produced by bytestream()
   */
  static Operator construct(const char *dump);
  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param bs The bytestream produced by bytestream()
   */
  static Operator construct(const class bytestream &bs) { return construct(&(bs.data()[0])); }
  /*!
   * \brief Construct an object from an FCIdump. If the FCIdump
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a FCIdump.
   */
  static Operator construct(const FCIdump &dump);
  static Operator construct(FCIdump &&dump) {
    return construct(dump);
  }


  /*!
     * \brief int1 Generate array of diagonal one-electron integrals
     * \param spin positive for alpha, negative for beta
     * \return one-dimensional array with h(i,i) at i-1
     */
  Eigen::VectorXd int1(int spin) const;

  /*!
     * \brief intJ Generate array of two-electron exchange integrals
     * \param spini positive for alpha, negative for beta, first index
     * \param spinj positive for alpha, negative for beta, second index
     * \return array with (ii|jj)
     */
  Eigen::MatrixXd intJ(int spini, int spinj) const;
  /*!
     * \brief intK Generate array of two-electron Coulomb integrals
     * \param spin positive for alpha, negative for beta
     * \return array with (ij|ji)
     */
  Eigen::MatrixXd intK(int spin) const;

};

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
}

#endif // GCIOPERATOR_H
