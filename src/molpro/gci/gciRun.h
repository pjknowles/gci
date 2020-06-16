#ifndef GCIRUN_H
#define GCIRUN_H
#include <cstdint>
#include <vector>
#include <molpro/FCIdump.h>
#include <molpro/SMat.h>

#include "molpro/gci/gci.h"
#include "molpro/gci/gciOptions.h"
#include "molpro/gci/gciState.h"
#include "molpro/gci/gciWavefunction.h"

namespace gci {
/*!
 * \brief The Run class encapsulates a complete calculation using gci
 */
class Run {
 public:
  /*!
   * \brief Construct Run object from an FCIdump
   * \param fcidump the file name of the FCIdump
   */
  explicit Run(std::string fcidump);
  ~Run();
  /*!
   * \brief Start the run
   */
  std::vector<double> run();
 private:
  /*!
   * \brief Perform a Rayleigh-Schroedinger perturbation theory calculation
   * \param hams A vector of pointers to the zero, first, second... order hamiltonians
   * \param prototype A State object specifying number of electrons, symmetry, spin
   * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
   * \param maxOrder The maximum order of perturbation theory
   * \param energyThreshold Convergence threshold for energy
   * \param maxIterations The maximum number of iterations to perform
   * \return the energies order by order
   */
  std::vector<double> RSPT(const std::vector<molpro::Operator *> &hams,
                           const State &prototype,
                           int maxOrder = -1,
                           double energyThreshold = -1,
                           int maxIterations = -1);
  std::vector<double> ISRSPT(const molpro::Operator &ham, const molpro::Operator &ham0,
                             const State &prototype,
                             int maxOrder = -1,
                             double energyThreshold = -1,
                             int maxIterations = -1);
  void IPT(const molpro::Operator &ham,
           const State &prototype, size_t referenceLocation);

  /*!
   * \brief Perform a variational calculation using the Davidson algorithm
   * \param ham The Hamiltonian
   * \param prototype A State object specifying number of electrons, symmetry, spin
   * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
   * \param energyThreshold Convergence threshold for energy
   * \param nState The number of eigenstates to be converged
   * \param maxIterations The maximum number of iterations to perform
   * \return the energies for each state. Note that only the first nState energies are considered converged
   */
  std::vector<double> Davidson(const molpro::Operator &ham,
                               const State &prototype,
                               double energyThreshold = (double) -1, int nState = -1, int maxIterations = -1);
  std::vector<double> CSDavidson(const molpro::Operator &ham,
                                 const State &prototype,
                                 double energyThreshold = (double) -1, int nState = -1, int maxIterations = -1);
  /*!
   * \brief Perform a variational calculation using the preconditioned stepest descent algorithm
   * \param ham The Hamiltonian
   * \param prototype A State object specifying number of electrons, symmetry, spin
   * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
   * \param energyThreshold Convergence threshold for energy
   * \param maxIterations The maximum number of iterations to perform
   * \return the energy of the state.
   */
  std::vector<double> DIIS(const molpro::Operator &ham,
                           const State &prototype,
                           double energyThreshold = (double) -1, int maxIterations = -1);

  void HamiltonianMatrixPrint(molpro::Operator &hamiltonian, const State &prototype, int verbosity = 0);

  molpro::Operator m_hamiltonian;
 public:
  std::vector<molpro::Operator> m_densityMatrices; // the individual state density matrices

  std::vector<std::shared_ptr<Wavefunction> > m_wavefunctions;
  gci::Options options;

};

/*!
   * \brief int1 Generate array of diagonal one-electron integrals
   * \param hamiltonian The hamiltonian
   * \param spin positive for alpha, negative for beta
   * \return one-dimensional array with h(i,i) at i-1
   */
Eigen::VectorXd int1(const molpro::Operator& hamiltonian, int spin);

/*!
   * \brief intJ Generate array of two-electron exchange integrals
   * \param hamiltonian The hamiltonian
   * \param spini positive for alpha, negative for beta, first index
   * \param spinj positive for alpha, negative for beta, second index
   * \return array with (ii|jj)
   */
Eigen::MatrixXd intJ(const molpro::Operator& hamiltonian, int spini, int spinj);
/*!
   * \brief intK Generate array of two-electron Coulomb integrals
   * \param hamiltonian The hamiltonian
   * \param spin positive for alpha, negative for beta
   * \return array with (ij|ji)
   */
Eigen::MatrixXd intK(const molpro::Operator& hamiltonian, int spin);



  /*!
   * \brief Construct an Operator from an FCIdump. If the FCIdump
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a FCIdump.
   */
  molpro::Operator constructOperator(const molpro::FCIdump &dump, bool collective=true);
//  inline SymmetryMatrix::Operator constructOperator(FCIdump &&dump) {
//    return constructOperator(dump);
//  }

/*!
 * \brief Build a Fock operator from the density arising from a single Slater determinant
 * \param hamiltonian The hamiltonian
 * \param reference The Slater determinant
 * \param description Descriptive text
 */
molpro::Operator fockOperator(const molpro::Operator& hamiltonian, const Determinant &reference, std::string description = "Fock");

void gsum(molpro::Operator& op);

/*!
 * \brief Build a same-spin operator from the density arising from a single Slater determinant
 * \param reference The Slater determinant
 * \param description Descriptive text
 */
molpro::Operator sameSpinOperator(const molpro::Operator& hamiltonian, const Determinant &reference,
                          std::string description = "Same Spin Hamiltonian");

/*!
 * @brief Write an Operator to an FCIdump
 * @param op The operator
 * @param filename
 * @param orbital_symmetries
 */
void FCIDump(const molpro::Operator& op, const std::string filename, std::vector<int> orbital_symmetries=std::vector<int>(0));

}

#endif // GCIRUN_H
