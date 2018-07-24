#ifndef GCIRUN_H
#define GCIRUN_H
#include <vector>
#include "gci.h"
#include "gciOperator.h"
#include "gciState.h"
#include "gciWavefunction.h"
#include "gciOptions.h"
#include <cstdint>

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
  std::vector<double> RSPT(const std::vector<Operator *> &hams,
                           const State &prototype,
                           int maxOrder = -1,
                           double energyThreshold = -1,
                           int maxIterations = -1);
  std::vector<double> ISRSPT(const gci::Operator &ham, const gci::Operator &ham0,
                             const State &prototype,
                             int maxOrder = -1,
                             double energyThreshold = -1,
                             int maxIterations = -1);
  void IPT(const gci::Operator &ham,
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
  std::vector<double> Davidson(const Operator &ham,
                               const State &prototype,
                               double energyThreshold = (double) -1, int nState = -1, int maxIterations = -1);
  std::vector<double> CSDavidson(const Operator &ham,
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
  std::vector<double> DIIS(const Operator &ham,
                           const State &prototype,
                           double energyThreshold = (double) -1, int maxIterations = -1);

  /*!
   * @brief Performs RHF calculation within the space of supplied MOs
   * @param hamiltonian The hamiltonian
   * @param prototype A State object specifying number of electrons, symmetry, spin
   * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
   * @param threshold Convergence threshold
   * @param maxIterations The maximum number of iterations to perform
   * @return the energy of the state.
   */
  double RHF(const Operator &hamiltonian, const State &prototype, double thresh=1.0e-4, int maxIterations=20);

  void HamiltonianMatrixPrint(Operator &hamiltonian, const State &prototype, int verbosity = 0);

  Operator m_hamiltonian;
 public:
  std::unique_ptr<Operator> m_densityMatrix; // the (state-averaged) density matrix
  std::vector<Operator> m_densityMatrices; // the individual state density matrices

  std::vector<std::shared_ptr<Wavefunction> > m_wavefunctions;
  gci::Options options;
};
}

using namespace gci;

#endif // GCIRUN_H
