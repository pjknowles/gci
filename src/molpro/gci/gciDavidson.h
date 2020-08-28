#ifndef GCI_GCIDAVIDSON_H
#define GCI_GCIDAVIDSON_H

#include <molpro/Operator.h>
#include <molpro/linalg/IterativeSolver.h>
#include <vector>

#include "molpro/gci/gci.h"
#include "molpro/gci/gciOptions.h"
#include "molpro/gci/gciWavefunction.h"

/*
 * Now that there are two wavefunction objects and two operators (pure electronic and mixed). There is no
 * point in having a Run class.
 * Each member function runs independently and they should be kept independent as templated functions in a common
 * namespace. I'll write a simple example to illustrate my point using with Davidson's algorithm.
 *
 * The `main` will initiate the calculation by reading in options and calling run::drive function.
 */

namespace molpro {
namespace gci {
namespace run {

/*!
 * @note Copied from gci::Run::Davidson and transformed into a template
 * \brief Perform a variational calculation using the Davidson algorithm
 *
 *
 * @brief Exposes key options for the Davidson's diagonalisation algorithm. Also provides access to the global options
 * class.
 * @note Design structure I have in my head assumes that there is only one options class and it is available through
 * out the calculation at the top level `main` routine.
 *
 * Releveant FCIdump keys:
 *      TOL     -- Energy threshold for convergence
 *      NSTATE  -- Number of target states
 *      MAXIT   -- Maximum number of iterations
 *      PSPACE_INITIAL  -- Size of the initial parameter sapce
 *      PSPACE  -- Maximum size of the parameter space
 *      PSPACE_REBUILD  -- Rebuild P space at every iteration
 *      SOLVER_VERBOSITY -- Verbosity of the output from solver (see IterativeSolver)
 *
 * \param ham The Hamiltonian
 * \param prototype A State object specifying number of electrons, symmetry, spin
 * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default
 * FCI \param energyThreshold Convergence threshold for energy \param nState The number of eigenstates to be converged
 * \param maxIterations The maximum number of iterations to perform
 * \return the energies for each state. Note that only the first nState energies are considered converged
 */
template <class t_Wavefunction, class t_Operator> class Davidson {
public:
  using ParameterVectorSet = std::vector<t_Wavefunction>;
  using value_type = typename t_Wavefunction::value_type;

  /*!
   * @brief Sets up the linear eigenvalue problem
   * @param prototype Prototype wavefunction; it defines the Fock space for the problem.
   * @param ham Hamiltonian operator that is to be diagonalised
   * @param options Global options defining the state of the whole calculation. Options relevant to the calculation
   *        are extracted into DavidsonOptions member.
   */
  explicit Davidson(const t_Wavefunction &prototype, const t_Operator &_ham, Options options);

  /*!
   * @brief Runs the calculation.
   */
  void run();

  t_Wavefunction prototype;                  //!< Prototype wavefunction
  const t_Operator &ham;                     //!< Hamiltonian operator that is being diagonalized
  Options options;                           //!< Options governing the type of calculation
  std::shared_ptr<t_Wavefunction> diagonalH; //!< Stored diagonal values of the hamiltonian operator
                                             //    std::vector<value_type> eigVal; //!< Solution eigenvalues
                                             //    std::vector<t_Wavefunction> eigVec; //!< Solution eigenvectors
  ParameterVectorSet ww;                     //!< Set of current solutions
  ParameterVectorSet gg;                     //!< Set of residual vectors
protected:
  std::shared_ptr<std::vector<Wavefunction>> ref_elec_states; //!< Set of reference electronic states
  void printMatrix(const std::string &fname) const;
  void message() const;
  //! Initializes containers necessary for running the calculation
  void initialize();
  //! Prepares an initial guess. For Mixed wavefunction it entails solving the electronic CI problem
  void prepareGuess();
  //! Prepares reference electronic states
  void reference_electronic_states();
  //! Run analysis of solutions
  void analysis();
  //! Seperate energetic contributions from different parts of the Hamiltonian
  void energy_decomposition();
  //! Apply the Hamiltonian on the current solution
  void action();
  //! Get the new vector, $r = A u - \lambda u$
  void update();
  //! Store the current solutions in a backup file
  void backup(std::vector<t_Wavefunction> &ww);
  // TODO reinstate the following, which doesn't compile until handler implementation for Wavefunction is done
//  molpro::linalg::LinearEigensystem<t_Wavefunction> solver; //!< Iterative solver
  molpro::linalg::LinearEigensystem<std::vector<double>> solver; //!< Iterative solver

  double energyThreshold;
  unsigned int nState;
  unsigned int maxIterations;
  int solverVerbosity;
  int parallel_stringset;
  std::string restart_file;
  std::string backup_file;
  std::vector<size_t> diag_minlocN;
  std::vector<double> diag_val_at_minlocN;
};

template <class t_Wavefunction>
void davidson_read_write_wfn(std::vector<t_Wavefunction> &ww, const std::string &fname, bool save);

} // namespace run
} // namespace gci
} // namespace molpro

#endif // GCI_GCIDAVIDSON_H
