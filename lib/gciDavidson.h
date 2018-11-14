#ifndef GCI_GCIDAVIDSON_H
#define GCI_GCIDAVIDSON_H

#include <vector>
#include <IterativeSolver.h>

#include "gci.h"
#include "gciOptions.h"
#include "gciIterativeSolverUtil.h"

/*
 * Now that there are two wavefunction objects and two operators (pure electronic and mixed). There is no
 * point in having a Run class.
 * Each member function runs independently and they should be kept independent as templated functions in a common
 * namespace. I'll write a simple example to illustrate my point using with Davidson's algorithm.
 *
 * The `main` will initiate the calculation by reading in options and calling run::drive function.
 */

namespace gci {
namespace run {


//FIXME A seperate options class is not necessary anymore. Absorb it into Davidson.
/*!
 * @brief Exposes key options for the Davidson's diagonalisation algorithm. Also provides access to the global options
 * class.
 * @note Design structure I have in my head assumes that there is only one options class and it is available through
 * out the calculation at the top level `main` routine.
 *
 * Releveant FCIdump keys:
 *      TOL
 *          -- energy threshold for convergence
 *      NSTATE
 *          --
 *      MAXIT
 *          --
 *      PSPACE_INITIAL
 *          --
 *      PSPACE
 *          --
 *      PSPACE_REBUILD
 *          --
 *      SOLVER_VERBOSITY
 *          --
 */
class DavidsonOptions {
public:
    explicit DavidsonOptions(Options &options) :
            externalOptions(options),
            energyThreshold(options.parameter("TOL", 1e-8)),
            nState(options.parameter("NSTATE", 1)),
            maxIterations(options.parameter("MAXIT", 1000)),
            pSpaceInitial(options.parameter("PSPACE_INITIAL", nState)),
            pSpace(options.parameter("PSPACE", 100)),
            solverVerbosity(options.parameter("SOLVER_VERBOSITY", 1)),
            pSpaceRebuild(options.parameter("PSPACE_REBUILD", 0)) {
        pSpaceInitial = std::max(pSpaceInitial, nState);
        pSpace = std::max(pSpace, nState);
    }

    /*!
     * @brief Gives access to the full set of options
     * @return
     */
    Options &global() const {return externalOptions;}

    double energyThreshold;
    int nState;
    int maxIterations;
    int pSpaceInitial;
    int pSpace;
    int pSpaceRebuild;
    int solverVerbosity;
private:
    Options &externalOptions;
};

//FIXME once DavidsonOptions is merged into Davidson make this a seperate function in Davidson.
/*!
 * @brief Writes the settings of calculation.
 */
std::ostream &operator<<(std::ostream &os, DavidsonOptions const &obj);

/*!
 * @note Copied from gci::Run::Davidson and transformed into a template
 * \brief Perform a variational calculation using the Davidson algorithm
 *
 * @warning The following are my attempts to understand Peter's algorithm. If anything is inaccurate, please correct it!!!
 *
 * As far as I can make out, this is Davidson's algorithm with an extended parameter space added at the start of the
 * calculation. This should improve convergence in a strongly coupled system, by including most relevant determinants.
 * When searching for a few lowest excited states, providing their approximations (i.e. single determinant) should
 * also improve their convergence. More importantly if the excited state is strictly orthogonal to the ground state guess
 * than in normal Davidson it would not be included, however by adding it into the parameter space it can be found.
 *
 * \param ham The Hamiltonian
 * \param prototype A State object specifying number of electrons, symmetry, spin
 * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
 * \param energyThreshold Convergence threshold for energy
 * \param nState The number of eigenstates to be converged
 * \param maxIterations The maximum number of iterations to perform
 * \return the energies for each state. Note that only the first nState energies are considered converged
 */
template<class t_Wavefunction, class t_Operator>
class Davidson {
public:
    using value_type = typename t_Wavefunction::value_type;

    /*!
     * @brief Sets up the linear eigenvalue problem
     * @param prototype Prototype wavefunction; it defines the Fock space for the problem.
     * @param ham Hamiltonian operator that is to be diagonalised
     * @param options Global options defining the state of the whole calculation. Options relevant to the calculation
     *        are extracted into DavidsonOptions member.
     */
    explicit Davidson(const std::shared_ptr<t_Wavefunction> &prototype, const std::shared_ptr<t_Operator> &ham,
                      Options &options);

    /*!
     * @brief Runs the calculation.
     */
    void run();

    DavidsonOptions options;//!< Options governing the type of calculation
    std::shared_ptr<t_Wavefunction> prototype; //!< Prototype wavefunction
    std::shared_ptr<t_Operator> ham; //!< Hamiltonian operator that is being diagonalized
    std::shared_ptr<t_Wavefunction> diagonalH; //!< Stored diagonal values of the hamiltonian operator
    LinearAlgebra::LinearEigensystem<t_Wavefunction> solver; //!< Iterative solver
    std::vector<value_type> eigVal; //!< Solution eigenvalues
    std::vector<t_Wavefunction> eigVec; //!< Solution eigenvectors
protected:
    /*!
     * @brief Initializes containers necessary for running the calculation
     */
    void initialize();
    /*!
     * @brief Builds up the P space
     */
    void buildPspace();
    int maxNP; //!<
    int initialNP;//!<
    updater<t_Wavefunction> update;//!< Applies the preconditioner to the working solutions, based on the residual
    residual<t_Wavefunction, t_Operator> resid;//!<
    std::vector<bool> active; //!< Active components of ww contribute to expanding the P space
    ParameterVectorSet<t_Wavefunction> gg; //!< Set of residual vectors
    ParameterVectorSet<t_Wavefunction> ww; //!< Set of current solutions
    std::vector<std::vector<double>> Pcoeff; //!< Target states represented in the P space
    std::vector<Pvector> P;
    Presidual<t_Wavefunction, t_Operator> Presid;
};


}  // namespace run
}  // namespace gci

#endif //GCI_GCIDAVIDSON_H
