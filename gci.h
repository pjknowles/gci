#ifndef GCI_H
#define GCI_H
#include <vector>

#ifndef MOLPRO
#define xout std::cout
#else
#include <cic/ItfCommon.h>
//#define xout std::cout
#endif

namespace gci {
class Hamiltonian;
class State;
/*!
 * \brief Perform a Rayleigh-Schroedinger perturbation theory calculation
 * \param hamiltonians A vector of pointers to the zero, first, second... order hamiltonians
 * \param prototype A State object specifying number of electrons, symmetry, spin
 * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
 * \param energyThreshold Convergence threshold for energy
 * \param maxOrder The maximum order of perturbation theory
 * \return the energies order by order
 */
std::vector<double> RSPT(const std::vector<Hamiltonian *> &hamiltonians,
                         const State& prototype,
                         double energyThreshold=1e-8, int maxOrder=20);

/*!
 * \brief Perform a variational calculation
 * \param hamiltonians A vector of pointers to the zero, first, second... order hamiltonians
 * \param prototype A State object specifying number of electrons, symmetry, spin
 * If it is a Wavefunction , then that will be used also for specifying the configuration space, otherwise the default FCI
 * \param energyThreshold Convergence threshold for energy
 * \param nState The number of eigenstates to be converged
 * \return the energies for each state. Note that only the first nState energies are considered converged
 */
std::vector<double> Davidson(const Hamiltonian &hamiltonian,
                             const State& prototype,
                             double energyThreshold=1e-8, int nState=1, int maxIterations=100);

void HamiltonianMatrixPrint (Hamiltonian& hamiltonian, const State &prototype, int verbosity=0);


}

//using namespace gci;

#endif // GCI_H
