#ifndef GCI_H
#define GCI_H
#include <vector>

#ifndef MOLPRO
#define xout std::cout
#else
#define xout std::cout
#endif

namespace gci {
class Hamiltonian;
class State;
/*!
 * \brief Perform a Rayleigh-Schroedinger perturbation theory calculation
 * \return the energies order by order
 */
std::vector<double> RSPT( std::vector<gci::Hamiltonian*>& hamiltonians , State& prototype, int maxOrder=4);

void HamiltonianPrint (Hamiltonian& hamiltonian, State& prototype, int verbosity=0);


}

//using namespace gci;

#endif // GCI_H
