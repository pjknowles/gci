#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction
{
public:
    Wavefunction(unsigned int Nelec=0,int MS2=0);
    void setNelec(unsigned int Nelec); //< set the number of electrons
    void setMS2(int MS2); //< set the Sz quantum number
    unsigned int Nelec(); //< report the number of electrons
    int MS2();//< report the Sz quantum number * 2
    double* buffer;
private:
    Hamiltonian hamiltonian; //< the underlying second-quantized hamiltonian operator
    unsigned int nelec; //< number of electrons
    int ms2; //< twice the spin quantum number, ie multiplicity minus one
};
}
using namespace gci;
#endif // GCIWAVEFUNCTION_H
