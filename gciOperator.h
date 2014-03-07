#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H

#include "gciHamiltonian.h"

namespace gci {
class Operator
{
public:
    Operator(Hamiltonian &hamiltonian, std::vector<int> ranks);
    std::vector<int> ranks;
    Hamiltonian hamiltonian;
};
}

#endif // GCIOPERATOR_H
