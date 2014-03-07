#include "gciOperator.h"

Operator::Operator(Hamiltonian &hamiltonian, std::vector<int> ranks)
{
    this->ranks = ranks;
    this->hamiltonian = hamiltonian;
}
