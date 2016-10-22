#include "gciOperator.h"

Operator::Operator(Hamiltonian &hamiltonian_, std::vector<int> ranks_)
{
  this->ranks = ranks_;
  this->hamiltonian = hamiltonian_;
}
