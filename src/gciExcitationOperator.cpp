#include "gciExcitationOperator.h"

ExcitationOperator::ExcitationOperator(Hamiltonian &hamiltonian_, std::vector<int> ranks_)
{
  this->ranks = ranks_;
  this->hamiltonian = hamiltonian_;
}
