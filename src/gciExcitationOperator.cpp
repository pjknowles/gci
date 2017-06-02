#include "gciExcitationOperator.h"

ExcitationOperator::ExcitationOperator(OldOperator &hamiltonian_, std::vector<int> ranks_)
{
  this->ranks = ranks_;
  this->hamiltonian = hamiltonian_;
}
