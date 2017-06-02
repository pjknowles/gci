#ifndef GCIEXCITATIONOPERATOR_H
#define GCIEXCITATIONOPERATOR_H

#include "gciOperator.h"

namespace gci {
/**
     * @brief
     * A quantum-mechanical operator
     *
     */
class ExcitationOperator
{
public:
  /*!
     * \brief Construct operator object linked to a hamiltonian
     * \param hamiltonian the hamiltonian defining orbital spaces
     * \param ranks Ranks of operator
     */
  ExcitationOperator(OldOperator &hamiltonian, std::vector<int> ranks);
  /*!
     \brief
    Ranks of operator
     \return std::string
    */

  std::vector<int> ranks;
  /*!
     \brief
     hamiltonian object that defines the second-quantised spaces for the operator
     \return std::string
    */
  OldOperator hamiltonian;
};
}

#endif // GCIEXCITATIONOPERATOR_H
