#ifndef GCIOPERATOR_H
#define GCIOPERATOR_H
#include "Operator.h"
#include "FCIdump.h"
#include "Eigen/Dense"


namespace gci {

  class Operator : SymmetryMatrix::Operator
  {
  public:
    using SymmetryMatrix::Operator::Operator;
    /*!
     * \brief Construct an object from an FCIdump. If the FCIdump
     * contains data, it will be loaded, otherwise the contents of the object are undefined,
     * and only the dimensions and parameters are loaded.
     * \param dump The raw buffer of a FCIdump.
     */
    static Operator construct(const FCIdump& dump);

    /*!
       * \brief int1 Generate array of diagonal one-electron integrals
       * \param spin positive for alpha, negative for beta
       * \return one-dimensional array with h(i,i) at i-1
       */
    Eigen::VectorXd int1(int spin) const;

    /*!
       * \brief intJ Generate array of two-electron exchange integrals
       * \param spini positive for alpha, negative for beta, first index
       * \param spinj positive for alpha, negative for beta, second index
       * \return array with (ii|jj)
       */
    Eigen::MatrixXd intJ(int spini, int spinj) const;
    /*!
       * \brief intK Generate array of two-electron Coulomb integrals
       * \param spin positive for alpha, negative for beta
       * \return array with (ij|ji)
       */
    Eigen::MatrixXd intK(int spin) const;

  private:
  std::vector<int> m_orbital_symmetries;

    /*!
     * \brief offset Return the number of orbitals of the same symmetry before the given one.
     * \param i Absolute number (starting with 1) of the orbital.
     * \return
     */
    size_t offset(unsigned int i) const {
      return i < 1 ? 0 : std::count(m_orbital_symmetries.begin(),m_orbital_symmetries.begin()+i-1,m_orbital_symmetries[i-1]);
    }
  };

}

#endif // GCIOPERATOR_H
