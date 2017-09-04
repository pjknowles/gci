#ifndef GCIDETERMINANT_H
#define GCIDETERMINANT_H
#include "gci.h"
#include "gciState.h"
#include "gciString.h"
#include <vector>
namespace gci {
/*!
 \brief
A Slater determinant
*/
class Determinant : public State
{
public:
  /*!
 \brief

 \param State some object from which to copy number of electrons etc for bound checking
 \param alpha Alpha string
 \param beta Beta string
*/
  Determinant(State* State=nullptr, String* alpha=nullptr, String* beta=nullptr);

  /*!
     \brief

     \param orbital Add an orbital to the determinant. Negative means beta spin.
     \return int On exit, the phase change required to bring the determinant into canonical form is returned (plus or minus 1), or else zero if the orbital was already present in the determinant.
    */
  int create(int orbital);
  /*!
     \brief

     \param orbital Remove an orbital from the determinant. Negative means beta spin
     \return int On exit, return 1 if successful (orbital was in the determinant originally) or 0 if not (it wasn't)
    */
  int destroy(int orbital);
  /*!
     * \brief Set to the canonically first determinant
     */
  void first();
  /*!
     \brief
    Advance to the canonically next determinant
     \return false if the end of the set is reached.
    */
  bool next();

  /*!
     * \brief Genenerate a printable representation of the object
     * \param verbosity How much to print
     * \param columns Page width
     * \return
     */
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;


  String stringAlpha; /*!< The String of alpha-spin orbitals that make up the determinant */
  String stringBeta; /*!< The String of beta-spin orbitals that make up the determinant */

private:
  friend class Operator;
};
}

using namespace gci;

#endif // GCIDETERMINANT_H
