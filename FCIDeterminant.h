#ifndef FCIDETERMINANT_H
#define FCIDETERMINANT_H
#include "FCIParameters.h"
#include "FCIString.h"
#include <vector>

namespace FCI {
/*!
 \brief
A Slater determinant
*/
class Determinant : public Parameters
{
public:
/*!
 \brief

 \param parameters some object from which to copy number of electrons etc for bound checking
*/
    Determinant(Parameters* parameters=NULL);
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
     \brief
    get the canonically next determinant
     \return Determinant
    */
    Determinant next();

private:
    String stringAlpha, stringBeta; /*!< The orbitals that make up the determinant */
};
}

using namespace FCI;

#endif // FCIDETERMINANT_H
