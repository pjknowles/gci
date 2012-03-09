#ifndef TREECISTRING_H
#define TREECISTRING_H

#include "TreeCIParameters.h"
#include <vector>

namespace TreeCI {
/*!
 \brief
A string, which is an ordered set of orbitals
*/
class String : public Parameters
{
public:
/*!
 \brief

 \param parameters some object from which to copy number of electrons etc for bound checking
*/
    String(Parameters* parameters=NULL);
    /*!
     \brief

     \param orbital Add an orbital to the string.
     \return int On exit, the phase change required to bring the determinant into canonical form is returned (plus or minus 1), or else zero if the orbital was already present in the determinant.
    */
    int create(int orbital);
    /*!
     \brief

     \param orbital Remove an orbital from the string.
     \return int On exit, return 1 if successful (orbital was in the determinant originally) or 0 if not (it wasn't)
    */
    int destroy(int orbital);
    /*!
     \brief
    get the canonically next string
     \return String
    */
    String next();
    std::vector<int> orbitals();  /*!< The orbitals that make up the string */
    /*!
     \brief
    printable form of the String.
     \return std::string
    */
    std::string printable();

private:
    std::vector<int> orbitals_; /*!< The orbitals that make up the string */
};
}

using namespace TreeCI;

#endif // TREECISTRING_H
