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
    int create(unsigned int orbital);
    /*!
     \brief

     \param orbital Remove an orbital from the string.
     \return int the phase change required to bring the determinant into canonical form before annihilation is returned (plus or minus 1), or else zero if the orbital was not present in the determinant.
    */
    int destroy(unsigned int orbital);
    /*!
     \brief
    advance to the canonically next string
     \return whether successful; false if you try to advance the canonically last string
    */
    bool next();
    /*!
     \brief
    Set to the first string with n objects
    \param n the number of objects. If 0, use whatever we have presently.
    */
    void first(int n=0);
    std::vector<unsigned int> orbitals();  /*!< The orbitals that make up the string */
    /*!
     \brief
    printable form of the String.
     \return std::string
    */
    std::string printable();
    static String exhausted; /*!< returned by next() when we're already on the last string */

private:
    std::vector<unsigned int> orbitals_; /*!< The orbitals that make up the string */
};
}

using namespace TreeCI;

#endif // TREECISTRING_H
