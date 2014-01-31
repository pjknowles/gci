#ifndef GCISTRING_H
#define GCISTRING_H
#include "gci.h"

#include "gciState.h"
#include <vector>

namespace gci {
/*!
 \brief
A string, which is an ordered set of orbitals
*/
class String : public State
{
public:
/*!
 \brief

 \param State Some State object from which to copy number of electrons etc for bound checking
 \param spin 1=alpha, -1=beta
*/
    String(State* State=NULL, int spin=1);
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
     * \brief Set the string to the vacuum
     */
    void nullify();    /*!
     \brief
    Set to the first string with n objects
    \param n the number of objects. If 0, use whatever we have presently.
    */
    void first(int n=0);
    std::vector<unsigned int> orbitals();  /*!< The orbitals that make up the string */
    /*!
     \brief
    printable form of the String.
    \param verbosity How much to print
     \return std::string
    */
    std::string printable(int verbosity=0);
    int spin; ///< \brief spin 1=alpha, -1=beta
    static String exhausted; /*!< returned by next() when we're already on the last string */

    /*!
     * \brief Build the complete set of Strings
     * \param prototype An object holding numbers of electrons and hamiltonian
     * \param strings Container to hold the complete set of Strings
     */
    static void buildStrings (State prototype, std::vector<String>* strings);

private:
    std::vector<unsigned int> orbitals_; /*!< The orbitals that make up the string */
    static vector< vector<int> > PartialWeightArray (int nitem,int nbox);
};
    long binomial_coefficient(unsigned long n, unsigned long k) ;
}

using namespace gci;

#endif // GCISTRING_H
