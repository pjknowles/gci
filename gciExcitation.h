#ifndef gciEXCITATION_H
#define gciEXCITATION_H
#include "gciDeterminant.h"
namespace gci {
/*!
 \brief

*/
class Excitation
{
public:
/*!
 \brief

*/
/*!
 \brief

*/
    Excitation();
    int index; /*!< The canonical index of the excitation */
    int h1,h2,p1,p2; /*!< The orbitals destroyed and created */
    int rank; /*!< How many electrons are moved */

    static Excitation emptygciExcitation;
};
}

using namespace gci;

#endif // gciEXCITATION_H
