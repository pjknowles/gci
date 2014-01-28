#ifndef FCIEXCITATION_H
#define FCIEXCITATION_H
#include "FCIDeterminant.h"
namespace FCI {
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

    static Excitation emptyFCIExcitation;
};
}

using namespace FCI;

#endif // FCIEXCITATION_H
