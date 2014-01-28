#ifndef gciNODE_H
#define gciNODE_H
#include <vector>
#include "gciExcitation.h"

using namespace gci;

namespace gci {
/*!
 \brief
A node in the tree representation of a wavefunction
*/
class Node
{
public:
/*!
 \brief

*/
    Node(Node* parent=NULL,Excitation excitation=Excitation::emptygciExcitation);

    /*!
     \brief
    Have children.
    */
    void procreate();

    Node* parent; /*!< the parent of this node */
    bool terminal; /*!< whether or not this is a terminal node */
    std::vector<Node> children; /*!< the children of this node */
    Excitation excitation; /*!< the excitation that produced this node from its parent */
    std::vector<Excitation> excitations;
    Determinant Phi;

};
}


#endif // gciNODE_H
