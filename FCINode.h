#ifndef FCINODE_H
#define FCINODE_H
#include <vector>
#include "FCIExcitation.h"

using namespace FCI;

namespace FCI {
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
    Node(Node* parent=NULL,Excitation excitation=Excitation::emptyFCIExcitation);

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


#endif // FCINODE_H
