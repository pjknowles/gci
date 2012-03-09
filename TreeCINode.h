#ifndef TREECINODE_H
#define TREECINODE_H
#include <vector>
#include "TreeCIExcitation.h"

using namespace TreeCI;

namespace TreeCI {
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
    Node(Node* parent=NULL,Excitation excitation=Excitation::emptyTreeCIExcitation);

    Node* parent; /*!< the parent of this node */
    bool terminal; /*!< whether or not this is a terminal node */
    std::vector<Node> children; /*!< the children of this node */
    Excitation excitation; /*!< the excitation that produced this node from its parent */
    std::vector<Excitation> excitations;
    Determinant Phi;

};
}


#endif // TREECINODE_H
