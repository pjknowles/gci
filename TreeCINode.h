#ifndef TREECINODE_H
#define TREECINODE_H
#include <vector>
#include "TreeCIExcitation.h"

namespace TreeCI {
/*!
 \brief

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

};
}

using namespace TreeCI;

#endif // TREECINODE_H
