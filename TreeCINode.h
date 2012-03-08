#ifndef TREECINODE_H
#define TREECINODE_H
#include <vector>
#include "TreeCIExcitation.h"
/*!
 \brief

*/
class TreeCINode
{
public:
/*!
 \brief

*/
    TreeCINode(TreeCINode* parent=NULL,TreeCIExcitation excitation=TreeCIExcitation::emptyTreeCIExcitation);

    TreeCINode* parent; /*!< the parent of this node */
    bool terminal; /*!< whether or not this is a terminal node */
    std::vector<TreeCINode> children; /*!< the children of this node */
    TreeCIExcitation excitation; /*!< the excitation that produced this node from its parent */

};

#endif // TREECINODE_H
