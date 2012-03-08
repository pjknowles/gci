#include "TreeCINode.h"


Node::Node(Node* myparent,Excitation myexcitation)
{
    parent=myparent;
    excitation=myexcitation;
    terminal=true;
    if (parent!=NULL){
        parent->terminal=false;
        parent->children.push_back(this);
    }
}

//TreeCINode::~TreeCINode()
//{
////    if (parent!=NULL)
//}
