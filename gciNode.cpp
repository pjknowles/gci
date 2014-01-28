#include "gciNode.h"


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

//gciNode::~gciNode()
//{
////    if (parent!=NULL)
//}

void Node::procreate() {
    std::vector<unsigned int> notHole, notParticle;
    for (Node* p=parent; p!=NULL; p=p->parent) { // do not repeat excitations made by parents
        if (p->excitation.h1) notHole.push_back(p->excitation.h1);
        if (p->excitation.h2) notHole.push_back(p->excitation.h2);
        if (p->excitation.p1) notParticle.push_back(p->excitation.p1);
        if (p->excitation.p2) notParticle.push_back(p->excitation.p2);

        for (unsigned int h1=1;;) ;
}
}
