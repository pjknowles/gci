#ifndef TREECIEXCITATION_H
#define TREECIEXCITATION_H
namespace TreeCI {
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
    int h1,h2,p1,p2;
    int rank;

    static Excitation emptyTreeCIExcitation;
};
}

using namespace TreeCI;

#endif // TREECIEXCITATION_H
