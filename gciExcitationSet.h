#ifndef GCIEXCITATIONSET_H
#define GCIEXCITATIONSET_H
#include "gciStringSet.h"
#include <vector>

/*!
 * \brief The ExcitationSet class
 */
class ExcitationSet
{
public:
    ExcitationSet(StringSet* from, StringSet* to, int annihilations, int creations);
    std::vector<int> stringIndex;
    std::vector<int> phase;
    std::vector<long> orbitalAddress;
};

#endif // GCIEXCITATIONSET_H
