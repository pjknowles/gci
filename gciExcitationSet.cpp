#include "gciExcitationSet.h"

ExcitationSet::ExcitationSet(StringSet* from, StringSet* to, int annihilations, int creations)
{
    stringIndex.clear();
    phase.clear();
    orbitalAddress.clear();
    int symexc=-1;
    if (from->symmetry>=0 && to->symmetry>=0) symexc = from->symmetry ^ to->symmetry; // use symmetry if we can
    for (StringSet::iterator f=from->begin(); f!=from->end(); f++) {

    }
}
