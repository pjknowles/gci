#include "gciExcitationSet.h"
Excitation::Excitation(long StringIndex, int Phase, long OrbitalAddress)
{
    stringIndex = StringIndex;
    phase = Phase;
    orbitalAddress = OrbitalAddress;
}

ExcitationSet::ExcitationSet(String &from, StringSet &to, int annihilations, int creations)
{
    int symexc=-1;
    if (to.symmetry>=0) symexc = from.symmetry ^ to.symmetry; // use symmetry if we can
    if (annihilations==1 && creations==0) {
        for (int i=0; i<(int)from.orbitals_.size(); i++) {
            if (from.hamiltonian->orbital_symmetries[i]==(unsigned int)symexc || symexc==-1) {
                String tt = from;
                int phase = tt.destroy(i);
                long ti=to.offset(tt);
                push_back(Excitation(ti,phase,i));
            }
        }
    }
}
