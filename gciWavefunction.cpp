#include "gciWavefunction.h"
#include <sstream>
#include <iostream>

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) {
    buildStrings();
}

Wavefunction::Wavefunction(std::string filename) : State(filename) {
    if (filename!="") buildStrings();
}
Wavefunction::Wavefunction(Hamiltonian *h, int n, int s, int m2) : State(h,n,s,m2) {
   buildStrings();
}

void Wavefunction::buildStrings()
{
    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ symmetry;
        String stringa(this,1);
        stringa.first((nelec+ms2)/2);
        alphaStrings[syma] = StringSet(stringa, true, syma);
        String stringb(this,-1);
        stringb.first((nelec-ms2)/2);
        betaStrings[symb] = StringSet(stringb, true, symb);
    }
}

std::string Wavefunction::printable(int verbosity)
{
    std::ostringstream s;
    if (verbosity >= 1) {
        for (unsigned int syma=0; syma<8; syma++) {
            unsigned int symb = syma ^ symmetry;
            if (alphaStrings[syma].size() && betaStrings[symb].size()) {
                s<<std::endl<< "Alpha strings of symmetry "<<syma+1<<":";
                for (StringSet::iterator i=alphaStrings[syma].begin(); i!=alphaStrings[syma].end(); i++) s <<std::endl<< i->printable();
                s<<std::endl<< "Beta strings of symmetry "<<symb+1<<":";
                for (StringSet::iterator i=betaStrings[symb].begin(); i!=betaStrings[symb].end(); i++) s <<std::endl<< i->printable();
            }
        }
    }
    return this->State::printable(verbosity)+s.str();
}
