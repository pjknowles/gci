#include "gciWavefunction.h"
#include <sstream>
#include <iostream>

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) {
    buildStrings();
}

Wavefunction::Wavefunction(std::string filename) : State(filename) {
    if (filename!="") buildStrings();
}

void Wavefunction::buildStrings()
{
    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ symmetry;
        String stringa(this,1);
        stringa.first((nelec+ms2)/2);
        StringSet stringsa(stringa, true, syma);
        alphaStrings[syma] = stringsa;
        //    for (StringSet::iterator s=stringsa.begin(); s!=stringsa.end(); s++) xout <<"Alpha string " << s->printable() <<std::endl;
        //    xout << "before build beta" << std::endl;
        String stringb(this,-1);
        stringb.first((nelec-ms2)/2);
        //    xout << "Wavefunction::buildStrings nelec,ms2 "<<nelec<<ms2<<std::endl;
        //    xout << "Wavefunction::buildStrings beta prototype "<<stringb.printable()<<std::endl;
        StringSet stringsb(stringb, true, symb);
        betaStrings[symb] = stringsb;
        //    for (StringSet::iterator s=stringsb.begin(); s!=stringsb.end(); s++) xout <<"beta string " << s->printable() <<std::endl;
        //    xout << "end of Wavefunction::buildStrings" << std::endl;
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
