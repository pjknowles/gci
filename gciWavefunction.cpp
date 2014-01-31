#include "gciWavefunction.h"
#include <sstream>
#include <iostream>

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) { }
Wavefunction::Wavefunction(std::string filename) : State(filename) { }

void Wavefunction::buildStrings()
{
    String stringa(this,1);
    stringa.first((nelec+ms2)/2);
    StringSet stringsa(stringa);
    alphaStrings = stringsa;
    for (StringSet::iterator s=stringsa.begin(); s!=stringsa.end(); s++) xout <<"Alpha string " << s->printable() <<std::endl;
//    xout << "before build beta" << std::endl;
    String stringb(this,-1);
    stringb.first((nelec-ms2)/2);
    xout << "Wavefunction::buildStrings nelec,ms2 "<<nelec<<ms2<<std::endl;
    xout << "Wavefunction::buildStrings beta prototype "<<stringb.printable()<<std::endl;
    StringSet stringsb(stringb);
    betaStrings = stringsb;
    for (StringSet::iterator s=stringsb.begin(); s!=stringsb.end(); s++) xout <<"beta string " << s->printable() <<std::endl;
//    xout << "end of Wavefunction::buildStrings" << std::endl;
}

std::string Wavefunction::printable(int verbosity)
{
    std::ostringstream s;
    if (verbosity >= 1) {
        s<<std::endl<< "Alpha strings:";
        for (StringSet::iterator i=alphaStrings.begin(); i!=alphaStrings.end(); i++) s <<std::endl<< i->printable();
        s<<std::endl<< "Beta strings:";
        for (StringSet::iterator i=betaStrings.begin(); i!=betaStrings.end(); i++) s <<std::endl<< i->printable();
    }
    return this->State::printable(verbosity)+s.str();
}
