#include "gciWavefunction.h"
#include <sstream>
#include <iostream>

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) { }
Wavefunction::Wavefunction(std::string filename) : State(filename) { }

void Wavefunction::buildStrings()
{
    String stringa(this,1);
    stringa.first((nelec+ms2)/2);
//    xout <<"buildStrings string"<<stringa.printable(1)<<std::endl;
//    xout <<"buildStrings string"<<stringa.State::printable(1)<<std::endl;
//    xout << "before build alpha" << std::endl;
    String::buildStrings(stringa,&alphaStrings);
//    for (std::vector<String>::iterator s=alphaStrings.begin(); s!=alphaStrings.end(); s++) xout <<"Alpha string " << (*s).printable() <<std::endl;
//    xout << "before build beta" << std::endl;
    String stringb(this,-1);
    stringb.first((nelec-ms2)/2);
    String::buildStrings(stringb,&betaStrings);
//    xout << "end of Wavefunction::buildStrings" << std::endl;
}

std::string Wavefunction::printable(int verbosity)
{
    std::ostringstream s;
    if (verbosity >= 1) {
        s<<std::endl<< "Alpha strings:";
        for (std::vector<String>::iterator i=alphaStrings.begin(); i!=alphaStrings.end(); s <<std::endl<< (*(i++)).printable() );
        s<<std::endl<< "Beta strings:";
        for (std::vector<String>::iterator i=betaStrings.begin(); i!=betaStrings.end(); s <<std::endl<< (*(i++)).printable() );
    }
    return this->State::printable(verbosity)+s.str();
}
