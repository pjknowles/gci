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
    dimension = 0;
    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ symmetry;
        String stringa(this,1);
        stringa.first((nelec+ms2)/2);
        alphaStrings[syma] = StringSet(stringa, true, syma);
        String stringb(this,-1);
        stringb.first((nelec-ms2)/2);
        betaStrings[symb] = StringSet(stringb, true, symb);
        dimension += alphaStrings[syma].size()*betaStrings[symb].size();
    }
}

size_t Wavefunction::size()
{
    return dimension;
}

void Wavefunction::allocate_buffer()
{
    buffer.resize(dimension,(double)0);
}


void Wavefunction::set(const double value)
{
    allocate_buffer();
    for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b=value;
}


Wavefunction& Wavefunction::operator*=(const double &value)
{
    for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b*=value;
    return *this;
}

Wavefunction& Wavefunction::operator*(const double &value)
{
//    return Wavefunction(*this) *= value;
    Wavefunction result = *this;
    result *= value;
    return result;
}

double Wavefunction::operator * ( const Wavefunction &ket)
{
    double result=(double)0;
    for (size_t i=0; i<buffer.size();
i++)
result += buffer[i]*ket.buffer[i];
    return result;
}

std::string Wavefunction::toString(int verbosity) const
{
    std::ostringstream s;
    if (verbosity >= 1) {
        size_t address=0;
        for (unsigned int syma=0; syma<8; syma++) {
            unsigned int symb = syma ^ symmetry;
            if (alphaStrings[syma].size() && betaStrings[symb].size()) {
                s<<std::endl<< "Alpha strings of symmetry "<<syma+1<<":";
                for (StringSet::const_iterator i=alphaStrings[syma].begin(); i!=alphaStrings[syma].end(); i++) s <<std::endl<< i->toString();
                s<<std::endl<< "Beta strings of symmetry "<<symb+1<<":";
                for (StringSet::const_iterator i=betaStrings[symb].begin(); i!=betaStrings[symb].end(); i++) s <<std::endl<< i->toString();
                if (buffer.size() == dimension && verbosity >=2) {
                        s<<std::endl<<"Coefficients:";
                    for (size_t i=0; i<alphaStrings[syma].size(); i++) {
                        s<<std::endl;
                        for (size_t j=0; j<betaStrings[symb].size(); j++) {
                            s << buffer[address++] << " ";
                        }
                    }
                }
            }
        }
    }
    return this->State::toString(verbosity)+s.str();
}
