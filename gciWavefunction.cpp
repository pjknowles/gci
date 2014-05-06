#include "gciWavefunction.h"
#include <sstream>
#include <iostream>
#include "gciStringSet.h"

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) {
    buildStrings();
}

Wavefunction::Wavefunction(std::string filename) : State(filename) {
    if (filename!="") buildStrings();
}
Wavefunction::Wavefunction(Hamiltonian *h, int n, int s, int m2) : State(h,n,s,m2) {
    buildStrings();
}

Wavefunction::Wavefunction(const Wavefunction &other) : State(other)
{
    for (int i=0;i<8;i++)
    {
        alphaStrings[i] = other.alphaStrings[i];
        betaStrings[i] = other.betaStrings[i];
    }
    dimension = other.dimension;
    buffer = other.buffer;
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

void Wavefunction::diagonalHamiltonian()
{
    std::vector<double> haa=hamiltonian->int1(1);
    std::vector<double> hbb=hamiltonian->int1(-1);
    std::vector<double> Jaa=hamiltonian->intJ(1,1);
    std::vector<double> Jab=hamiltonian->intJ(1,-1);
    std::vector<double> Jbb=hamiltonian->intJ(1,1);
    std::vector<double> Kaa=hamiltonian->intK(1);
    std::vector<double> Kbb=hamiltonian->intK(-1);
    xout << "Jbb" <<std::endl;
    for (int j=0; j<hamiltonian->basisSize; j++) {
        for (int i=0; i<hamiltonian->basisSize; i++)
            xout << Jbb[i+j*hamiltonian->basisSize] << " ";
        xout <<std::endl;
    }
    size_t offset=0;
    for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symb = syma ^ symmetry;
        std::vector<double> ona = alphaStrings[syma].occupationNumbers();
        std::vector<double> onb = betaStrings[symb].occupationNumbers();
        if (hamiltonian->spinUnrestricted) { // UHF
        } else { // RHF
            for (size_t ia=0; ia < alphaStrings[syma].size(); ia++) {
                std::vector<double> on(onb);
                for (size_t ib=0; ib < betaStrings[symb].size(); ib++)
                    on[ib] += ona[ia];
            }
        }
        offset += alphaStrings[syma].size()*betaStrings[symb].size();
    }
}

Wavefunction& Wavefunction::operator*=(const double &value)
{
    for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b*=value;
    return *this;
}

//Wavefunction& Wavefunction::operator = ( const Wavefunction &other)
//{
//    if (this != &other) {
//    if (! compatible(other)) throw "attempt to copy between incompatible Wavefunction objects";
//        *this = other;
//    }
//    return *this;
//}

Wavefunction& Wavefunction::operator+=(const Wavefunction &other)
{
    if (! compatible(other)) throw "attempt to add incompatible Wavefunction objects";
    xout << "Wavefunction::operator += &this=" << this <<" &other="<<&other <<std::endl;
    for (size_t i=0; i<buffer.size(); i++)  buffer[i] += other.buffer[i];
    return *this;
}


Wavefunction& Wavefunction::operator-=(const Wavefunction &other)
{
    if (! compatible(other)) throw "attempt to add incompatible Wavefunction objects";
    for (size_t i=0; i<buffer.size(); i++)  buffer[i] -= other.buffer[i];
    return *this;
}


Wavefunction gci::operator+(const Wavefunction &w1, const Wavefunction &w2)
{
    Wavefunction result = w1;
    return result += w2;
}

Wavefunction gci::operator-(const Wavefunction &w1, const Wavefunction &w2)
{
    Wavefunction result = w1;
    return result -= w2;
}

Wavefunction gci::operator*(const Wavefunction &w1, const double &value)
{
    Wavefunction result = w1;
    return result *= value;
}

Wavefunction gci::operator*(const double &value, const Wavefunction &w1)
{
    Wavefunction result = w1;
    return result *= value;
}

double gci::operator *(const Wavefunction &w1, const Wavefunction &w2)
{
    if (! w1.compatible(w2)) throw "attempt to form scalar product between incompatible Wavefunction objects";
    double result=(double)0;
    for (size_t i=0; i<w1.buffer.size(); i++) result += w1.buffer[i]*w2.buffer[i];
    return result;
}

std::string Wavefunction::str(int verbosity) const
{
    std::ostringstream s;
    if (verbosity >= 1) {
        s<<std::endl<<"Wavefunction object at address " << this ;
        s<<std::endl<<"Coefficients at address " << &buffer <<" and of length " << buffer.size();
        size_t address=0;
        for (unsigned int syma=0; syma<8; syma++) {
            unsigned int symb = syma ^ symmetry ;
            if (alphaStrings[syma].size() && betaStrings[symb].size()) {
                s<<std::endl<< "Alpha strings of  "<<syma+1<<":";
                for (StringSet::const_iterator i=alphaStrings[syma].begin(); i!=alphaStrings[syma].end(); i++) s <<std::endl<< i->str();
                s<<std::endl<< "Beta strings of symmetry "<<symb+1<<":";
                for (StringSet::const_iterator i=betaStrings[symb].begin(); i!=betaStrings[symb].end(); i++) s <<std::endl<< i->str();
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
    return this->State::str(verbosity)+s.str();
}

bool Wavefunction::compatible(const Wavefunction &other) const
{
    return dimension==other.dimension && buffer.size() == other.buffer.size();
}
