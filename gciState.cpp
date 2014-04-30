#include "gciState.h"
#include "FCIdump.h"
#include <fstream>
#include <iostream>
#include <sstream>

State::State(std::string filename)
{
    hamiltonian=NULL;
    if (filename!="") load(filename);
}

State::State(FCIdump* dump)
{
    hamiltonian=NULL;
    load(dump);
}

State::State(Hamiltonian *h, int n, int s, int m2)
{
    hamiltonian=h;
    nelec = n;
    symmetry = s;
    ms2 = m2;
}

State::State(State *s)
{
    hamiltonian = s->hamiltonian;
    nelec = s->nelec;
    symmetry = s->symmetry;
    ms2 = s->ms2;
}

State::~State()
{
}


void State::load(std::string filename) {
    FCIdump dump(filename);
    load(&dump);
}

void State::load(FCIdump* dump)
{
    nelec = dump->parameter("NELEC").at(0);
    ms2 = dump->parameter("MS2").at(0);
    symmetry = dump->parameter("ISYM",std::vector<int>(1,1)).at(0)-1; // D2h symmetries are 0-7 internally, 1-8 externally
//    xout <<"nelec="<<nelec<<endl;
//    xout <<"ms2="<<ms2<<endl;
    hamiltonian = new Hamiltonian(dump);
//    xout << "State::load hamiltonian=" << (hamiltonian != NULL) << std::endl;
//    xout << "basisSize=" << hamiltonian->basisSize <<std::endl;
}

std::string State::toString(int verbosity)
{
//    xout << "State::printable hamiltonian=" << (hamiltonian != NULL) << verbosity << std::endl;
//    xout << "basisSize=" << hamiltonian->basisSize <<std::endl;
    std::ostringstream s;
    if (verbosity >= 0) {
        s<< "nelec="<<nelec<<" ms2="<<ms2<<" symmetry="<<symmetry+1;
        if (hamiltonian!=NULL)
            s << std::endl << "Hamiltonian: " << hamiltonian->toString(1);
    }
//    xout << "basisSize=" << hamiltonian->basisSize <<std::endl;
    return s.str();
}
