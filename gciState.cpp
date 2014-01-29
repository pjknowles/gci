#include "gciState.h"
#include "FCIdump.h"
using namespace std;
#include <fstream>
#include <iostream>
#include <sstream>

State::State(string filename)
{
    hamiltonian=NULL;
    if (filename!="") load(filename);
}

State::State(FCIdump* dump)
{
    hamiltonian=NULL;
    load(dump);
}

State::~State()
{
}


void State::load(string filename) {
    FCIdump dump(filename);
    load(&dump);
}

void State::load(FCIdump* dump)
{
    nelec = dump->parameter("NELEC").at(0);
    ms2 = dump->parameter("MS2").at(0);
    symmetry = dump->parameter("ISYM",std::vector<int>(1,1)).at(0);
    cout <<"nelec="<<nelec<<endl;
    cout <<"ms2="<<ms2<<endl;
}
