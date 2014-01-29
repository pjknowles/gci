#include "gciHamiltonian.h"
#include "FCIdump.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <ostream>
#include <iterator>

using namespace std;

Hamiltonian::Hamiltonian(string filename)
{
    loaded = false;
    verbosity=1;
    if (filename != "") load(filename);
}

Hamiltonian::Hamiltonian(FCIdump* dump)
{
    loaded = false;
    verbosity=1;
    load(dump);
}

Hamiltonian::~Hamiltonian() {

}

void Hamiltonian::load(string filename) {
    FCIdump d(filename);
    load(&d);
}

void Hamiltonian::load(FCIdump* dump) {
    if (loaded) unload();
    if (verbosity) cout <<"Load hamiltonian from " << dump->fileName <<endl;
//    State::load(filename);

    basisSize = dump->parameter("NORB").at(0);
    spinUnrestricted=false;

    ijSize = (basisSize*(basisSize+1))/2;
    ijklSize =ijSize*ijSize;
    integrals_a = new vector<double>(ijSize,0.0);
    if (spinUnrestricted)
        integrals_b = new vector<double>(ijSize,0.0);
    else
        integrals_b = integrals_a;
    integrals_aa = new vector<double>(ijklSize,0.0);
    if (spinUnrestricted) {
        integrals_ab = new vector<double>(ijklSize,0.0);
        integrals_bb = new vector<double>(ijklSize,0.0);
    } else {
        integrals_ab = integrals_aa;
        integrals_bb = integrals_aa;
    }


    ifstream s;
    s.open(dump->fileName.c_str());
    string ss;
    double value;
    int i,j,k,l,ij,kl,ijkl;
    while (s >> ss && ss != "&END") ;
    while (s >> value) {
        s >> i; s >> j; s >> k; s >> l;
        ij = i > j ? (i*(i-1))/2+j : (j*(j-1))/2+i;
        kl = k > l ? (k*(k-1))/2+l : (l*(l-1))/2+k;
        if (kl) {
         ijkl = (ij-1)*ijSize+kl-1;
        if (verbosity>2) cout << "("<< i << j <<"|"<< k << l <<") = " << value <<endl;
        integrals_aa->at(ijkl)=value;
        } else if (ij) {
            if (verbosity>1) cout << "h("<< i <<","<< j <<") = " << value <<endl;
            integrals_a->at(ij-1)=value;
        } else
            coreEnergy = value;
    }
    s.close();
    loaded=true;
    if (verbosity>3) {
        cout << "integrals_a: ";copy(integrals_a->begin(), integrals_a->end(), ostream_iterator<double>(cout, ", "));cout <<endl;
        cout << "integrals_aa: ";copy(integrals_aa->begin(), integrals_aa->end(), ostream_iterator<double>(cout, ", "));cout <<endl;
    }

    }

void Hamiltonian::unload() {
    if (loaded) {
    delete [] integrals_a;
    delete [] integrals_aa;
    if (spinUnrestricted) {
        delete [] integrals_b;
        delete [] integrals_ab;
        delete [] integrals_bb;
    }
    }
    loaded=false;
}
