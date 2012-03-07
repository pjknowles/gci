#include "TreeCIHamiltonian.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>

using namespace std;

TreeCIHamiltonian::TreeCIHamiltonian()
{
}

TreeCIHamiltonian::~TreeCIHamiltonian() {

}

void TreeCIHamiltonian::loadHamiltonian(string filename) {
    cout <<"Load hamiltonian from " << filename <<endl;
    loadParameters(filename);
    spinUnrestricted=false;

    ijSize = (basisSize*(basisSize+1))/2;
    ijklSize =ijSize*ijSize;
    integrals_a = new double[ijSize];
    if (spinUnrestricted)
        integrals_b = new double[ijSize];
    else
        integrals_b = integrals_a;
    integrals_aa = new double[ijklSize];
    if (spinUnrestricted) {
        integrals_ab = new double[ijklSize];
        integrals_bb = new double[ijklSize];
    } else {
        integrals_ab = integrals_aa;
        integrals_bb = integrals_aa;
    }


    ifstream s;
    s.open(filename.c_str());
    string ss;
    double value;
    int i,j,k,l,ij,kl,ijkl;
    while (s >> ss && ss != "&END") ;
    while (s >> value) {
        s >> i;
        s >> j;
        s >> k;
        s >> l;
        ij = i > j ? (i*(i-1))/2+j : (j*(j-1))/2+i;
        kl = k > l ? (k*(k-1))/2+l : (l*(l-1))/2+k;
        if (kl) {
         ijkl = (ij-1)*ijSize+kl-1;
        cout << "("<< i << j <<"|"<< k << l <<") = " << value <<endl;
        integrals_aa[ijkl]=value;
        } else if (ij) {
            cout << "h("<< i <<","<< j <<") = " << value <<endl;
            integrals_a[ij-1]=value;
        } else
            coreEnergy = value;
    }
    s.close();
    loaded=true;
}

void TreeCIHamiltonian::unloadHamiltonian() {
    delete [] integrals_a;
    delete [] integrals_aa;
    if (spinUnrestricted) {
        delete [] integrals_b;
        delete [] integrals_ab;
        delete [] integrals_bb;
    }
    loaded=false;
}
