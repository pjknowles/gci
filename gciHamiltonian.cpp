#include "gciHamiltonian.h"
#include "FCIdump.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>

Hamiltonian::Hamiltonian(std::string filename)
{
    loaded = false;
    if (filename != "") load(filename);
}

Hamiltonian::Hamiltonian(FCIdump* dump)
{
    loaded = false;
    load(dump);
}

Hamiltonian::~Hamiltonian() {

}

void Hamiltonian::load(std::string filename, int verbosity) {
    FCIdump d(filename);
    load(&d, verbosity);
}

void Hamiltonian::load(FCIdump* dump, int verbosity) {
    int debug=1;
    if (loaded) unload();
    if (verbosity) xout <<"Load hamiltonian from " << dump->fileName <<std::endl;
    //    State::load(filename);

    basisSize = dump->parameter("NORB").at(0);
    std::vector<int> syms = dump->parameter("ORBSYM");
    orbital_symmetries = std::vector<unsigned int>(basisSize,0);
//    orbital_symmetries = dump->parameter("ORBSYM");;
    symmetry_dimensions = SymmetryOffset("Numbers of orbitals in each symmetry");
    symmetric_pair_dimensions = SymmetryOffset("Numbers of orbital pairs in each symmetry");
    for (std::vector<int>::iterator s=syms.begin(); s!=syms.end(); s++) {
        orbital_symmetries[s-syms.begin()]=(*s)-1; // convert 1-8 to 0-7
        symmetry_dimensions[(*s)-1]++;
    }
    for (std::vector<int>::iterator s=syms.begin(); s!=syms.end(); s++) {
        for (std::vector<int>::iterator t=syms.begin(); t<=s; t++) {
            unsigned int stsym = orbital_symmetries[s-syms.begin()]^orbital_symmetries[t-syms.begin()];
            symmetric_pair_dimensions[stsym]++;
        }
    }
    unsigned int off1=0;
    unsigned int off2=0;
    symmetry_offsets_pairs = SymmetryOffset("Orbital pair symmetry offsets");
    symmetry_offsets_2e_ints = SymmetryOffset("2-electron integral symmetry offsets");
    for (unsigned int i=0; i<8; i++) {
        symmetry_offsets_pairs[i]=off1;
        symmetry_offsets_2e_ints[i]=off2;
        off1+=symmetry_dimensions[i]^2;
        off2+=symmetric_pair_dimensions[i]^2;
//        xout <<"symmetry="<<i+1<<", symmetry_dimensions="<<symmetry_dimensions[i]<<std::endl;
//        xout <<"symmetry="<<i+1<<", symmetric_pair_dimensions="<<symmetric_pair_dimensions[i]<<std::endl;
//        xout <<"symmetry="<<i+1<<", symmetry_offsets_2e_ints="<<symmetry_offsets_2e_ints[i]<<std::endl;
    }
    spinUnrestricted=false;

    if (debug) {//debugging
        for (unsigned int i=0; i<basisSize; i++) {
            xout <<"Orbital "<<i+1<<" symmetry="<<orbital_symmetries[i]+1<<", index="<<orbitalIndex(i+1)<<std::endl;
            if (debug>1)
                for (unsigned int j=0; j<basisSize; j++) {
                    xout <<"Pair "<<i+1<<","<<j+1<<" symmetry="<<(orbital_symmetries[i]^orbital_symmetries[j])+1<<", index="<<pairIndex(i+1,j+1)<<std::endl;
                }
        }
    }

    ijSize = (basisSize*(basisSize+1))/2;
    ijklSize =ijSize*ijSize;
    integrals_a = new std::vector<double>(ijSize,0.0);
    if (spinUnrestricted)
        integrals_b = new std::vector<double>(ijSize,0.0);
    else
        integrals_b = integrals_a;
    integrals_aa = new std::vector<double>(ijklSize,0.0);
    if (spinUnrestricted) {
        integrals_ab = new std::vector<double>(ijklSize,0.0);
        integrals_bb = new std::vector<double>(ijklSize,0.0);
    } else {
        integrals_ab = integrals_aa;
        integrals_bb = integrals_aa;
    }


    std::ifstream s;
    s.open(dump->fileName.c_str());
    std::string ss;
    double value;
    int i,j,k,l,ij,kl,ijkl;
    while (s >> ss && ss != "&END") ;
    while (s >> value) {
        s >> i; s >> j; s >> k; s >> l;
        ij = i > j ? (i*(i-1))/2+j : (j*(j-1))/2+i;
        kl = k > l ? (k*(k-1))/2+l : (l*(l-1))/2+k;
        if (kl) {
            ijkl = (ij-1)*ijSize+kl-1;
            if (verbosity>2) xout << "("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
            integrals_aa->at(int2Index(i,j,k,l))=value;
        } else if (ij) {
            if (verbosity>1) xout << "h("<< i <<","<< j <<") = " << value <<std::endl;
            integrals_a->at(int1Index(i,j))=value;
        } else
            coreEnergy = value;
    }
    s.close();
    loaded=true;
    if (verbosity>3) {
        xout << "integrals_a: ";copy(integrals_a->begin(), integrals_a->end(), std::ostream_iterator<double>(xout, ", "));xout <<std::endl;
        xout << "integrals_aa: ";copy(integrals_aa->begin(), integrals_aa->end(), std::ostream_iterator<double>(xout, ", "));xout <<std::endl;
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

std::string Hamiltonian::printable(int verbosity)
{
    std::ostringstream o;
    if (verbosity>=0) {
        o << "Basis size="<<basisSize<<" Spin unrestricted? "<<spinUnrestricted<<" Loaded? "<<loaded;
        o << std::endl << symmetry_dimensions.printable();
        o << std::endl << symmetric_pair_dimensions.printable();
        o << std::endl << symmetry_offsets_pairs.printable();
        o << std::endl << symmetry_offsets_2e_ints.printable();
    }
    if (verbosity>=1) {
        o << std::endl << "Orbital symmetries";
        for (std::vector<unsigned int>::iterator i=orbital_symmetries.begin(); i!=orbital_symmetries.end(); i++)
            o<<" "<<*i+1;
    }
    if (verbosity>=2) {
        int precision=6;
        o<<std::setprecision(precision);
        o << std::endl << "Core energy " << coreEnergy;
        if (integrals_a!=NULL && integrals_b != NULL) {
            o<<std::endl << "1-electron integrals:";
            for (unsigned int i=1; i<=basisSize; i++) {
                for (unsigned int j=1; j<=i; j++) {
                    unsigned int ij = i > j ? (i*(i-1))/2+j-1 : (j*(j-1))/2+i-1;
                    if (integrals_a->at(ij) != (double)0 || integrals_b->at(ij) != (double)0) o<<std::endl<<std::setw(4)<<i<<" "<<j<<" "<<std::setw(precision+7)<<integrals_a->at(ij)<<" "<<integrals_b->at(ij);
                }
            }
        }
        if (integrals_aa!=NULL && integrals_ab != NULL && integrals_bb != NULL) {
            o<<std::endl << "2-electron integrals:";
            for (unsigned int i=1; i<=basisSize; i++) {
                for (unsigned int j=1; j<=i; j++) {
                    unsigned int ij = i > j ? (i*(i-1))/2+j-1 : (j*(j-1))/2+i-1;
                    for (unsigned int k=1; k<=i; k++) {
                        for (unsigned int l=1; l<=i; l++) {
                            unsigned int kl = k > l ? (k*(k-1))/2+l-1 : (l*(l-1))/2+k-1;
                            if (kl>ij) break;
                            unsigned int ijkl = ij*ijSize+kl;
                            if (integrals_aa->at(ijkl) != (double)0 || integrals_ab->at(ijkl) != (double)0 || integrals_bb->at(ijkl) != (double)0) o<<std::endl<<std::setw(4)<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<std::setw(precision+7)<<integrals_aa->at(ijkl)<<" "<<integrals_ab->at(ijkl)<<" "<<integrals_bb->at(ijkl);
                        }
                    }
                }
            }
        }
    }
    return o.str();
}

unsigned int Hamiltonian::orbitalIndex(unsigned int i) {
    unsigned int n=0;
    for (unsigned int j=1; j<i; j++)
        if (orbital_symmetries[j-1] == orbital_symmetries[i-1]) n++;
    return n;
}

unsigned int Hamiltonian::pairIndex(unsigned int i, unsigned int j) {
    unsigned int ii = orbitalIndex(i);
    unsigned int jj = orbitalIndex(j);
    unsigned int ijsym = orbital_symmetries[i-1] ^ orbital_symmetries[j-1];
    return symmetry_offsets_pairs[ijsym] + ((ii>jj) ? (ii*(ii+1))/2+jj : (jj*(jj+1))/2+ii);
}

unsigned int Hamiltonian::int1Index(unsigned int i, unsigned int j) {
    return pairIndex(i,j);
}

unsigned int Hamiltonian::int2Index(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    unsigned int ijsym = orbital_symmetries[i-1]^orbital_symmetries[j-1];
    unsigned int ij = pairIndex(i,j);
    unsigned int kl = pairIndex(k,l);
    xout << "int2Index "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<ij<<" "<<kl<<" "<<symmetry_offsets_2e_ints[ijsym]+ij*symmetric_pair_dimensions[ijsym]+kl<<std::endl;
    return symmetry_offsets_2e_ints[ijsym]+ij*symmetric_pair_dimensions[ijsym]+kl;
}
