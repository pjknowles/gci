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

Hamiltonian::Hamiltonian(std::string filename) : OrbitalSpace(filename)
{
    loaded = false;
    if (filename != "") load(filename);
}

Hamiltonian::Hamiltonian(FCIdump* dump) : OrbitalSpace(dump)
{
    loaded = false;
    load(dump,0);
}

Hamiltonian::~Hamiltonian() {

}

void Hamiltonian::load(std::string filename, int verbosity) {
    FCIdump d(filename);
    load(&d, verbosity);
}

void Hamiltonian::load(FCIdump* dump, int verbosity) {
    if (loaded) unload();
    if (verbosity) xout <<"Load hamiltonian from " << dump->fileName <<std::endl;
    //    State::load(filename);

    basisSize = dump->parameter("NORB").at(0);

    ijSize = total(0,1);
    ijklSize = pairSpace[1].total(0);
    xout << "ijklSize=" << ijklSize <<std::endl;
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
            if (verbosity>2) xout << "("<< i << j <<"|"<< k << l <<") [" << int2Index(i,j,k,l) << "]= " << value <<std::endl;
            if (verbosity>2) xout << "("<< k << l <<"|"<< i << j <<") [" << int2Index(k,l,i,j) << "]= " << value <<std::endl;
            integrals_aa->at(int2Index(i,j,k,l))=value;
            integrals_aa->at(int2Index(k,l,i,j))=value;
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

    // construct <ik||jl> = (ij|kl) - (il|kj)
    for (unsigned int symik=0; symik<8; symik++) {

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

std::string Hamiltonian::str(int verbosity) const
{
    std::ostringstream o;
    o << OrbitalSpace::str(verbosity);
    if (verbosity>=2) {
        int precision=6;
        o<<std::setprecision(precision);
        o << std::endl << "Core energy " << coreEnergy;
        if (integrals_a!=NULL && integrals_b != NULL) {
            o<<std::endl << "1-electron integrals:";
            for (unsigned int i=1; i<=basisSize; i++) {
                for (unsigned int j=1; j<=i; j++) {
                    if (orbital_symmetries[i-1]^orbital_symmetries[j-1]) continue;
                    unsigned int ij = int1Index(i,j);
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
                            if (orbital_symmetries[i-1]^orbital_symmetries[j-1]^orbital_symmetries[k-1]^orbital_symmetries[l-1]) continue;
                            unsigned int kl = k > l ? (k*(k-1))/2+l-1 : (l*(l-1))/2+k-1;
                            if (kl>ij) break;
                            unsigned int ijkl = int2Index(i,j,k,l);
                            if (integrals_aa->at(ijkl) != (double)0 || integrals_ab->at(ijkl) != (double)0 || integrals_bb->at(ijkl) != (double)0) o<<std::endl<<std::setw(4)<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<std::setw(precision+7)<<integrals_aa->at(ijkl)<<" "<<integrals_ab->at(ijkl)<<" "<<integrals_bb->at(ijkl);
                        }
                    }
                }
            }
        }
    }
    return o.str();
}

unsigned int Hamiltonian::int1Index(unsigned int i, unsigned int j) const {
    return pairIndex(i,j,1);
}

unsigned int Hamiltonian::int2Index(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
{
    return quadIndex(i,j,k,l,1,0);
}

std::vector<double> Hamiltonian::int1(int spin)
{
    std::vector<double> result(basisSize,(double)0);
    std::vector<double> * integrals = spin < 0 ? integrals_b : integrals_a;
    for (unsigned int i=0; i < basisSize; i++) {
        result[i] = integrals->at(int1Index(i+1,i+1));
    }
    return result;
}


std::vector<double> Hamiltonian::intJ(int spini, int spinj)
{
    std::vector<double> result(basisSize*basisSize,(double)0);
    std::vector<double> * integrals = spini < 0 ? (spinj < 0 ? integrals_bb : integrals_ab) : integrals_aa;
    for (unsigned int j=0; j < basisSize; j++) {
        for (unsigned int i=0; i < basisSize; i++) {
            result[i+j*basisSize] = integrals->at(int2Index(i+1,i+1,j+1,j+1));
            //            xout <<"intJ["<<i<<","<<j<<"]="<<result[i+j*basisSize]<<std::endl;
        }
    }
    return result;
}

std::vector<double> Hamiltonian::intK(int spin)
{
    std::vector<double> result(basisSize*basisSize,(double)0);
    std::vector<double> * integrals = spin < 0 ? integrals_bb : integrals_aa;
    for (unsigned int j=0; j < basisSize; j++) {
        for (unsigned int i=0; i < basisSize; i++) {
            result[i+j*basisSize] = integrals->at(int2Index(i+1,j+1,j+1,i+1));
        }
    }
    return result;
}

Hamiltonian Hamiltonian::FockHamiltonian(Determinant &reference)
{
    Hamiltonian f;
    for (int i=0; i<8; i++)
        f[i]=at(i);
    f.calculateOffsets();
    f.spinUnrestricted = spinUnrestricted;
    f.coreEnergy = coreEnergy;
    f.basisSize = basisSize;
    f.ijklSize = ijklSize;
    f.ijSize = ijSize;
    f.orbital_symmetries = orbital_symmetries;
    f.integrals_a = integrals_a;
    for (std::vector<unsigned int>::const_iterator o=reference.stringAlpha.orbitals().begin(); o != reference.stringAlpha.orbitals().end(); o++)
    {
        for (unsigned int i=1; i<=basisSize; i++)
            for (unsigned int j=1; j<=i; j++) {
                if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
                (*f.integrals_a)[int1Index(i,j)] += (*integrals_aa)[int2Index(i,j,*o,*o)] - (*integrals_aa)[int2Index(i,*o,*o,j)];
            }
    }
    for (std::vector<unsigned int>::const_iterator o=reference.stringBeta.orbitals().begin(); o != reference.stringBeta.orbitals().end(); o++)
    {
        for (unsigned int i=1; i<=basisSize; i++)
            for (unsigned int j=1; j<=i; j++) {
                if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
                (*f.integrals_a)[int1Index(i,j)] += (*integrals_ab)[int2Index(i,j,*o,*o)];
            }
    }
    f.integrals_aa = NULL;
    f.integrals_ab = NULL;
    f.integrals_bb = NULL;
    if (spinUnrestricted) {
        f.integrals_b = integrals_b;
        for (std::vector<unsigned int>::const_iterator o=reference.stringBeta.orbitals().begin(); o != reference.stringBeta.orbitals().end(); o++)
        {
            for (unsigned int i=1; i<=basisSize; i++)
                for (unsigned int j=1; j<=i; j++) {
                    if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
                    (*f.integrals_a)[int1Index(i,j)] += (*integrals_bb)[int2Index(i,j,*o,*o)] - (*integrals_bb)[int2Index(i,*o,*o,j)];
                }
        }
        for (std::vector<unsigned int>::const_iterator o=reference.stringAlpha.orbitals().begin(); o != reference.stringAlpha.orbitals().end(); o++)
        {
            for (unsigned int i=1; i<=basisSize; i++)
                for (unsigned int j=1; j<=i; j++) {
                    if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
                    (*f.integrals_a)[int1Index(i,j)] += (*integrals_ab)[int2Index(*o,*o,i,j)];
                }
        }
    } else {
        f.integrals_b = f.integrals_a;
    }
    f.loaded = true;
    return f;
}
