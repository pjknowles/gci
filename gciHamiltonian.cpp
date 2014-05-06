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
    load(dump,3);
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
    nt = SymmetryOffset("Numbers of orbitals in each symmetry");
    for (std::vector<int>::iterator s=syms.begin(); s!=syms.end(); s++) {
        orbital_symmetries[s-syms.begin()]=(*s)-1; // convert 1-8 to 0-7
        nt[(*s)-1]++;
    }
//    xout << "orbital_symmetries:"; for (int i=0; i<basisSize; i++) xout <<" "<<orbital_symmetries[i]; xout <<std::endl;
    nt.calculateOffsets(); // set up the rest of the nt object
//    xout << "orbital_symmetries:"; for (int i=0; i<basisSize; i++) xout <<" "<<orbital_symmetries[i]; xout <<std::endl;
    xout << nt.toString(2) << std::endl;
//    symmetric_pair_dimensions = SymmetryOffset("Numbers of orbital pairs in each symmetry");
//    for (std::vector<int>::iterator s=syms.begin(); s!=syms.end(); s++) {
//        for (std::vector<int>::iterator t=syms.begin(); t<=s; t++) {
//            unsigned int stsym = orbital_symmetries[s-syms.begin()]^orbital_symmetries[t-syms.begin()];
//            symmetric_pair_dimensions[stsym]++;
//        }
//    }
    symmetry_offsets_2e_ints = SymmetryOffset("2-electron integral symmetry offsets");
//    unsigned int off1=0;
//    unsigned int off2=0;
//    symmetry_offsets_pairs = SymmetryOffset("Orbital pair symmetry offsets");
//    for (unsigned int i=0; i<8; i++) {
//        symmetry_offsets_pairs[i]=off1;
//        symmetry_offsets_2e_ints[i]=off2;
//        off1+=nt[i]^2;
//        off2+=symmetric_pair_dimensions[i]^2;
////        xout <<"symmetry="<<i+1<<", symmetry_dimensions="<<symmetry_dimensions[i]<<std::endl;
////        xout <<"symmetry="<<i+1<<", symmetric_pair_dimensions="<<symmetric_pair_dimensions[i]<<std::endl;
////        xout <<"symmetry="<<i+1<<", symmetry_offsets_2e_ints="<<symmetry_offsets_2e_ints[i]<<std::endl;
//    }
    spinUnrestricted=false;

//    ijSize = (basisSize*(basisSize+1))/2;
//    ijklSize =ijSize*ijSize;
    ijSize = nt.total(0,1);
    ijklSize = 0;
    for (int isym=0; isym<8; isym++) {
        symmetry_offsets_2e_ints[isym]=ijklSize;
        ijklSize += nt.total(isym,1) * nt.total(isym,1);
    }
//    xout << symmetry_offsets_2e_ints << std::endl;
//    xout << "ijklSize=" << ijklSize <<std::endl;
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

std::string Hamiltonian::toString(int verbosity)
{
    std::ostringstream o;
    if (verbosity>=0) {
        o << "Basis size="<<basisSize<<" Spin unrestricted? "<<spinUnrestricted<<" Loaded? "<<loaded;
        o << std::endl << nt.toString();
        o << std::endl << symmetry_offsets_2e_ints.toString();
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
//                    unsigned int ij = i > j ? (i*(i-1))/2+j-1 : (j*(j-1))/2+i-1;
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
                            if (orbital_symmetries[i-1]^orbital_symmetries[j-1]^orbital_symmetries[k-1]^orbital_symmetries[l-1]) break;
                            unsigned int kl = k > l ? (k*(k-1))/2+l-1 : (l*(l-1))/2+k-1;
                            if (kl>ij) break;
//                            unsigned int ijkl = ij*ijSize+kl;
                    unsigned int ijkl = int2Index(i,j,k,l);
//                    xout << "i j k l =" <<i<<j<<k<<l<<" index="<<ijkl<<std::endl;
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
    unsigned int isym = orbital_symmetries[i-1];
    unsigned int jsym = orbital_symmetries[j-1];
    if (isym == jsym)
        return nt.offset(0,isym,1) + ((ii>jj) ? (ii*(ii+1))/2+jj : (jj*(jj+1))/2+ii);
    else if (isym > jsym ) {
//        xout << "in pairIndex(" << i << "," << j <<")="<< nt.offset(isym^jsym,isym,1) + jj*nt[isym] + ii <<"; offset=" << nt.offset(isym^jsym,isym,1) <<std::endl;
        return nt.offset(isym^jsym,isym,1) + jj*nt[isym] + ii;
    }
    else
        return nt.offset(isym^jsym,jsym,1) + ii*nt[jsym] + jj;
}

unsigned int Hamiltonian::int1Index(unsigned int i, unsigned int j) {
    return pairIndex(i,j);
}

unsigned int Hamiltonian::int2Index(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
{
    unsigned int ijsym = orbital_symmetries[i-1]^orbital_symmetries[j-1];
    unsigned int ij = pairIndex(i,j);
    unsigned int kl = pairIndex(k,l);
//    xout << "orbital_symmetries:"; for (int i=0; i<basisSize; i++) xout <<" "<<orbital_symmetries[i]; xout <<std::endl;
//    xout << "int2Index - symmetries";  xout<<orbital_symmetries[i-1]<<":"<<orbital_symmetries[j-1]<<":"<<ijsym<<":"<<std::endl;
//    xout << "int2Index - offset array "; for (int i=0; i<8; i++) xout <<symmetry_offsets_2e_ints[i]<<" "; xout<<std::endl;
//    xout << "int2Index "<<i<<" "<<j<<" "<<k<<" "<<l<<"orbital indices"<< orbitalIndex(i)<<orbitalIndex(j)<<orbitalIndex(k)<<orbitalIndex(l)<<" ij="<<ij<<" kl="<<kl<<" offset="<<symmetry_offsets_2e_ints[ijsym]<<" address" <<symmetry_offsets_2e_ints[ijsym]+ij*nt.total(ijsym)+kl<<std::endl;
    return symmetry_offsets_2e_ints[ijsym]+ij*nt.total(ijsym,1)+kl;
}
std::vector<double> Hamiltonian::int1(int spin)
{
    std::vector<double> result(basisSize*basisSize,(double)0);
    std::vector<double> * integrals = spin < 0 ? integrals_b : integrals_a;
    for (unsigned int j=0; j < basisSize; j++) {
        for (unsigned int i=0; i < basisSize; i++) {
            if (orbital_symmetries[i]==orbital_symmetries[j])
                result[i+j*basisSize] = integrals->at(int1Index(i+1,j+1));
        }
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
