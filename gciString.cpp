#include "gciString.h"

#include <iostream>
#include <sstream>

String::String(State* State, int Spin)
{
    if (State == NULL) {
        hamiltonian=NULL;
    } else {
        hamiltonian=State->hamiltonian;
    }
    nullify();
    spin=Spin;
}

int String::create(unsigned int orbital) {
//        xout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
//        xout << "hamiltonian "<<(hamiltonian!=NULL)<<std::endl;
        if (hamiltonian==NULL)
            throw "String::create missing hamiltonian";
//        xout << "basisSize "<<hamiltonian->basisSize<<std::endl;
    if (hamiltonian==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) hamiltonian->basisSize) throw "invalid orbital";
//    xout <<"make iterator "<<std::endl;
    std::vector<unsigned int>::iterator ilast=orbitals_.begin();
//    xout <<"iterator OK"<<std::endl;

    int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
//    xout <<"phase="<<phase<<std::endl;
//    xout <<"spin="<<spin<<std::endl;
    ms2+=spin;
    nelec++;
    if (orbitals_.size() == 0 || *ilast > orbital) {
//        xout <<"first "<<orbital<<std::endl;
        orbitals_.insert(orbitals_.begin(),orbital);
//        xout <<"orbitals_[0]="<<orbitals_[0]<<std::endl;
        return phase;
    }
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital) return 0; // exclusion principle
        if (*ilast < orbital && *i > orbital){
            orbitals_.insert(++ilast,orbital);
            return phase;
        }
        phase=-phase;
        ilast=i;
    }
    orbitals_.insert(orbitals_.end(),orbital);
    return phase;
}

int String::destroy(unsigned int orbital) {
    if (hamiltonian==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) hamiltonian->basisSize ) throw "invalid orbital";
    if (orbitals_.size() <= 0) throw "too few electrons in String";
    int phase=1;
    ms2-=spin;
    nelec--;
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital)  {
            orbitals_.erase(i);
            return phase;
        }
        phase=-phase;
    }
    return (int)0; // exclusion principle
}

void String::nullify()
{
    orbitals_.clear();
    ms2=0;
    nelec=0;
}

std::vector<unsigned int> String::orbitals() {
    return orbitals_;
}

std::string String::printable(int verbosity) {
    std::string result;
//    xout <<"String::printable orbitals_[0]" <<orbitals_[0]<<std::endl;
    if (verbosity >=0) {
        for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
            if (i!=orbitals_.begin()) result.append(",");
            std::stringstream ss;
            int ispin=(int)(*i)*(int)spin;
            ss << ispin;
            std::string rr;
            ss >> rr;
            result.append(rr);
        }
        std::stringstream ss;
        ss << " ["<< symmetry()+1 <<"]"; // internally symmetries are implemented 0-7, but externally as 1-8
        std::string rr;
        ss >> rr;
        result.append(rr)
;    }
    return result;
}

unsigned int String::symmetry()
{
    unsigned int s=0;
    for (int i=0; i<(int)orbitals_.size(); i++) {
        s ^= hamiltonian->orbital_symmetries[orbitals_[i]-1];
//        xout <<"orbital symmetry="<<hamiltonian->orbital_symmetries[orbitals_[i]-1]<<" total symmetry now "<<s<<std::endl;
    }
    return s;
}

bool String::next(int sym) {
    if (sym<0)
    { // generate the next string of any symmetry
    unsigned int limit=hamiltonian->basisSize;
    std::vector<unsigned int>::iterator k;
    unsigned int floor;
    for (std::vector<unsigned int>::reverse_iterator i = orbitals_.rbegin(); i!=orbitals_.rend(); ++i) {
        floor=++(*i);
        k=i.base();
        if (*i <= limit) break;
        limit--;
    }
    if (limit <= hamiltonian->basisSize-orbitals_.size()) return false; // we ran out of boxes to put the objects into
    while (k!=orbitals_.rbegin().base())
        *(k++)=++floor;
    return true;
    }
    else { // call myself until we get the symmetry required
        bool notexhausted;
        while ((notexhausted=next())&&symmetry()!=(unsigned int)sym) ;
        return notexhausted;
    }
}

bool String::first(int n, int sym) {
//    xout << "String::first " << n <<" nelec="<<nelec<< std::endl;
    if (n <=0 ) n=nelec;
    if (n <=0 ) n=orbitals_.size();
    nullify();
//    xout << n <<std::endl;
    for (unsigned int i=1;i<=(unsigned int)n;i++)
        create(i);
    nelec=n;
//    xout << "String::first first go, sym, symmetry() "<<sym<<symmetry()<<std::endl;
    if (sym<0 || symmetry() == (unsigned int) sym) return true;
    return next(sym);
}

String String::exhausted;


