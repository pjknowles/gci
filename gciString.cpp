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
    nelec=0;
    ms2=0;
    spin=Spin;
}

int String::create(unsigned int orbital) {
//        xout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
//        xout << "hamiltonian "<<(hamiltonian!=NULL)<<std::endl;
//        xout << "basisSize "<<hamiltonian->basisSize<<std::endl;
    if (hamiltonian==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) hamiltonian->basisSize) throw "invalid orbital";
//    xout <<"make iterator "<<std::endl;
    std::vector<unsigned int>::iterator ilast=orbitals_.begin();
//    xout <<"iterator OK"<<std::endl;

    int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
//    xout <<"phase="<<phase<<std::endl;
    ms2+=spin;
    nelec++;
    if (orbitals_.size() == 0 || *ilast > orbital) {
//        xout <<"first"<<std::endl;
        orbitals_.insert(orbitals_.begin(),orbital);
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

std::vector<unsigned int> String::orbitals() {
    return orbitals_;
}

std::string String::printable() {
    std::string result;
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (i!=orbitals_.begin()) result.append(",");
        std::stringstream ss;
        ss << *i * spin;
        std::string rr;
        ss >> rr;
        result.append(rr);
    }
    return result;
}

bool String::next() {
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

void String::first(int n) {
    if (n <=0 ) n=orbitals_.size();
    orbitals_.clear();
    for (int i=1;i<=n;i++)
        orbitals_.push_back(i);
}

String String::exhausted;

void String::buildStrings(State prototype, std::vector<String> *strings)
{
    String string(&prototype);
    string.first(prototype.nelec);
    strings->erase(strings->begin(),strings->end());
    do {
        strings->push_back(string);
    } while (string.next());
}
