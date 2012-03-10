#include "TreeCIString.h"

#include <iostream>
#include <sstream>

String::String(Parameters* parameters)
{
    if (parameters == NULL) {
        basisSize=99999999;
    } else {
        basisSize=parameters->basisSize;
    }
}

int String::create(unsigned int orbital) {
    //    std::cout << "orbitals_.size()" << orbitals_.size() << "orbital "<<orbital <<std::endl;
    if (orbital==(unsigned int)0 || orbital > (unsigned int) basisSize) throw "invalid orbital";
    std::vector<unsigned int>::iterator ilast=orbitals_.begin();

    int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
    if (orbitals_.size() == 0 || *ilast > orbital) {
        orbitals_.insert(orbitals_.begin(),orbital);
        return phase;
    }
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital) return 0; // exclusion principle
        if (*ilast < orbital && *i > orbital){
            orbitals_.insert(ilast,orbital);
            return phase;
        }
        phase=-phase;
        ilast=i;
    }
    orbitals_.insert(orbitals_.end(),orbital);
    return phase;
}

int String::destroy(unsigned int orbital) {
    if (orbital==(unsigned int)0 || orbital > (unsigned int) basisSize ) throw "invalid orbital";
    if (orbitals_.size() <= 0) throw "too few electrons in String";
    int phase=1;
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
        ss << *i;
        std::string rr;
        ss >> rr;
        result.append(rr);
    }
    return result;
}

int String::next() {
    unsigned int limit=basisSize;
    std::vector<unsigned int>::reverse_iterator j;
    for (std::vector<unsigned int>::reverse_iterator i = orbitals_.rbegin(); i!=orbitals_.rend(); ++i) {
        j=i;
        if (++(*i) <= limit) break;
        limit--;
    }
    if (limit < basisSize-orbitals_.size()) return 0;
    limit=*j;
    for (--j; j!=orbitals_.rend();++j)
        *(j)=(--limit);
    return 1;
}

String String::exhausted;
