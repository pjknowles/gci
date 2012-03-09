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

int String::create(int orbital) {
    //    std::cout << "orbitals_.size()" << orbitals_.size() << "orbital "<<orbital <<std::endl;
    if (orbital==(int)0 || orbital > (int) basisSize) throw "invalid orbital";
    std::vector<int>::iterator ilast=orbitals_.begin();

    int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
    if (orbitals_.size() == 0 || *ilast > orbital) {
        orbitals_.insert(orbitals_.begin(),orbital);
        return phase;
    }
    for (std::vector<int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
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

int String::destroy(int orbital) {
    if (orbital==(int)0 || orbital > (int) basisSize || orbital < -(int)basisSize) throw "invalid orbital";
    if (orbitals_.size() <= 0) throw "too few electrons in String";
    int phase=1;
    for (std::vector<int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital)  {
            orbitals_.erase(i);
            return phase;
        }
        phase=-phase;
    }
    return 0; // exclusion principle
}

std::vector<int> String::orbitals() {
    return orbitals_;
}

std::string String::printable() {
    std::string result;
    for (std::vector<int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (i!=orbitals_.begin()) result.append(",");
        std::stringstream ss;
        ss << *i;
        std::string rr;
        ss >> rr;
        result.append(rr);
    }
    return result;
}
