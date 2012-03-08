#include "TreeCIDeterminant.h"
#include <iostream>

TreeCIDeterminant::TreeCIDeterminant(TreeCIParameters* parameters)
{
    if (parameters == NULL) {
        nelec=999999999;
        basisSize=99999999;
        ms2=0;
    } else {
        nelec=parameters->nelec;
        basisSize=parameters->basisSize;
        ms2=parameters->ms2;
    }
}

int TreeCIDeterminant::create(int orbital) {
//    std::cout << "orbitals.size()" << orbitals.size() << "orbital "<<orbital <<std::endl;
    if (orbital==(int)0 || orbital > (int) basisSize || orbital < -(int)basisSize) throw "invalid orbital";
    int phase=1;
    if (orbital > 0) {
        if (orbitalsAlpha.size() >= (nelec+ms2)/2) throw "too many electrons in determinant";
        for (std::vector<int>::iterator i = orbitalsAlpha.begin(); i!=orbitalsAlpha.end(); ++i) {
            if (*i==orbital) return 0; // exclusion principle
            if (*i < orbital){
                orbitalsAlpha.insert(i,orbital);
                return phase;
            }
            phase=-phase;
        }
        orbitalsAlpha.insert(orbitalsAlpha.end(),orbital);
    } else {
        if (orbitalsBeta.size() >= (nelec-ms2)/2) throw "too many electrons in determinant";
        for (std::vector<int>::iterator i = orbitalsBeta.begin(); i!=orbitalsBeta.end(); ++i) {
            if (*i==-orbital) return 0; // exclusion principle
            if (*i < -orbital){
                orbitalsBeta.insert(i,-orbital);
                return phase;
            }
            phase=-phase;
        }
    orbitalsBeta.insert(orbitalsBeta.end(),-orbital);
    }
    return phase;
}

int TreeCIDeterminant::destroy(int orbital) {
    if (orbital==(int)0 || orbital > (int) basisSize || orbital < -(int)basisSize) throw "invalid orbital";
    if (orbital > 0) {
        if (orbitalsAlpha.size() <= 0) throw "too few electrons in determinant";
        for (std::vector<int>::iterator i = orbitalsAlpha.begin(); i!=orbitalsAlpha.end(); ++i) {
            if (*i==orbital)  {
                orbitalsAlpha.erase(i);
                return 1;
            }
        }
    } else {
        if (orbitalsBeta.size() <= 0) throw "too few electrons in determinant";
        for (std::vector<int>::iterator i = orbitalsBeta.begin(); i!=orbitalsBeta.end(); ++i) {
            if (*i==-orbital)  {
                orbitalsBeta.erase(i);
                return 1;
            }
        }
    }
    return 0; // exclusion principle
}
