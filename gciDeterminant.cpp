#include "gciDeterminant.h"
#include <iostream>

Determinant::Determinant(Parameters* parameters)
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

int Determinant::create(int orbital) {
//    xout << "orbitals.size()" << orbitals.size() << "orbital "<<orbital <<std::endl;
    unsigned int orbabs = orbital > 0 ? orbital : -orbital;
    if (orbital==(int)0 || orbital > (int) basisSize || orbital < -(int)basisSize) throw "invalid orbital";
    if (orbital > 0) {
        if (stringAlpha.orbitals().size() >= (nelec+ms2)/2) throw "too many electrons in determinant";
        return stringAlpha.create(orbabs);
    } else {
        if (stringBeta.orbitals().size() >= (nelec-ms2)/2) throw "too many electrons in determinant";
        return stringBeta.create(orbabs);
    }
}

int Determinant::destroy(int orbital) {
    if (orbital==(int)0 || orbital > (int) basisSize || orbital < -(int)basisSize) throw "invalid orbital";
    unsigned int orbabs = orbital > 0 ? orbital : -orbital;
    String* string = orbital > 0 ? &stringAlpha : &stringBeta;
    if (string->orbitals().size() <= 0) throw "too few electrons in determinant";
    return string->destroy(orbabs);

}
