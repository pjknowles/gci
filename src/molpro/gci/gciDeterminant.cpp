#include "gciDeterminant.h"
#include <memory>

namespace gci {

Determinant::Determinant(const State *State, const String *alpha, const String *beta) {
    if (State == nullptr) {
        nelec = 999999999;
        orbitalSpace = nullptr;
        ms2 = 0;
    } else {
        nelec = State->nelec;
        orbitalSpace = std::make_shared<OrbitalSpace>(*State->orbitalSpace);
        ms2 = State->ms2;
    }
    if (alpha != nullptr) stringAlpha = *alpha;
    if (beta != nullptr) stringBeta = *beta;
    if (State != nullptr) {
        stringAlpha.orbitalSpace = std::make_shared<OrbitalSpace>(*State->orbitalSpace);
        stringBeta.orbitalSpace = std::make_shared<OrbitalSpace>(*State->orbitalSpace);
    }
    //    xout << "determinant constructor, hamiltonian="<<(hamiltonian!=nullptr)<<hamiltonian->total()<<std::endl;
}

int Determinant::create(int orbital) {
    //    xout << "create orbital "<<orbital <<std::endl;
    auto orbabs = static_cast<String::orbital_type >(orbital > 0 ? orbital : -orbital);
    if (orbitalSpace == nullptr || orbital == (int) 0 || orbital > (int) orbitalSpace->total()
        || orbital < -(int) orbitalSpace->total())
        throw std::range_error("invalid orbital");
    if (orbital > 0) {
        if (stringAlpha.orbitals().size() >= (nelec + ms2) / 2)
            throw std::range_error("too many electrons in determinant");
        //        xout <<"try to populate stringAlpha"<<std::endl;
        return stringAlpha.create(orbabs);
    } else {
        if (stringBeta.orbitals().size() >= (nelec - ms2) / 2)
            throw std::range_error("too many electrons in determinant");
        return stringBeta.create(orbabs);
    }
}

int Determinant::destroy(int orbital) {
    if (orbitalSpace == nullptr || orbital == (int) 0 || orbital > (int) orbitalSpace->total()
        || orbital < -(int) orbitalSpace->total())
        throw std::range_error("invalid orbital");
    auto orbabs = static_cast<String::orbital_type>(orbital > 0 ? orbital : -orbital);
    String *string = orbital > 0 ? &stringAlpha : &stringBeta;
//  if (string->orbitals().size() <= 0) throw std::range_error("too few electrons in determinant");
    return string->destroy(orbabs);

}

//void Determinant::first()
//{
//  //    xout <<"Determiant::first nelec="<<nelec<<", ms2="<<ms2<<(nelec+ms2)/2<<(nelec-ms2)/2<<std::endl;
//  stringAlpha.first((nelec+ms2)/2);
//  stringBeta.first((nelec-ms2)/2);
//}

//bool Determinant::next()
//{
//  if (stringBeta.next()) return true;
//  //    xout << "Determinant::next needs to make a new alpha string"<<std::endl;
//  stringBeta.first((nelec-ms2)/2);
//  return stringAlpha.next();
//}

std::string Determinant::str(int verbosity, unsigned int columns) const {
    return verbosity >= 0 ? stringAlpha.str() + "|" + stringBeta.str() : std::string("");
}
}//  namespace gci
