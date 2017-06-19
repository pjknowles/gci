#include "gciOldOperator.h"
#include "FCIdump.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <cmath>
#include <numeric>
static gci::Operator dummyOperator{dim_t{},0};

template <class T>
std::vector<double>* vecdup(const T& source) {
  std::vector<double>* result = new std::vector<double>(source.size());
  std::copy(source.begin(), source.end(),result->begin());
  return result;
}

OldOperator::OldOperator(const gci::Operator &source)
  : m_Operator(source)
{
  std::vector<int> syms;
  for (auto& s : source.orbital_symmetries()) syms.push_back(s+1);
//  xout << "OldOperator from Operator: syms ="; for (auto& s: syms) xout << " "<<s; xout <<std::endl;
  OrbitalSpace::load(syms,source.m_uhf);
  basisSize = std::accumulate(source.m_dimensions[0].begin(),source.m_dimensions[0].end(),0);
  ijSize = total(0,1);
  ijklSize = pairSpace[1].total(0);
  integrals_a = new std::vector<double>(ijSize,0.0);
  integrals_aa = new std::vector<double>(ijklSize,0.0);
  if (spinUnrestricted) {
      integrals_b = new std::vector<double>(ijSize,0.0);
      integrals_ab = new std::vector<double>(ijklSize,0.0);
      integrals_bb = new std::vector<double>(ijklSize,0.0);
    } else {
      integrals_b = integrals_a;
      integrals_ab = integrals_aa;
      integrals_bb = integrals_aa;
    }

  std::copy(source.O1(true).data()->begin(),source.O1(true).data()->end(),integrals_a->begin());
  std::copy(source.O1(false).data()->begin(),source.O1(false).data()->end(),integrals_b->begin());
  std::copy(source.O2(true,true).data()->begin(),source.O2(true,true).data()->end(),integrals_aa->begin());
  std::copy(source.O2(true,false).data()->begin(),source.O2(true,false).data()->end(),integrals_ab->begin());
  std::copy(source.O2(false,false).data()->begin(),source.O2(false,false).data()->end(),integrals_bb->begin());
  coreEnergy = source.m_O0;

  loaded=true;
  bracket_integrals_a=bracket_integrals_b=NULL;
  bracket_integrals_aa=bracket_integrals_ab=bracket_integrals_bb=NULL;
//  constructBraKet();
//  std::cout << "bracket_integrals_aa:\n"; for(auto& s : *bracket_integrals_aa) std::cout << " "<<s ; std::cout << std::endl;
  bracket_integrals_a = vecdup(*source.O1(true).data());
  bracket_integrals_aa = vecdup(*source.O2(true,true,false).data());
  bracket_integrals_ab = vecdup(*source.O2(true,false,false).data());
  if (spinUnrestricted) {
      bracket_integrals_b = vecdup(*source.O1(false).data());
      bracket_integrals_bb = vecdup(*source.O2(false,false,false).data());
    }
  else {
      bracket_integrals_bb = bracket_integrals_aa;
      bracket_integrals_b = bracket_integrals_a;
    }
//  std::cout << "bracket_integrals_aa:\n"; for(auto& s : *bracket_integrals_aa) std::cout << " "<<s ; std::cout << std::endl;
}

OldOperator::OldOperator()
  : m_Operator(dummyOperator)
{
  bracket_integrals_a=bracket_integrals_b=NULL;
  bracket_integrals_aa=bracket_integrals_ab=bracket_integrals_bb=NULL;
  loaded = false;
}

OldOperator::OldOperator(const OldOperator &source)
  : OrbitalSpace(source)
  , m_Operator(source.m_Operator)
  ,loaded(source.loaded)
  , coreEnergy(source.coreEnergy)
  , basisSize(source.basisSize), ijSize(source.ijSize), ijklSize(source.ijklSize)
{
 this->_copy(source);
}

OldOperator::OldOperator(const OldOperator &source, const bool forceSpinUnrestricted, const bool oneElectron, const bool twoElectron)
  : OrbitalSpace(source)
  , m_Operator(source.m_Operator)
  ,loaded(source.loaded)
  , coreEnergy(source.coreEnergy)
  , basisSize(source.basisSize), ijSize(source.ijSize), ijklSize(source.ijklSize)
{
  this->_copy(source,forceSpinUnrestricted,oneElectron,twoElectron);
}

OldOperator::OldOperator(const std::string special, const OldOperator &source, const bool forceSpinUnrestricted)
  : OrbitalSpace(source)
  , m_Operator(source.m_Operator)
  ,loaded(false)
  , coreEnergy(0)
  , basisSize(source.basisSize), ijSize(source.ijSize), ijklSize(source.ijklSize)
{
  this->_copy(source,forceSpinUnrestricted,false,false);
  if (special=="P" || special=="Q") {
      double fill=(special=="Q"? 0 : 1);
      integrals_a = new std::vector<double>(ijSize,0.0);
      if (spinUnrestricted)
        integrals_b = new std::vector<double>(ijSize,0.0);
      else
        integrals_b = integrals_a;
      for (auto aa=integrals_b->begin(); aa!= integrals_b->end(); aa++) *aa=0;
      for (auto aa=integrals_a->begin(); aa!= integrals_a->end(); aa++) *aa=0;
      for (unsigned int sym=0; sym<8; sym++) {
          for (size_t k=0; k<this->at(sym); k++)
            integrals_a->at(int1Index(offset(sym)+1+k,offset(sym)+1+k))=
            integrals_b->at(int1Index(offset(sym)+1+k,offset(sym)+1+k))=fill;
        }
      integrals_aa = NULL;
      integrals_ab = NULL;
      integrals_bb = NULL;
      bracket_integrals_aa = NULL;
      bracket_integrals_ab = NULL;
      bracket_integrals_bb = NULL;
//      constructBraKet(0,1);
      // determine the uncoupled orbital
      size_t uncoupled_orbital=0;
      double min_rowsum=1e50;
      for (unsigned int sym=0; sym<8; sym++) {
          for (size_t k=0; k<this->at(sym); k++) {
              double rowsum=0;
              for (size_t l=0; l<this->at(sym); l++)
                rowsum+=source.integrals_a->at(int1Index(k+offset(sym)+1,l+offset(sym)+1));
//              xout << "k="<<k<<", rowsum="<<rowsum<<std::endl;
              if (std::fabs(rowsum) < min_rowsum) {
                  min_rowsum=std::fabs(rowsum);
                  uncoupled_orbital=k+offset(sym)+1;
                }
            }
        }
      xout << "non-interacting orbital is "<<uncoupled_orbital<<std::endl;
      integrals_b->at(int1Index(uncoupled_orbital,uncoupled_orbital))=1-fill;
      loaded=true;
      bracket_integrals_a = new std::vector<double>(*integrals_a) ;
      bracket_integrals_b = new std::vector<double>(*integrals_b) ;
//        xout << "integrals_b"<<std::endl;
//        for (size_t k=0; k<6; k++) xout <<integrals_b->at(k)<<std::endl;
//        xout << "bracket_integrals_b"<<std::endl;
//        for (size_t k=0; k<6; k++) xout <<bracket_integrals_b->at(k)<<std::endl;
    }
}

OldOperator& OldOperator::operator=(const OldOperator &source)
{
 this->_copy(source);
 return *this;
}

void OldOperator::_copy(const OldOperator &source, const bool forceSpinUnrestricted, const bool oneElectron, const bool twoElectron)
{
//  xout << "Operator::_copy"<<std::endl<<source.str(2)<<std::endl;
  if (forceSpinUnrestricted) spinUnrestricted = true;
  bracket_integrals_a = bracket_integrals_b = NULL;
  if (loaded) {
       integrals_a = oneElectron ? new std::vector<double>(*source.integrals_a) : NULL;
       bracket_integrals_a = bracket_integrals_b = (oneElectron && source.bracket_integrals_a != NULL) ? new std::vector<double>(*source.bracket_integrals_a) : NULL;
     if (source.integrals_aa != NULL || spinUnrestricted) {
       integrals_aa = (twoElectron && source.integrals_aa != NULL) ? new std::vector<double>(*source.integrals_aa) : NULL;
       bracket_integrals_aa = (twoElectron && source.bracket_integrals_aa != NULL) ? new std::vector<double>(*source.bracket_integrals_aa) : NULL;
       bracket_integrals_ab = (twoElectron && source.bracket_integrals_ab != NULL)
           ?  new std::vector<double>(*source.bracket_integrals_ab)
           : NULL;
     }
     if (spinUnrestricted && ! source.spinUnrestricted) {
       integrals_b = oneElectron ? new std::vector<double>(*source.integrals_a): NULL;
       if (source.integrals_ab != NULL || spinUnrestricted) {
         integrals_ab = (twoElectron && source.integrals_aa != NULL)
             ? new std::vector<double>(*source.integrals_aa)
             : NULL;
         integrals_bb = (twoElectron && source.integrals_aa != NULL) ? new std::vector<double>(*source.integrals_aa) : NULL;
         bracket_integrals_bb = (twoElectron && source.bracket_integrals_aa != NULL) ? new std::vector<double>(*source.bracket_integrals_aa) : NULL;
         bracket_integrals_b = (oneElectron && source.bracket_integrals_a != NULL) ? new std::vector<double>(*source.bracket_integrals_a) : NULL;
//         xout << "Operator::_copy unrestrict bracket_integrals.b "<<bracket_integrals_b<<std::endl;
//         if (bracket_integrals_b!=NULL)
//         for (size_t i=0; i<bracket_integrals_b->size();i++) xout << " " <<(*bracket_integrals_b)[i]; xout <<std::endl;
       }
     } else  if (spinUnrestricted) {
       integrals_b = oneElectron ? new std::vector<double>(*source.integrals_b): NULL;
       if (source.integrals_ab != NULL || spinUnrestricted) {
         integrals_ab = (twoElectron && source.integrals_ab != NULL)
             ? new std::vector<double>(*source.integrals_ab)
             : NULL;
         integrals_bb = (twoElectron && source.integrals_bb != NULL) ? new std::vector<double>(*source.integrals_bb) : NULL;
         bracket_integrals_bb = (twoElectron && source.bracket_integrals_bb != NULL) ? new std::vector<double>(*source.bracket_integrals_bb) : NULL;
         bracket_integrals_b = (oneElectron && source.bracket_integrals_b != NULL) ? new std::vector<double>(*source.bracket_integrals_b) : NULL;
//         xout << "Operator::_copy bracket_integrals.b "<<bracket_integrals_b<<std::endl;
//         if (bracket_integrals_b!=NULL)
//         for (size_t i=0; i<bracket_integrals_b->size();i++) xout << " " <<(*bracket_integrals_b)[i]; xout <<std::endl;
       }
     } else {
       integrals_b = integrals_a;
       integrals_ab = integrals_aa;
       integrals_bb = integrals_aa;
       bracket_integrals_bb = bracket_integrals_aa;
       bracket_integrals_b = bracket_integrals_a;
     }
  }
}

OldOperator::~OldOperator() {
  deconstructBraKet();
}

#define del(x) if (x != NULL && x->size()) delete x; x=NULL;
#define del2(x,y) if (x != y) del(x); del(y); del(x);
void OldOperator::constructBraKet(int neleca, int nelecb)
{
  throw std::logic_error("OldOperator::constructBraKet obsolete");
  deconstructBraKet();
  // construct <ik||jl> = (ij|kl) - (il|kj)
  bracket_integrals_aa = new std::vector<double>(pairSpace[-1].total(0,0),0.0);
  bracket_integrals_ab = new std::vector<double>(pairSpace[0].total(0,0),0.0);
  if (spinUnrestricted) {
    bracket_integrals_bb = new std::vector<double>(pairSpace[-1].total(0,0),0.0);
  } else {
    bracket_integrals_bb = bracket_integrals_aa;
  }
  for (unsigned int symik=0; symik<8; symik++) {
    for (unsigned int symi=0; symi<8; symi++) {
      for (unsigned int symj=0; symj<8; symj++) {
        unsigned int symk = symi^symik;
        unsigned int syml = symj^symik;
        unsigned int symij = symi^symj;
        if ((*this)[symi]==0) continue;
        if ((*this)[symj]==0) continue;
        if ((*this)[symk]==0) continue;
        if ((*this)[syml]==0) continue;
        //                xout << "symi symj symk syml" << symi<<symj<<symk<<syml;
        //                xout << "; pairSpace offset="<<pairSpace[0].offset(0,symik)<<"; within-pair offsets: "
        //                     <<offset(symik,symj,0) <<" "
        //                    <<offset(symik,symj,0)
        //                   <<std::endl;
        for (size_t i=0; i< (*this)[symi] ; i++) {
          for (size_t j=0; j< (*this)[symj] ; j++) {
            for (size_t k=0; k< (*this)[symk] ; k++) {
              for (size_t l=0; l< (*this)[syml] ; l++) {
                if (symij==0)
                {
                  //                                    xout << "i,j,k,l,toaddress, fromaddress "<<i<<" "<<j<<" "<<k<<" "<<l<<" " <<
                  //                                            pairSpace[0].offset(0,symik)
                  //                                            + (offset(symik,symi,0)+i+k*at(symi))*total(symik,0)
                  //                                            + (offset(symik,symj,0)+j+l*at(symj))
                  //                                         <<" " <<
                  //                                           pairSpace[1].offset(0,symij)
                  //                                           + (offset(symij,symi,1)+((i>j) ? (i*(i+1))/2+j : (j*(j+1))/2+i)) * total(symij,1)
                  //                                           + (offset(symij,symk,1)+((k>l) ? (k*(k+1))/2+l : (l*(l+1))/2+k))
                  //                                        <<"="<<
                  //                                          integrals_ab->at(pairSpace[0].offset(0,symij)
                  //                                          + (offset(symij,symi,1)+((i>j) ? (i*(i+1))/2+j : (j*(j+1))/2+i)) * total(symij,1)
                  //                                          + (offset(symij,symk,1)+((k>l) ? (k*(k+1))/2+l : (l*(l+1))/2+k))
                  //                                          )
                  //                                            <<std::endl;
                  size_t address =
                      pairSpace[0].offset(0,symik)
                      + (offset(symik,symi,0)+i+k*at(symi))*total(symik,0)
                      + (offset(symik,symj,0)+j+l*at(symj));
                  bracket_integrals_ab->at(address) =
                      integrals_ab->at(pairSpace[1].offset(0,symij)
                      + (offset(symij,symi,1)+((i>j) ? (i*(i+1))/2+j : (j*(j+1))/2+i)) * total(symij,1)
                      + (offset(symij,symk,1)+((k>l) ? (k*(k+1))/2+l : (l*(l+1))/2+k))
                      );
                  if (nelecb !=0 && k==l) { // absorb 1-electron integrals into two-electron for a given number of electrons
                      bracket_integrals_ab->at(address) +=
                          integrals_a->at(
                            offset(0,symi,1) + (i>j ? (i*(i+1))/2+j : (j*(j+1))/2+i)
                            ) / (.5*nelecb);
//                      xout << "i j k l hij"<<i<<j<<k<<l<<" "<< integrals_a->at( offset(0,symi,1) + (i>j ? (i*(i+1))/2+j : (j*(j+1))/2+i) ) <<std::endl;
                  }
                  if (neleca !=0 && i==j) { // absorb 1-electron integrals into two-electron for a given number of electrons
                      bracket_integrals_ab->at(address) +=
                          integrals_a->at(
                            offset(0,symk,1) + (k>l ? (k*(k+1))/2+l : (l*(l+1))/2+k)
                            ) / (.5*neleca);
//                      xout << "i j k l hkl"<<i<<j<<k<<l<<" "<< integrals_a->at( offset(0,symk,1) + (k>l ? (k*(k+1))/2+l : (l*(l+1))/2+k) )  <<std::endl;
                  }
                }
                else
                {
                  //                                    xout << "i,j,k,l,toaddress, fromaddress "<<i<<" "<<j<<" "<<k<<" "<<l<<" " <<
                  //                                            pairSpace[0].offset(0,symik)
                  //                                            + (offset(symik,symi,0)+i+k*at(symi))*total(symik,0)
                  //                                            + (offset(symik,symj,0)+j+l*at(symj))
                  //                                         <<" " <<
                  //                                           pairSpace[1].offset(0,symij)
                  //                                           + ((symi > symj) ?
                  //                                                  offset(symij,symi,1)+j*at(symi)+i
                  //                                                : offset(symij,symj,1)+i*at(symj)+j
                  //                                                  ) * total(symij,1)
                  //                                           + ((symk > syml) ?
                  //                                                  offset(symij,symk,1)+l*at(symk)+k
                  //                                                : offset(symij,syml,1)+k*at(syml)+l)
                  //                                        <<"="<<
                  //                                          integrals_ab->at(pairSpace[1].offset(0,symij)
                  //                                          + ((symi > symj) ?
                  //                                                 offset(symij,symi,1)+j*at(symi)+i
                  //                                               : offset(symij,symj,1)+i*at(symj)+j
                  //                                                 ) * total(symij,1)
                  //                                          + ((symk > syml) ?
                  //                                                 offset(symij,symk,1)+l*at(symk)+k
                  //                                               : offset(symij,syml,1)+k*at(syml)+l)
                  //                                          )
                  //                                            <<std::endl;
                  size_t address =
                      pairSpace[0].offset(0,symik)
                      + (offset(symik,symi,0)+i+k*at(symi))*total(symik,0)
                      + (offset(symik,symj,0)+j+l*at(symj));
                  bracket_integrals_ab->at( address ) =
                      integrals_ab->at(pairSpace[1].offset(0,symij)
                      + ((symi > symj) ?
                           offset(symij,symi,1)+j*at(symi)+i
                         : offset(symij,symj,1)+i*at(symj)+j
                           ) * total(symij,1)
                      + ((symk > syml) ?
                           offset(symij,symk,1)+l*at(symk)+k
                         : offset(symij,syml,1)+k*at(syml)+l
                           )
                      );
                }
              }
            }
          }
        }
        del(bracket_integrals_a);
	// xout << "considering 1-e bracket " << neleca<<nelecb<<std::endl;
        if (nelecb == 0 && integrals_a != NULL) bracket_integrals_a = new std::vector<double>(*integrals_a);
        del(bracket_integrals_b);
        if (neleca == 0 && integrals_b != NULL) bracket_integrals_b = new std::vector<double>(*integrals_b);
	// xout << "bracket_integrals_a "<<bracket_integrals_a<<std::endl;
	// xout << "bracket_integrals_b "<<bracket_integrals_b<<std::endl;
//        if (neleca != 0 && integrals_b != NULL) {
//          xout << "copied into bracket_integrals_b "<<std::endl;
//          for (unsigned int i=0; i<bracket_integrals_b->size(); i++) xout << " " << (*bracket_integrals_b)[i]; xout <<std::endl;
//        }
        // alpha-alpha and beta-beta
        unsigned int symil = symi^syml;
        unsigned int symjl = symj^syml;
        if (symi >= symk && symj >=syml) {
          for (size_t i=0; i< (*this)[symi] ; i++) {
            size_t klimit = symi == symk ? i : (*this)[symk];
            for (size_t j=0; j< (*this)[symj] ; j++) {
              size_t llimit = symj == syml ? j : (*this)[syml];
              for (size_t k=0; k< klimit ; k++) {
                for (size_t l=0; l< llimit ; l++) {
                  size_t ijkl =
                      pairSpace[1].offset(0,symij) +
                      ((symij==0) ?
                         ((offset(symij,symi,1)+((i>j) ? (i*(i+1))/2+j : (j*(j+1))/2+i)) * total(symij,1)
                          +(offset(symij,symk,1)+((k>l) ? (k*(k+1))/2+l : (l*(l+1))/2+k)))
                       : (
                         ((symi > symj) ?
                            (offset(symij,symi,1)+j*at(symi)+i)
                          : (offset(symij,symj,1)+i*at(symj)+j)
                            ) * total(symij,1)
                         +((symk > syml) ?
                             (offset(symij,symk,1)+l*at(symk)+k)
                           : (offset(symij,syml,1)+k*at(syml)+l))
                         )) ;
                  size_t ilkj =
                      pairSpace[1].offset(0,symil) +
                      ((symil==0) ?
                         ((offset(symil,symi,1)+((i>l) ? (i*(i+1))/2+l : (l*(l+1))/2+i)) * total(symil,1)
                          +(offset(symil,symk,1)+((k>j) ? (k*(k+1))/2+j : (j*(j+1))/2+k)))
                       : (
                         ((symi > syml) ?
                            (offset(symil,symi,1)+l*at(symi)+i)
                          : (offset(symil,syml,1)+i*at(syml)+l)
                            ) * total(symil,1)
                      + ((symk > symj) ?
                          (offset(symil,symk,1)+j*at(symk)+k)
                        : (offset(symil,symj,1)+k*at(symj)+j))
                      ));
                  size_t ik = offset(symik,symi,-1) +  ((symi == symk) ? ((i*(i-1))/2+k) : (i+k*at(symi))) ;
                  size_t jl = offset(symjl,symj,-1) +  ((symj == syml) ? ((j*(j-1))/2+l) : (j+l*at(symj))) ;
//                  xout << "i,j,k,l,ijkl,ilkj,ik,jl: "<<i<<" "<<j<<" "<<k<<" "<<l<<" "
//                          <<ijkl<<" "
//                          <<ilkj<<" "
//                          <<ik<<" "
//                          <<jl<<" destination: "
//                                                       <<pairSpace[-1].offset(0,symik,0) + ik*total(symik,-1) + jl  << " value: "
//                                                       <<integrals_aa->at(ijkl) -integrals_aa->at(ilkj)<<std::endl;
                  size_t address = pairSpace[-1].offset(0,symik,0) + ik*total(symik,-1) + jl ;
                  bracket_integrals_aa->at(address) = integrals_aa->at(ijkl) -integrals_aa->at(ilkj);
                  bracket_integrals_bb->at(address) = integrals_bb->at(ijkl) -integrals_bb->at(ilkj);
                }
              }
            }
          }
        }
      }
    }
  }

}


void OldOperator::deconstructBraKet()
{
  if (! loaded) return;
  del2(bracket_integrals_b, bracket_integrals_a);
  del2(bracket_integrals_aa, bracket_integrals_bb);
  del(bracket_integrals_ab);
}

void OldOperator::unload() {
  if (loaded) {
    delete integrals_a;
    delete integrals_aa;
    delete bracket_integrals_aa;
    delete bracket_integrals_ab;
    if (spinUnrestricted) {
      delete integrals_b;
      delete integrals_ab;
      delete integrals_bb;
      delete bracket_integrals_bb;
    }
  }
  loaded=false;
}

std::string OldOperator::str(int verbosity, unsigned int columns) const
{
  std::ostringstream o;
  o << OrbitalSpace::str(verbosity>3 ? verbosity : 0);
  if (verbosity>=0) {
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
    if (integrals_aa!=NULL || integrals_ab != NULL || integrals_bb != NULL) {
      std::vector<double>* zero =  new std::vector<double>(ijklSize,(double)0);
      std::vector<double>* iaa = integrals_aa==NULL ? zero : integrals_aa;
      std::vector<double>* iab = integrals_ab==NULL ? zero : integrals_ab;
      std::vector<double>* ibb = integrals_bb==NULL ? zero : integrals_bb;
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
              if (iaa->at(ijkl) != (double)0 || iab->at(ijkl) != (double)0 || ibb->at(ijkl) != (double)0) o<<std::endl<<std::setw(4)<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<std::setw(precision+7)<<iaa->at(ijkl)<<" "<<iab->at(ijkl)<<" "<<ibb->at(ijkl);
            }
          }
        }
      }
    }
    if (bracket_integrals_ab != NULL) {
      o<<std::endl << "2-electron integrals (bracket form, ab):";
      for (int symij=0; symij<8; symij++) {
        if (total(symij,0))
          o<<std::endl <<"symmetry block " << symij;
        for (size_t ij=0; ij< total(symij,0); ij++) {
          o<<std::endl;
          for (size_t kl=0; kl< total(symij,0); kl++) {
            o << bracket_integrals_ab->at( pairSpace.find(0)->second.offset(0,symij,0) + ij * total(symij,0) + kl ) <<" ";
          }
        }
      }
    }
    if (bracket_integrals_aa != NULL && bracket_integrals_bb != NULL) {
      o<<std::endl << "2-electron integrals (bracket form, antisymmetrised, {aa,bb}):";
      for (int symij=0; symij<8; symij++) {
        if (total(symij,-1))
          o<<std::endl <<"symmetry block " << symij;
        for (size_t ij=0; ij< total(symij,-1); ij++) {
          o<<std::endl;
          for (size_t kl=0; kl< total(symij,-1); kl++) {
            o <<"{"
             << bracket_integrals_aa->at( pairSpace.find(-1)->second.offset(0,symij,0) + ij * total(symij,-1) + kl ) <<","
             << bracket_integrals_bb->at( pairSpace.find(-1)->second.offset(0,symij,0) + ij * total(symij,-1) + kl ) <<"} ";
          }
        }
      }
    }
  }
  return o.str();
}

size_t OldOperator::int1Index(unsigned int i, unsigned int j) const {
  return pairIndex(i,j,1);
}

size_t OldOperator::int2Index(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
{
  return quadIndex(i,j,k,l,1,0);
}

std::vector<double> OldOperator::int1(int spin) const
{
  std::vector<double> result(basisSize,(double)0);
  std::vector<double> * integrals = spin < 0 ? integrals_b : integrals_a;
  for (unsigned int i=0; i < basisSize; i++) {
    result[i] = integrals->at(int1Index(i+1,i+1));
  }
  return result;
}


std::vector<double> OldOperator::intJ(int spini, int spinj) const
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

std::vector<double> OldOperator::intK(int spin) const
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

OldOperator OldOperator::FockOperator(const Determinant &reference) const
{
  OldOperator f;
  for (int i=0; i<8; i++)
    f[i]=at(i);
  f.calculateOffsets();
  // xout << "FockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout << "FockOperator Reference beta: "<<reference.stringBeta<<std::endl;
  bool closed = reference.stringAlpha==reference.stringBeta;
//  xout << "FockOperator Reference alpha=beta: "<<closed<<std::endl;
  f.spinUnrestricted = spinUnrestricted || ! closed;
//  xout << "FockOperator spinUnrestricted="<<f.spinUnrestricted<<std::endl;
  auto refAlphaOrbitals=reference.stringAlpha.orbitals();
  auto refBetaOrbitals=reference.stringBeta.orbitals();
  f.coreEnergy = coreEnergy;
  f.basisSize = basisSize;
  f.ijklSize = ijklSize;
  f.ijSize = ijSize;
  f.orbital_symmetries = orbital_symmetries;
  f.integrals_a = new std::vector<double>(ijSize,0.0);
  *f.integrals_a = *integrals_a;
  // xout <<"reference.stringAlpha.orbitals ";for (size_t i=0; i < reference.stringAlpha.orbitals().size(); i++) xout <<reference.stringAlpha.orbitals()[i]<<" ";xout <<std::endl;
  // for (std::vector<unsigned int>::const_iterator o=reference.stringAlpha.orbitals().begin(); o != reference.stringAlpha.orbitals().end(); o++)
  for (auto o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
  {
    // xout << "FockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout<< "f alpha, alpha occ: " <<*o << std::endl;
    for (unsigned int i=1; i<=basisSize; i++)
      for (unsigned int j=1; j<=i; j++) {
        if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
        (*f.integrals_a)[int1Index(i,j)] += (*integrals_aa)[int2Index(i,j,*o,*o)] - (*integrals_aa)[int2Index(i,*o,*o,j)];
      }
  }
  for (auto o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
  {
    // xout<< "f alpha, beta occ: " <<*o << std::endl;
    for (unsigned int i=1; i<=basisSize; i++)
      for (unsigned int j=1; j<=i; j++) {
        if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
        (*f.integrals_a)[int1Index(i,j)] += (*integrals_ab)[int2Index(i,j,*o,*o)];
      }
  }
  f.integrals_aa = NULL;
  f.integrals_ab = NULL;
  f.integrals_bb = NULL;
  f.bracket_integrals_aa = NULL;
  f.bracket_integrals_ab = NULL;
  f.bracket_integrals_bb = NULL;
  f.bracket_integrals_a = (f.integrals_a != NULL) ? new std::vector<double>(*f.integrals_a) : NULL;
  // xout << "in FockOperator, after alpha f="; for (size_t ij=0; ij< f.integrals_a->size(); ij++) xout <<" "<<(*f.integrals_a)[ij]; xout <<std::endl;
  if (f.spinUnrestricted) {
    f.integrals_b = new std::vector<double>(ijSize,0.0);
    *f.integrals_b = *integrals_b;
    for (auto o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
    {
      // xout<< "f beta, beta occ: " <<*o << std::endl;
      for (unsigned int i=1; i<=basisSize; i++)
        for (unsigned int j=1; j<=i; j++) {
          if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
          (*f.integrals_b)[int1Index(i,j)] += (*integrals_bb)[int2Index(i,j,*o,*o)] - (*integrals_bb)[int2Index(i,*o,*o,j)];
        }
    }
    for (auto o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
    {
      // xout<< "f beta, alpha occ: " <<*o << std::endl;
      for (unsigned int i=1; i<=basisSize; i++)
        for (unsigned int j=1; j<=i; j++) {
          if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
          (*f.integrals_b)[int1Index(i,j)] += (*integrals_ab)[int2Index(*o,*o,i,j)];
        }
    }
    f.bracket_integrals_b = f.integrals_b != NULL ? new std::vector<double>(*f.integrals_b) : NULL;
  } else {
    f.integrals_b = f.integrals_a;
    f.bracket_integrals_b = f.bracket_integrals_a;
  }
  f.loaded = true;
  return f;
}

OldOperator OldOperator::sameSpinOperator(const Determinant &reference) const
{
  OldOperator result = *this;
//  xout << "result when initialized: "<<result.str(2)<<std::endl;
  result.spinUnrestricted = true;
  if (!spinUnrestricted) *(result.integrals_b = new std::vector<double>(integrals_a->size())) = *result.integrals_a;
//  xout << "sameSpinOperator, old integrals_a="<<&integrals_a[0]<<", new integrals_a ="<<&result.integrals_a[0]<<std::endl;
//  xout << "sameSpinOperator, old integrals_b="<<&integrals_b[0]<<", new integrals_b ="<<&result.integrals_b[0]<<std::endl;
  Determinant ra = reference; ra.stringBeta.nullify();
//  xout << "this before alpha fock: "<<str(2)<<std::endl;
//  xout << "result before alpha fock: "<<result.str(2)<<std::endl;
  OldOperator f = this->FockOperator(ra);
//  xout << "this after alpha fock: "<<str(2)<<std::endl;
//  xout << "result before alpha: "<<result.str(2)<<std::endl;
  for (size_t i=0; i<integrals_a->size(); i++) {
    result.integrals_a->at(i) = this->integrals_a->at(i) - f.integrals_a->at(i);
//    xout << " result a = " << result.integrals_a->at(i) <<"=" << this->integrals_a->at(i) <<"-"<< f.integrals_a->at(i)<<std::endl;
  }
//  xout << "this after alpha: "<<str(2)<<std::endl;
//  xout << "result after alpha: "<<result.str(2)<<std::endl;
  ra = reference; ra.stringAlpha.nullify();
  f = this->FockOperator(ra);
//  xout << "sameSpinOperator, after beta Fock, f at "<<&f<<", f.bracket_integrals_a at "<<f.bracket_integrals_a<<std::endl;
//  xout << "f.bracket_integrals_a->size() "<<f.bracket_integrals_a->size()<<std::endl;
//  xout << "sameSpinOperator, fock integrals_a="<<&f.integrals_a[0]<<", fock integrals_b ="<<&f.integrals_b[0]<<std::endl;
  for (size_t i=0; i<integrals_b->size(); i++) {
    result.integrals_b->at(i) = this->integrals_b->at(i) - f.integrals_b->at(i);
//    xout << " result b = " << result.integrals_b->at(i) <<"=" << this->integrals_b->at(i) <<"-"<< f.integrals_b->at(i)<<std::endl;
}
  if (result.bracket_integrals_ab != NULL) {
//    xout << "result.bracket_integrals_ab isn't null; size="<<result.bracket_integrals_ab->size()<<std::endl;
    // not very satisfactory
    result.bracket_integrals_ab->assign(result.bracket_integrals_ab->size(),(double)0);
//    delete result.bracket_integrals_ab;
//    result.bracket_integrals_ab = NULL;
  }
  if (spinUnrestricted) delete result.integrals_ab;
  result.integrals_ab = NULL;
//  xout << "result on return: "<<result.str(2)<<std::endl;
  return result;
}

#include <assert.h>
OldOperator& OldOperator::plusminusOperator(const OldOperator &other, const char operation)
{
//  if (! compatible(other)) throw std::logic_error("attempt to add incompatible Operator objects");
  assert(this->spinUnrestricted || ! other.spinUnrestricted);
  if (operation == '+')
    coreEnergy += other.coreEnergy;
  else if (operation == '-')
    coreEnergy -= other.coreEnergy;
  plusminusEqualsHelper(this->integrals_a, other.integrals_a,operation);
  if (integrals_a != integrals_b)
    plusminusEqualsHelper(this->integrals_b, other.integrals_b,operation);
  plusminusEqualsHelper(this->integrals_aa, other.integrals_aa,operation);
  if (integrals_aa != integrals_ab)
    plusminusEqualsHelper(this->integrals_ab, other.integrals_ab,operation);
  plusminusEqualsHelper(this->bracket_integrals_a, other.bracket_integrals_a,operation);
  if (bracket_integrals_a != bracket_integrals_b)
    plusminusEqualsHelper(this->bracket_integrals_b, other.bracket_integrals_b,operation);
  if (integrals_aa != integrals_bb)
    plusminusEqualsHelper(this->integrals_bb, other.integrals_bb,operation);
  plusminusEqualsHelper(this->bracket_integrals_aa, other.bracket_integrals_aa,operation);
  plusminusEqualsHelper(this->bracket_integrals_ab, other.bracket_integrals_ab,operation);
  if (bracket_integrals_aa != bracket_integrals_bb)
    plusminusEqualsHelper(this->bracket_integrals_bb, other.bracket_integrals_bb,operation);
  return *this;
}
OldOperator& OldOperator::operator+=(const OldOperator &other)
{
  return plusminusOperator(other,'+');
}
OldOperator& OldOperator::operator-=(const OldOperator &other)
{
  return plusminusOperator(other,'-');
}

#include <cmath>
void OldOperator::plusminusEqualsHelper(std::vector<double> *&me,
                                    std::vector<double> * const &other,
                                        const char operation)
{
  if (other == NULL) return;
  size_t n = other->size();
  double valmax=(double)0;
  for (size_t i=0; i<n; i++)
    if (valmax < std::abs(other->at(i))) valmax = std::abs(other->at(i));
  if (valmax == (double)0) return;
  if (me == NULL) me = new std::vector<double>(n,(double)0);
  if (operation == '+')
    for (size_t i=0; i<n; i++)
      me->at(i) += other->at(i);
  else if (operation == '-')
    for (size_t i=0; i<n; i++)
      me->at(i) -= other->at(i);
}
void OldOperator::starEqualsHelper(std::vector<double> *&me,
                                    const double factor)
{
  if (factor == (double)1) return;
  if (me == NULL) return;
  size_t n = me->size();
  for (size_t i=0; i<n; i++)
    me->at(i) *= factor ;
}

OldOperator& OldOperator::operator*=(const double factor)
{
  coreEnergy *= factor;
  starEqualsHelper(this->integrals_a, factor);
  if (integrals_a != integrals_b)
    starEqualsHelper(this->integrals_b, factor);
  starEqualsHelper(this->integrals_aa, factor);
  if (integrals_aa != integrals_ab)
    starEqualsHelper(this->integrals_ab, factor);
  if (integrals_aa != integrals_bb)
    starEqualsHelper(this->integrals_bb, factor);
  starEqualsHelper(this->bracket_integrals_aa, factor);
  starEqualsHelper(this->bracket_integrals_ab, factor);
  if (bracket_integrals_aa != bracket_integrals_bb)
    starEqualsHelper(this->bracket_integrals_bb, factor);
  return *this;

}

OldOperator gci::operator+(const OldOperator &h1, const OldOperator &h2)
{
  OldOperator result = h1;
  return result += h2;
}
OldOperator gci::operator-(const OldOperator &h1, const OldOperator &h2)
{
  OldOperator result = h1;
  return result -= h2;
}
OldOperator gci::operator*(const OldOperator &h1, const double factor)
{
  OldOperator result = h1;
  return result *= factor;
}

void OldOperator::rotate1(std::vector<double>* integrals,std::vector<double> const * rot)
{
  if (integrals == NULL) return;
//  xout << "rotate1: input"; for (std::vector<double>::const_iterator i=integrals->begin(); i!=integrals->end(); i++) xout <<" "<<*i; xout <<std::endl;
//  xout << "rotate1: rotator"; for (std::vector<double>::const_iterator i=rot->begin(); i!=rot->end(); i++) xout <<" "<<*i; xout <<std::endl;
  // for (unsigned int i=1; i<=basisSize; i++)
    // for (unsigned int k=1; k<=basisSize; k++)
      // xout << k<<" "<<i<<" "<<pairIndex(k,i)<<" "<< (*rot)[pairIndex(k,i)]<<std::endl;
  std::vector<double> t1(total(0,0));
  for (unsigned int i=1; i<=basisSize; i++) {
    for (unsigned int k=1; k<=basisSize; k++) {
      if (orbital_symmetries[i-1]!=orbital_symmetries[k-1]) continue;
      t1[pairIndex(k,i)]=(double)0;
      for (unsigned int j=1; j<=basisSize; j++) {
        if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
        t1[pairIndex(k,i)]+=(*rot)[pairIndex(j,k)] * (*integrals)[int1Index(j,i)];
      }
    }
  }
//  xout << "rotate1: t1"; for (std::vector<double>::const_iterator i=t1.begin(); i!=t1.end(); i++) xout <<" "<<*i; xout <<std::endl;
  for (unsigned int i=1; i<=basisSize; i++) {
    for (unsigned int k=1; k<=i; k++) {
      if (orbital_symmetries[i-1]!=orbital_symmetries[k-1]) continue;
      (*integrals)[int1Index(k,i)]=(double)0;
      for (unsigned int j=1; j<=basisSize; j++) {
        if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
        (*integrals)[int1Index(k,i)] +=  (*rot)[pairIndex(j,k)] * t1[pairIndex(i,j)];
      }
    }
  }
//  xout << "rotate1: output"; for (std::vector<double>::const_iterator i=integrals->begin(); i!=integrals->end(); i++) xout <<" "<<*i; xout <<std::endl;
}

void OldOperator::rotate2(std::vector<double>* integrals,std::vector<double> const * rot1, std::vector<double> const * rot2)
{
  if (integrals == NULL) return;
  std::vector<double> t1(std::pow(this->at(0),4));
  std::vector<double> t2(std::pow(this->at(0),4));
  // for (unsigned int isym=0; isym<8; isym++) {
  //   xout <<"rotation for symmetry "<<isym<<std::endl;
  //   for (unsigned int i=0; i<this->at(isym); i++) {
  //     for (unsigned int j=0; j<this->at(isym); j++)
  // 	xout <<" "<<(*rot1)[offset(0,isym)+i+j*this->at(isym)];
  //     xout <<std::endl;
  //   }
  // }
  for (unsigned int ijsym=0; ijsym<8; ijsym++) {
    for (unsigned int isym=0; isym<8; isym++) {
      unsigned int ni=this->at(isym);
      unsigned int jsym=ijsym^isym;
      if (jsym > isym) continue;
      unsigned int nj=this->at(jsym);
      for (unsigned int ksym=0; ksym<8; ksym++) {
        unsigned int nk=this->at(ksym);
        unsigned int lsym=ijsym^ksym;
        unsigned int nl=this->at(lsym);
        if (lsym <= ksym && ni*nj*nk*nl > 0) {
        for (unsigned int l=0; l<nl; l++) {
          for (unsigned int k=0; k<nk; k++) {
            for (unsigned int j=0; j<nj; j++) {
              for (unsigned int i=0; i<ni; i++) {
                t1[i+j*ni+k*ni*nj+l*ni*nj*nk]=(double)0;
                for (unsigned int i1=0; i1<ni; i1++) {
                  t1[i+j*ni+k*ni*nj+l*ni*nj*nk] += (*rot1)[offset(0,isym)+i1+i*ni]
                      * (*integrals)[int2Index(i1+offset(isym)+1,j+offset(jsym)+1,k+offset(ksym)+1,l+offset(lsym)+1)];
                }
              }
            }
          }
        }
        for (unsigned int l=0; l<nl; l++) {
          for (unsigned int k=0; k<nk; k++) {
            for (unsigned int j=0; j<nj; j++) {
              for (unsigned int i=0; i<ni; i++) {
                t2[i+j*ni+k*ni*nj+l*ni*nj*nk]=(double)0;
                for (unsigned int j1=0; j1<nj; j1++) {
                  t2[i+j*ni+k*ni*nj+l*ni*nj*nk] += (*rot1)[offset(0,jsym)+j1+j*nj]
                   * t1[i+j1*ni+k*ni*nj+l*ni*nj*nk];
                }
              }
            }
          }
        }
        for (unsigned int l=0; l<nl; l++) {
          for (unsigned int k=0; k<nk; k++) {
            for (unsigned int j=0; j<nj; j++) {
              for (unsigned int i=0; i<ni; i++) {
                t1[i+j*ni+k*ni*nj+l*ni*nj*nk]=(double)0;
                for (unsigned int k1=0; k1<nk; k1++) {
                  t1[i+j*ni+k*ni*nj+l*ni*nj*nk] += (*rot2)[offset(0,ksym)+k1+k*nk]
                   * t2[i+j*ni+k1*ni*nj+l*ni*nj*nk];
                }
              }
            }
          }
        }
        for (unsigned int l=0; l<nl; l++) {
          for (unsigned int k=0; k<(ijsym ? nk : l+1); k++) {
            for (unsigned int j=0; j<nj; j++) {
              for (unsigned int i=0; i<(ijsym ? ni : j+1); i++) {
                (*integrals)[int2Index(i+offset(isym)+1,j+offset(jsym)+1,k+offset(ksym)+1,l+offset(lsym)+1)] = (double)0;
                for (unsigned int l1=0; l1<nl; l1++) {
                (*integrals)[int2Index(i+offset(isym)+1,j+offset(jsym)+1,k+offset(ksym)+1,l+offset(lsym)+1)] +=
                  (*rot2)[offset(0,lsym)+l1+l*nl]
                   * t1[i+j*ni+k*ni*nj+l1*ni*nj*nk];
                }
//                xout << "new integral "<<i<<j<<k<<l<<" "<<int2Index(i+offset(isym)+1,j+offset(jsym)+1,k+offset(ksym)+1,l+offset(lsym)+1)<<(*integrals)[int2Index(i+offset(isym)+1,j+offset(jsym)+1,k+offset(ksym)+1,l+offset(lsym)+1)]<<std::endl;
              }
            }
          }
        }
        }
      }
    }
  }
}

void OldOperator::rotate(std::vector<double> const * rota, std::vector<double> const * rotb)
{
  if (rotb == NULL) rotb=rota;
  xout << "Operator::rotate"<<std::endl;
  // if (integrals_a != NULL) xout << "integrals_a "<<integrals_a<<std::endl;
  // if (true || integrals_b != NULL) xout << "integrals_b "<<integrals_b<<std::endl;
  // if (integrals_aa != NULL) xout << "integrals_aa "<<integrals_aa<<std::endl;
  // if (integrals_ab != NULL) xout << "integrals_ab "<<integrals_ab<<std::endl;
  // if (integrals_bb != NULL) xout << "integrals_bb "<<integrals_bb<<std::endl;
  // if (true || bracket_integrals_a != NULL) xout << "bracket_integrals_a "<<bracket_integrals_a<<std::endl;
  // if (true || bracket_integrals_b != NULL) xout << "bracket_integrals_b "<<bracket_integrals_b<<std::endl;
  // if (bracket_integrals_aa != NULL) xout << "bracket_integrals_aa "<<bracket_integrals_aa<<std::endl;
  // if (bracket_integrals_ab != NULL) xout << "bracket_integrals_ab "<<bracket_integrals_ab<<std::endl;
  // if (bracket_integrals_bb != NULL) xout << "bracket_integrals_bb "<<bracket_integrals_bb<<std::endl;

//  xout << "Unrotated hamiltonian" <<std::endl << str() << std::endl;
  xout << "rotate: rota"; for (std::vector<double>::const_iterator i=rota->begin(); i!=rota->end(); i++) xout <<" "<<*i; xout <<std::endl;
//  xout << "rotate: rotb"; for (std::vector<double>::const_iterator i=rotb->begin(); i!=rotb->end(); i++) xout <<" "<<*i; xout <<std::endl;
  rotate1(integrals_a,rota);
  if (integrals_a != integrals_b) rotate1(integrals_b,rotb);
  deconstructBraKet();
//  xout << "Rotated 1-electron hamiltonian" <<std::endl << str() << std::endl;
  rotate2(integrals_aa,rota,rota);
//  xout << "Rotated hamiltonian after aa" <<std::endl << str() << std::endl;
  if (integrals_aa != integrals_ab) rotate2(integrals_ab,rota,rotb);
//  xout << "Rotated hamiltonian after ab" <<std::endl << str() << std::endl;
  if (integrals_aa != integrals_bb) rotate2(integrals_bb,rotb,rotb);
//  xout << "Rotated hamiltonian" <<std::endl << str() << std::endl;
  constructBraKet();
  return;
}

#include "memory.h"
void OldOperator::rotate(SMat const * rota, SMat const * rotb)
{
  SMat rotan = rota->desymmetrise();
  std::vector<double> rotanv;
  memory::vector<double>* dd;
  dd=rotan.data();
  for (size_t i=0; i<rotan.size(); i++)
    rotanv.push_back((*dd)[i]);
  if (rotb != NULL) {
    SMat rotbn = rotb->desymmetrise();
    std::vector<double> rotbnv;
    dd=rotbn.data();
    for (size_t i=0; i<rotbn.size(); i++)
      rotbnv.push_back((*dd)[i]);
    rotate(&rotanv,&rotbnv);
  }
  else
    rotate(&rotanv,NULL);
}
