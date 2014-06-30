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
  bracket_integrals_a=bracket_integrals_b=NULL;
  bracket_integrals_aa=bracket_integrals_ab=bracket_integrals_bb=NULL;
  loaded = false;
  if (filename != "") load(filename);
  xout <<"Hamiltonian filename constructor this="<<this<<", filename="<<filename<<", loaded="<<loaded<<std::endl;
}

Hamiltonian::Hamiltonian(FCIdump* dump) : OrbitalSpace(dump)
{
  xout <<"Hamiltonian FCIdump constructor this="<<this<<std::endl;
  bracket_integrals_a=bracket_integrals_b=NULL;
  bracket_integrals_aa=bracket_integrals_ab=bracket_integrals_bb=NULL;
  loaded = false;
  load(dump,0);
}

Hamiltonian::Hamiltonian(const Hamiltonian &source, const bool forceSpinUnrestricted, const bool oneElectron, const bool twoElectron)
  : OrbitalSpace(source)
  ,loaded(source.loaded)
  , coreEnergy(source.coreEnergy)
  , basisSize(source.basisSize), ijSize(source.ijSize), ijklSize(source.ijklSize)
{
  xout <<"Hamiltonian copy constructor source="<<&source<<", this="<<this<<std::endl;
  if (forceSpinUnrestricted) spinUnrestricted = true;
  bracket_integrals_a = bracket_integrals_b = NULL;
  if (loaded) {
       integrals_a = oneElectron ? new std::vector<double>(*source.integrals_a) : NULL;
       bracket_integrals_a = (oneElectron && source.bracket_integrals_a != NULL) ? new std::vector<double>(*source.bracket_integrals_a) : NULL;
     if (source.integrals_aa != NULL || spinUnrestricted) {
       integrals_aa = twoElectron ? new std::vector<double>(*source.integrals_aa) : NULL;
       bracket_integrals_aa = twoElectron ? new std::vector<double>(*source.bracket_integrals_aa) : NULL;
       bracket_integrals_ab = (twoElectron && source.bracket_integrals_ab != NULL)
           ?  new std::vector<double>(*source.bracket_integrals_ab)
           : NULL;
     }
     if (spinUnrestricted) {
       integrals_b = oneElectron ? new std::vector<double>(*source.integrals_b): NULL;
       if (source.integrals_ab != NULL || spinUnrestricted) {
         integrals_ab = (twoElectron && source.integrals_ab != NULL)
             ? new std::vector<double>(*source.integrals_ab)
             : NULL;
         integrals_bb = twoElectron ? new std::vector<double>(*source.integrals_bb) : NULL;
         bracket_integrals_bb = twoElectron ? new std::vector<double>(*source.bracket_integrals_bb) : NULL;
         bracket_integrals_b = (oneElectron && source.bracket_integrals_b != NULL) ? new std::vector<double>(*source.bracket_integrals_b) : NULL;
       }
     } else {
       integrals_b = integrals_a;
       integrals_ab = integrals_aa;
       integrals_bb = integrals_aa;
       bracket_integrals_bb = bracket_integrals_aa;
       bracket_integrals_b = bracket_integrals_a;
     }
  }
//  xout << "Hamiltonian copy constructor, old integrals_a="<<&source.integrals_a[0]<<", new integrals_a ="<<&integrals_a[0]<<std::endl;
//  xout << "Hamiltonian copy constructor, old integrals_b="<<&source.integrals_b[0]<<", new integrals_b ="<<&integrals_b[0]<<std::endl;
}

Hamiltonian::~Hamiltonian() {
  xout <<"Hamiltonian deconstructor this="<<this<<std::endl;
  deconstructBraKet();
}

void Hamiltonian::load(std::string filename, int verbosity) {
  FCIdump d(filename);
  load(&d, verbosity);
}

void Hamiltonian::load(FCIdump* dump, int verbosity) {
  xout <<"Hamiltonian::load, this="<<this<<std::endl;
  profiler.start("Hamiltonian::load");
  if (loaded) unload();
  if (verbosity) xout <<"Load hamiltonian from " << dump->fileName() <<std::endl;
  //    State::load(filename);

  basisSize = dump->parameter("NORB").at(0);

  ijSize = total(0,1);
  ijklSize = pairSpace[1].total(0);
//  xout << "ijklSize=" << ijklSize <<std::endl;
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


  {
    dump->rewind();
//  std::ifstream s;
//  s.open(dump->fileName().c_str());
//  std::string ss;
  double value;
  FCIdump::integralType type;
  int i,j,k,l,ij,kl;
  while ((type=dump->nextIntegral(i,j,k,l,value))!=FCIdump::endOfFile) {
    ij = i > j ? (i*(i-1))/2+j : (j*(j-1))/2+i;
    kl = k > l ? (k*(k-1))/2+l : (l*(l-1))/2+k;
    if (type == FCIdump::I2aa) {
      if (verbosity>2) xout << "aa("<< i << j <<"|"<< k << l <<") [" << int2Index(i,j,k,l) << "]= " << value <<std::endl;
      if (verbosity>2) xout << "aa("<< k << l <<"|"<< i << j <<") [" << int2Index(k,l,i,j) << "]= " << value <<std::endl;
      integrals_aa->at(int2Index(i,j,k,l))=value;
      integrals_aa->at(int2Index(k,l,i,j))=value;
    } else if (type == FCIdump::I2ab) {
      if (verbosity>2) xout << "ab("<< i << j <<"|"<< k << l <<") [" << int2Index(i,j,k,l) << "]= " << value <<std::endl;
      integrals_ab->at(int2Index(i,j,k,l))=value;
    } else if (type == FCIdump::I2bb) {
      if (verbosity>2) xout << "bb("<< i << j <<"|"<< k << l <<") [" << int2Index(i,j,k,l) << "]= " << value <<std::endl;
      if (verbosity>2) xout << "bb("<< k << l <<"|"<< i << j <<") [" << int2Index(k,l,i,j) << "]= " << value <<std::endl;
      integrals_bb->at(int2Index(i,j,k,l))=value;
      integrals_bb->at(int2Index(k,l,i,j))=value;
    } else if (type == FCIdump::I1a) {
      if (verbosity>1) xout << "ha("<< i <<","<< j <<") = " << value <<std::endl;
      integrals_a->at(int1Index(i,j))=value;
    } else if (type == FCIdump::I1b) {
      if (verbosity>1) xout << "hb("<< i <<","<< j <<") = " << value <<std::endl;
      integrals_b->at(int1Index(i,j))=value;
    } else if (type == FCIdump::I0)
      coreEnergy = value;
  }
  }
  loaded=true;
  xout <<"load sets loaded=true this="<<this<<std::endl;
  if (verbosity>3) {
    xout << "integrals_a: ";copy(integrals_a->begin(), integrals_a->end(), std::ostream_iterator<double>(xout, ", "));xout <<std::endl;
    xout << "integrals_aa: ";copy(integrals_aa->begin(), integrals_aa->end(), std::ostream_iterator<double>(xout, ", "));xout <<std::endl;
  }
  constructBraKet();
  //    xout <<str(3) <<std::endl;exit(0);
  profiler.stop("Hamiltonian::load");
}

#define del(x) if (x != NULL && x->size()) delete x; x=NULL;
#define del2(x,y) if (x != y) del(x); del(y); del(x);
void Hamiltonian::constructBraKet(int neleca, int nelecb)
{
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
        if (nelecb == 0 && integrals_a != NULL) bracket_integrals_a = new std::vector<double>(*integrals_a);
        del(bracket_integrals_b);
        if (neleca != 0 && integrals_b != NULL) bracket_integrals_b = new std::vector<double>(*integrals_b);
        xout <<"constructBraKet bracket_integrals_{a,b}="<<bracket_integrals_a<<", "<<bracket_integrals_b<<std::endl;
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


void Hamiltonian::deconstructBraKet()
{
  xout << "deconstructBraKet for object at "<<this<<", loaded="<<loaded<<std::endl;
  if (! loaded) return;
  xout << "integrals_a "<<integrals_a<<std::endl;
  xout << "bracket_integrals_a "<<bracket_integrals_a<<std::endl;
  if (bracket_integrals_a != NULL) xout << "bracket_integrals_a->size() "<<bracket_integrals_a->size()<<std::endl;
  xout << "bracket_integrals_b "<<bracket_integrals_b<<std::endl;
  if (bracket_integrals_b != NULL) xout << "bracket_integrals_b->size() "<<bracket_integrals_b->size()<<std::endl;
  del2(bracket_integrals_b, bracket_integrals_a);
  xout <<"del2 done"<<std::endl;
  del2(bracket_integrals_aa, bracket_integrals_bb);
  del(bracket_integrals_ab);
}

void Hamiltonian::unload() {
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

std::string Hamiltonian::str(int verbosity) const
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

size_t Hamiltonian::int1Index(unsigned int i, unsigned int j) const {
  return pairIndex(i,j,1);
}

size_t Hamiltonian::int2Index(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
{
  return quadIndex(i,j,k,l,1,0);
}

std::vector<double> Hamiltonian::int1(int spin) const
{
  std::vector<double> result(basisSize,(double)0);
  std::vector<double> * integrals = spin < 0 ? integrals_b : integrals_a;
  for (unsigned int i=0; i < basisSize; i++) {
    result[i] = integrals->at(int1Index(i+1,i+1));
  }
  return result;
}


std::vector<double> Hamiltonian::intJ(int spini, int spinj) const
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

std::vector<double> Hamiltonian::intK(int spin) const
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

Hamiltonian Hamiltonian::FockHamiltonian(const Determinant &reference) const
{
  Hamiltonian f;
  for (int i=0; i<8; i++)
    f[i]=at(i);
  f.calculateOffsets();
  // xout << "FockHamiltonian Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout << "FockHamiltonian Reference beta: "<<reference.stringBeta<<std::endl;
  bool closed = reference.stringAlpha==reference.stringBeta;
//  xout << "FockHamiltonian Reference alpha=beta: "<<closed<<std::endl;
  f.spinUnrestricted = spinUnrestricted || ! closed;
//  xout << "FockHamiltonian spinUnrestricted="<<f.spinUnrestricted<<std::endl;
  std::vector<unsigned int> refAlphaOrbitals=reference.stringAlpha.orbitals();
  std::vector<unsigned int> refBetaOrbitals=reference.stringBeta.orbitals();
  f.coreEnergy = coreEnergy;
  f.basisSize = basisSize;
  f.ijklSize = ijklSize;
  f.ijSize = ijSize;
  f.orbital_symmetries = orbital_symmetries;
  f.integrals_a = new std::vector<double>(ijSize,0.0);
  *f.integrals_a = *integrals_a;
  // xout <<"reference.stringAlpha.orbitals ";for (size_t i=0; i < reference.stringAlpha.orbitals().size(); i++) xout <<reference.stringAlpha.orbitals()[i]<<" ";xout <<std::endl;
  // for (std::vector<unsigned int>::const_iterator o=reference.stringAlpha.orbitals().begin(); o != reference.stringAlpha.orbitals().end(); o++)
  for (std::vector<unsigned int>::const_iterator o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
  {
    // xout << "FockHamiltonian Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout<< "f alpha, alpha occ: " <<*o << std::endl;
    for (unsigned int i=1; i<=basisSize; i++)
      for (unsigned int j=1; j<=i; j++) {
        if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
        (*f.integrals_a)[int1Index(i,j)] += (*integrals_aa)[int2Index(i,j,*o,*o)] - (*integrals_aa)[int2Index(i,*o,*o,j)];
      }
  }
  for (std::vector<unsigned int>::const_iterator o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
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
  xout << "in FockHamiltonian f.bracket_integrals_a created at "<<&f.bracket_integrals_a <<std::endl;
  // xout << "in FockHamiltonian, after alpha f="; for (size_t ij=0; ij< f.integrals_a->size(); ij++) xout <<" "<<(*f.integrals_a)[ij]; xout <<std::endl;
  if (f.spinUnrestricted) {
    f.integrals_b = new std::vector<double>(ijSize,0.0);
    *f.integrals_b = *integrals_b;
    for (std::vector<unsigned int>::const_iterator o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
    {
      // xout<< "f beta, beta occ: " <<*o << std::endl;
      for (unsigned int i=1; i<=basisSize; i++)
        for (unsigned int j=1; j<=i; j++) {
          if (orbital_symmetries[i-1]!=orbital_symmetries[j-1]) continue;
          (*f.integrals_b)[int1Index(i,j)] += (*integrals_bb)[int2Index(i,j,*o,*o)] - (*integrals_bb)[int2Index(i,*o,*o,j)];
        }
    }
    for (std::vector<unsigned int>::const_iterator o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
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
  xout <<"FockHamiltonian at "<<&f<<" sets loaded=true"<<std::endl;
  return f;
}

Hamiltonian Hamiltonian::sameSpinHamiltonian(const Determinant &reference) const
{
  Hamiltonian result = *this;
//  xout << "result when initialized: "<<result.str(2)<<std::endl;
  result.spinUnrestricted = true;
  if (!spinUnrestricted) *(result.integrals_b = new std::vector<double>(integrals_a->size())) = *result.integrals_a;
//  xout << "sameSpinHamiltonian, old integrals_a="<<&integrals_a[0]<<", new integrals_a ="<<&result.integrals_a[0]<<std::endl;
//  xout << "sameSpinHamiltonian, old integrals_b="<<&integrals_b[0]<<", new integrals_b ="<<&result.integrals_b[0]<<std::endl;
  Determinant ra = reference; ra.stringBeta.nullify();
//  xout << "this before alpha fock: "<<str(2)<<std::endl;
//  xout << "result before alpha fock: "<<result.str(2)<<std::endl;
  Hamiltonian f = this->FockHamiltonian(ra);
  xout << "sameSpinHamiltonian, after alpha Fock, f at "<<&f<<", f.bracket_integrals_a at "<<f.bracket_integrals_a<<std::endl;
//  xout << "this after alpha fock: "<<str(2)<<std::endl;
//  xout << "result before alpha: "<<result.str(2)<<std::endl;
  for (size_t i=0; i<integrals_a->size(); i++) {
    result.integrals_a->at(i) = this->integrals_a->at(i) - f.integrals_a->at(i);
//    xout << " result a = " << result.integrals_a->at(i) <<"=" << this->integrals_a->at(i) <<"-"<< f.integrals_a->at(i)<<std::endl;
  }
//  xout << "this after alpha: "<<str(2)<<std::endl;
//  xout << "result after alpha: "<<result.str(2)<<std::endl;
  ra = reference; ra.stringAlpha.nullify();
  f = this->FockHamiltonian(ra);
  xout << "sameSpinHamiltonian, after beta Fock, f at "<<&f<<", f.bracket_integrals_a at "<<f.bracket_integrals_a<<std::endl;
  xout << "f.bracket_integrals_a->size() "<<f.bracket_integrals_a->size()<<std::endl;
//  xout << "sameSpinHamiltonian, fock integrals_a="<<&f.integrals_a[0]<<", fock integrals_b ="<<&f.integrals_b[0]<<std::endl;
  for (size_t i=0; i<integrals_b->size(); i++) {
    result.integrals_b->at(i) = this->integrals_b->at(i) - f.integrals_b->at(i);
//    xout << " result b = " << result.integrals_b->at(i) <<"=" << this->integrals_b->at(i) <<"-"<< f.integrals_b->at(i)<<std::endl;
}
  xout << "f.bracket_integrals_a->size() "<<f.bracket_integrals_a->size()<<std::endl;
  xout << "result.bracket_integrals_ab "<<result.bracket_integrals_ab<<std::endl;
  if (result.bracket_integrals_ab != NULL) {
//    xout << "result.bracket_integrals_ab isn't null; size="<<result.bracket_integrals_ab->size()<<std::endl;
    // not very satisfactory
    result.bracket_integrals_ab->assign(result.bracket_integrals_ab->size(),(double)0);
//    delete result.bracket_integrals_ab;
//    result.bracket_integrals_ab = NULL;
  }
  xout << "f.bracket_integrals_a->size() "<<f.bracket_integrals_a->size()<<std::endl;
  if (spinUnrestricted) delete result.integrals_ab; result.integrals_ab = NULL;
//  xout << "result on return: "<<result.str(2)<<std::endl;
  xout << "f.bracket_integrals_a->size() "<<f.bracket_integrals_a->size()<<std::endl;
  return result;
}

#include <assert.h>
Hamiltonian& Hamiltonian::plusminusOperator(const Hamiltonian &other, const char operation)
{
//  if (! compatible(other)) throw "attempt to add incompatible Hamiltonian objects";
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
Hamiltonian& Hamiltonian::operator+=(const Hamiltonian &other)
{
  return plusminusOperator(other,'+');
}
Hamiltonian& Hamiltonian::operator-=(const Hamiltonian &other)
{
  return plusminusOperator(other,'-');
}

#include <cmath>
void Hamiltonian::plusminusEqualsHelper(std::vector<double> *&me,
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
void Hamiltonian::starEqualsHelper(std::vector<double> *&me,
                                    const double factor)
{
  if (factor == (double)1) return;
  if (me == NULL) return;
  size_t n = me->size();
  for (size_t i=0; i<n; i++)
    me->at(i) *= factor ;
}

Hamiltonian& Hamiltonian::operator*=(const double factor)
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

Hamiltonian gci::operator+(const Hamiltonian &h1, const Hamiltonian &h2)
{
  Hamiltonian result = h1;
  return result += h2;
}
Hamiltonian gci::operator-(const Hamiltonian &h1, const Hamiltonian &h2)
{
  Hamiltonian result = h1;
  return result -= h2;
}
Hamiltonian gci::operator*(const Hamiltonian &h1, const double factor)
{
  Hamiltonian result = h1;
  return result *= factor;
}
