#include "gciWavefunction.h"
#include <sstream>
#include <iostream>
#include "gciStringSet.h"
#include "gciTransitionDensity.h"
//#include "mkl.h"

Wavefunction::Wavefunction(FCIdump *dump) : State(dump) {
  buildStrings();
}

Wavefunction::Wavefunction(std::string filename) : State(filename) {
  if (filename!="") buildStrings();
}
Wavefunction::Wavefunction(OrbitalSpace* h, int n, int s, int m2) : State(h,n,s,m2) {
  buildStrings();
}

Wavefunction::Wavefunction(const Wavefunction &other) : State(other)
{
  alphaStrings.resize(8); betaStrings.resize(8);
  for (int i=0;i<8;i++)
  {
    alphaStrings[i] = other.alphaStrings[i];
    betaStrings[i] = other.betaStrings[i];
  }
  dimension = other.dimension;
  _blockOffset = other._blockOffset;
  buffer = other.buffer;
}

void Wavefunction::buildStrings()
{
  alphaStrings.resize(8); betaStrings.resize(8);
  dimension = 0;
  _blockOffset.resize(8);
  for (unsigned int syma=0; syma<8; syma++) {
    unsigned int symb = syma ^ symmetry;
    String stringa(this,1);
    stringa.first((nelec+ms2)/2);
    alphaStrings[syma] = StringSet(stringa, true, syma);
    String stringb(this,-1);
    stringb.first((nelec-ms2)/2);
    betaStrings[symb] = StringSet(stringb, true, symb);
    _blockOffset[syma]=dimension;
    dimension += alphaStrings[syma].size()*betaStrings[symb].size();
  }
}

size_t Wavefunction::size()
{
  return dimension;
}

void Wavefunction::allocate_buffer()
{
  buffer.resize(dimension,(double)0);
}

void Wavefunction::set(size_t offset, const double val)
{
  buffer.at(offset) = val;
}

void Wavefunction::set(const double value)
{
  allocate_buffer();
  for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b=value;
}

void Wavefunction::diagonalHamiltonian(Hamiltonian &hamiltonian)
{
  std::vector<double> ha=hamiltonian.int1(1);
  std::vector<double> hbb=hamiltonian.int1(-1);
  std::vector<double> Jaa=hamiltonian.intJ(1,1);
  std::vector<double> Jab=hamiltonian.intJ(1,-1);
  std::vector<double> Jbb=hamiltonian.intJ(1,1);
  std::vector<double> Kaa=hamiltonian.intK(1);
  std::vector<double> Kbb=hamiltonian.intK(-1);
  //    xout << "ha" <<std::endl;
  //        for (size_t i=0; i<hamiltonian->basisSize; i++)
  //            xout << ha[i] << " ";
  //        xout <<std::endl;
  //    xout << "Jaa" <<std::endl;
  //    for (size_t j=0; j<hamiltonian->basisSize; j++) {
  //        for (size_t i=0; i<hamiltonian->basisSize; i++)
  //            xout << Jaa[i+j*hamiltonian->basisSize] << " ";
  //        xout <<std::endl;
  //    }
  //    xout << "Kaa" <<std::endl;
  //    for (size_t j=0; j<hamiltonian->basisSize; j++) {
  //        for (size_t i=0; i<hamiltonian->basisSize; i++)
  //            xout << Kaa[i+j*hamiltonian->basisSize] << " ";
  //        xout <<std::endl;
  //    }
  size_t offset=0;
  set(hamiltonian.coreEnergy);
  for (unsigned int syma=0; syma<8; syma++) {
    unsigned int symb = syma ^ symmetry;
    size_t nsa = alphaStrings[syma].size();
    size_t nsb = betaStrings[symb].size();
    size_t nact = orbitalSpace->total();
    if (! nsa || ! nsb) continue;
    std::vector<double> ona = alphaStrings[syma].occupationNumbers();
    std::vector<double> onb = betaStrings[symb].occupationNumbers();
    if (false && orbitalSpace->spinUnrestricted) { // UHF
    } else { // RHF
      for (size_t ia=0; ia < nsa; ia++) {
        std::vector<double> on(onb);
        for (size_t i=0; i<nact; i++) {
          for (size_t ib=0; ib < nsb; ib++)
            on[ib+i*nsb] += ona[ia+i*nsa];
        }
        //                const double one=(double) 1;
        //                const MKL_INT ione = 1;
        //                const MKL_INT m = betaStrings[symb].size();
        //                const MKL_INT n = hamiltonian->basisSize;
        //                DGEMV("N",&m,&n, &one, &on[0], &m, &ha[0], &ione, &one, &buffer[offset+ia*m], &ione);
        for (size_t i=0; i<nact; i++)
          for (size_t ib=0; ib < nsb; ib++)
            buffer[offset+ib] += on[ib+i*nsb] * ha[i];
        for (size_t i=0; i<nact; i++) {
          for (size_t j=0; j<=i; j++) {
            double zz = Jaa[j+i*nact] - (double)0.5 * Kaa[j+i*nact];
            if (i == j) zz *= (double)0.5;
            for (size_t ib=0; ib < nsb; ib++)
              buffer[offset+ib] += on[ib+i*nsb] * on[ib+j*nsb] * zz;
          }
        }
        double vv = (double) ms2*ms2;
        std::vector<double> f(nsb,(double)0);
        for (size_t i=0; i<nact; i++)
          for (size_t ib=0; ib < nsb; ib++) {
            on[ib+i*nsb] *= ((double)2 - on[ib+i*nsb]); // on becomes a mask for singly-occupied orbitals
            f[ib] += on[ib+i*nsb];
          }
        for (size_t ib=0; ib < nsb; ib++)
          f[ib] = f[ib] < (double) 1.1 ? (double) 1.1 : f[ib]; // mask off singularities (this code does nothing for 0 or 1 open-shell orbitals)
        for (size_t ib=0; ib < nsb; ib++)
          f[ib] = (vv-f[ib]) / (f[ib]*(f[ib]-(double)1));
        for (size_t i=0; i<nact; i++) {
          for (size_t j=0; j<i; j++) {
            double zz = -(double)0.5 * Kaa[i*nact+j];
            for (size_t ib=0; ib < nsb; ib++)
              buffer[offset+ib] += f[ib] * on[ib+i*nsb] * on[ib+j*nsb] * zz;
          }
          double zz = -(double)0.25 * Kaa[i*nact+i];
          for (size_t ib=0; ib < nsb; ib++)
            buffer[offset+ib] += on[ib+i*nsb] * zz;
        }
        offset += nsb;
      }
    }
  }
}

Wavefunction& Wavefunction::operator*=(const double &value)
{
  for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b*=value;
  return *this;
}

//Wavefunction& Wavefunction::operator = ( const Wavefunction &other)
//{
//    if (this != &other) {
//    if (! compatible(other)) throw "attempt to copy between incompatible Wavefunction objects";
//        *this = other;
//    }
//    return *this;
//}

Wavefunction& Wavefunction::operator+=(const Wavefunction &other)
{
  if (! compatible(other)) throw "attempt to add incompatible Wavefunction objects";
  xout << "Wavefunction::operator += &this=" << this <<" &other="<<&other <<std::endl;
  for (size_t i=0; i<buffer.size(); i++)  buffer[i] += other.buffer[i];
  return *this;
}


Wavefunction& Wavefunction::operator-=(const Wavefunction &other)
{
  if (! compatible(other)) throw "attempt to add incompatible Wavefunction objects";
  for (size_t i=0; i<buffer.size(); i++)  buffer[i] -= other.buffer[i];
  return *this;
}


Wavefunction gci::operator+(const Wavefunction &w1, const Wavefunction &w2)
{
  Wavefunction result = w1;
  return result += w2;
}

Wavefunction gci::operator-(const Wavefunction &w1, const Wavefunction &w2)
{
  Wavefunction result = w1;
  return result -= w2;
}

Wavefunction gci::operator*(const Wavefunction &w1, const double &value)
{
  Wavefunction result = w1;
  return result *= value;
}

Wavefunction gci::operator*(const double &value, const Wavefunction &w1)
{
  Wavefunction result = w1;
  return result *= value;
}

double gci::operator *(const Wavefunction &w1, const Wavefunction &w2)
{
  if (! w1.compatible(w2)) throw "attempt to form scalar product between incompatible Wavefunction objects";
  double result=(double)0;
  for (size_t i=0; i<w1.buffer.size(); i++) result += w1.buffer[i]*w2.buffer[i];
  return result;
}

std::string Wavefunction::str(int verbosity) const
{
  std::ostringstream s;
  if (verbosity >= 1) {
    s<<std::endl<<"Wavefunction object at address " << this ;
    s<<std::endl<<"Values at address " << &buffer <<" and of length " << buffer.size();
    size_t address=0;
    for (unsigned int syma=0; syma<8; syma++) {
      unsigned int symb = syma ^ symmetry ;
      if (alphaStrings[syma].size() && betaStrings[symb].size()) {
        s<<std::endl<< "Alpha strings of symmetry "<<syma+1<<":";
        for (StringSet::const_iterator i=alphaStrings[syma].begin(); i!=alphaStrings[syma].end(); i++) s <<std::endl<< i->str();
        s<<std::endl<< "Beta strings of symmetry "<<symb+1<<":";
        for (StringSet::const_iterator i=betaStrings[symb].begin(); i!=betaStrings[symb].end(); i++) s <<std::endl<< i->str();
        if (buffer.size() == dimension && verbosity >=2) {
          s<<std::endl<<"Values:";
          for (size_t i=0; i<alphaStrings[syma].size(); i++) {
            s<<std::endl;
            for (size_t j=0; j<betaStrings[symb].size(); j++) {
              s << buffer[address++] << " ";
            }
          }
        }
      }
    }
  }
  return this->State::str(verbosity)+s.str();
}

bool Wavefunction::compatible(const Wavefunction &other) const
{
  return dimension==other.dimension && buffer.size() == other.buffer.size();
}

size_t Wavefunction::minloc()
{
  size_t result=0;
  for (size_t offset=0; offset<buffer.size(); offset++)
    if (buffer[offset] < buffer[result]) result=offset;
  return result;
}

size_t Wavefunction::maxloc()
{
  size_t result=0;
  for (size_t offset=0; offset<buffer.size(); offset++)
    if (buffer[offset] > buffer[result]) result=offset;
  return result;
}

double Wavefunction::at(size_t offset)
{
  return buffer.at(offset);
}

Determinant Wavefunction::determinantAt(size_t offset)
{
  size_t address=0;
  for (unsigned int syma=0; syma<8; syma++) {
    unsigned int symb = syma ^ symmetry ;
    size_t newaddress = address +  alphaStrings[syma].size() * betaStrings[symb].size();
    if (offset >= address && offset < newaddress) {
      size_t a=(offset-address)/betaStrings[symb].size();
      size_t b=offset-address-a*betaStrings[symb].size();
      return Determinant(this,&alphaStrings[syma][a],&betaStrings[symb][b]);
    }
    address=newaddress;
  }
  throw "Wavefunction::determinantAt cannot find";
}

size_t Wavefunction::blockOffset(const unsigned int syma) const
{
  return _blockOffset.at(syma);
}

void Wavefunction::hamiltonianOnWavefunction(Hamiltonian &h, const Wavefunction &w)
{
  for (size_t i=0; i<buffer.size(); i++)
    buffer[i] = h.coreEnergy * w.buffer[i];


  size_t offset=0;
  for (unsigned int syma=0; syma<8; syma++) {
    unsigned int symb = w.symmetry^syma;
    size_t nsa = alphaStrings[syma].size();
    size_t nsb = betaStrings[symb].size();
    xout << "syma="<<syma<<", symb="<<symb<<std::endl;
    TransitionDensity d(w,
                        w.alphaStrings[syma].begin(),
                        w.alphaStrings[syma].end(),
                        w.betaStrings[symb].begin(),
                        w.betaStrings[symb].end(),
                        1,true, !h.spinUnrestricted);
    xout <<"Transition density: "<<d<<std::endl;
    for (size_t iab=0; iab<nsa*nsb; iab++)
      for (size_t ij=0; ij<w.orbitalSpace->total(0,1); ij++)
        buffer[offset+iab]+=d[ij*nsa*nsb+iab] * (*h.integrals_a)[ij];
    if (h.spinUnrestricted) {
    TransitionDensity d(w,
                        w.alphaStrings[syma].begin(),
                        w.alphaStrings[syma].end(),
                        w.betaStrings[symb].begin(),
                        w.betaStrings[symb].end(),
                        1,false, true);
    xout <<"Transition density: "<<d<<std::endl;
    for (size_t iab=0; iab<nsa*nsb; iab++)
      for (size_t ij=0; ij<w.orbitalSpace->total(0,1); ij++)
        buffer[offset+iab]+=d[ij*nsa*nsb+iab] * (*h.integrals_b)[ij];
    }
//    size_t offa = offset;
//    for (StringSet::iterator s = alphaStrings[syma].begin(); s != alphaStrings[syma].end(); s++) {
//      //            xout << "alpha string "<<*s<<std::endl;
//      ExcitationSet ee(*s,alphaStrings[syma],1,1);
//      //            xout << "alpha excitations " << ee.str() <<std::endl;
//      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
//        //                xout << "alpha excitation " << e->orbitalAddress <<"="<<(*h.integrals_a)[e->orbitalAddress]<<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
//        for (size_t ib=0; ib<nsb; ib++)
//          buffer[offa+ib] += (*h.integrals_a)[e->orbitalAddress] * e->phase * w.buffer[offset+e->stringIndex*nsb+ib];
//      }
//      offa += nsb;
//    }
//    size_t offb = offset;
//    for (StringSet::iterator s = betaStrings[symb].begin(); s != betaStrings[symb].end(); s++) {
//      ExcitationSet ee(*s,betaStrings[symb],1,1);
//      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
//        //                xout << "beta excitation " << e->orbitalAddress <<"="<<(*h.integrals_b)[e->orbitalAddress]<<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
//        for (size_t ia=0; ia<nsa; ia++)
//          buffer[offb+ia*nsb] += (*h.integrals_b)[e->orbitalAddress] * e->phase * w.buffer[offset+e->stringIndex+ia*nsb];
//      }
//      offb ++;
//    }
    offset += nsa*nsb;
  }

  xout <<"residual after 1-electron:"<<std::endl<<str(2)<<std::endl;

  if (h.bracket_integrals_aa != NULL) { // two-electron contribution, alpha-alpha
    size_t nsbbMax = 64; // temporary static
    for (unsigned int syma=0; syma<8; syma++) {
      StringSet aa(w.alphaStrings,2,0,syma);
      if (aa.size()==0) continue;
      xout <<"StringSet aa: " <<aa.str(2)<<std::endl;
      for (unsigned int symb=0; symb<8; symb++) {
        unsigned int symexc = syma^symb^w.symmetry;
        size_t nexc = h.pairSpace[-1][symexc];
        size_t nsb = betaStrings[symb].size(); if (nsb==0) continue;
        for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsbbMax > aa.end() ? aa.end() : aa0+nsbbMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
          size_t nsa = aa1-aa0;
          xout <<"nsa= "<<nsa<<std::endl;
          TransitionDensity d(w,aa0,aa1,w.betaStrings[symb].begin(),w.betaStrings[symb].end(),-1,false,false);
          xout <<"Transition density aa: "<<d<<std::endl;
          TransitionDensity e(d); e.assign(d.size(),(double)0);
          xout << "nexc="<<nexc<<", d.size()="<<d.size()<<", nsa="<<nsa<<", nsb="<<nsb<<std::endl;
          xout << "h.pairSpace[-1].at(symexc)"<<h.pairSpace[-1][symexc]<<std::endl;
          if (nexc * nsa * nsb != d.size()) throw "nexc";
          for (size_t excd=0; excd<nexc; excd++)
            for (size_t exce=0; exce<nexc; exce++)
              for (size_t ab=0; ab < nsa*nsb; ab++)
                e[ab+exce*nsa*nsb] += d[ab+excd*nsa*nsb]
                    * (*h.bracket_integrals_aa)[h.pairSpace[-1].offset(0,symexc,0)+excd*nexc+exce];
          e.action(*this);
          xout <<"residual after aa:"<<std::endl<<str(2)<<std::endl;
        }
      }
    }
  }

  if (h.bracket_integrals_bb != NULL) { // two-electron contribution, beta-beta
    size_t nsbbMax = 64; // temporary static
    for (unsigned int symb=0; symb<8; symb++) {
      StringSet bb(w.betaStrings,2,0,symb);
      if (bb.size()==0) continue;
      xout <<"StringSet bb: " <<bb.str(2)<<std::endl;
      for (unsigned int syma=0; syma<8; syma++) {
        unsigned int symexc = symb^syma^w.symmetry;
        size_t nexc = h.pairSpace[-1][symexc];
        size_t nsa = alphaStrings[syma].size(); if (nsa==0) continue;
        for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
          size_t nsb = bb1-bb0;
          xout <<"nsb= "<<nsb<<std::endl;
          TransitionDensity d(w,w.alphaStrings[syma].begin(),w.alphaStrings[syma].end(),bb0,bb1,-1,false,false);
          xout <<"Transition density bb: "<<d<<std::endl;
          TransitionDensity e(d); e.assign(d.size(),(double)0);
          xout << "nexc="<<nexc<<", d.size()="<<d.size()<<", nsb="<<nsb<<", nsa="<<nsa<<std::endl;
          xout << "h.pairSpace[-1].at(symexc)"<<h.pairSpace[-1][symexc]<<std::endl;
          if (nexc * nsb * nsa != d.size()) throw "nexc";
          for (size_t excd=0; excd<nexc; excd++)
            for (size_t exce=0; exce<nexc; exce++)
              for (size_t ab=0; ab < nsb*nsa; ab++)
                e[ab+exce*nsb*nsa] += d[ab+excd*nsb*nsa]
                    * (*h.bracket_integrals_bb)[h.pairSpace[-1].offset(0,symexc,0)+excd*nexc+exce];
          e.action(*this);
          xout <<"residual after bb:"<<std::endl<<str(2)<<std::endl;
        }
      }
    }
  }

  if (h.bracket_integrals_ab != NULL) { // two-electron contribution, alpha-beta
    size_t nsaaMax = 16; // temporary static
    size_t nsbbMax = 16; // temporary static
    for (unsigned int symb=0; symb<8; symb++) {
      StringSet bb(w.betaStrings,1,0,symb);
      if (bb.size()==0) continue;
      xout <<"StringSet bb: " <<bb.str(2)<<std::endl;
      for (unsigned int syma=0; syma<8; syma++) {
      StringSet aa(w.alphaStrings,1,0,syma);
      if (aa.size()==0) continue;
        unsigned int symexc = symb^syma^w.symmetry;
        size_t nexc = h.pairSpace[0][symexc];
        size_t nsa = alphaStrings[syma].size(); if (nsa==0) continue;
        for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsaaMax > aa.end() ? aa.end() : aa0+nsaaMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
          size_t nsa = aa1-aa0;
          for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
            size_t nsb = bb1-bb0;
            xout <<"nsb= "<<nsb<<std::endl;
            TransitionDensity d(w,aa0,aa1, bb0,bb1,0,false,false);
            xout <<"Transition density ab: "<<d<<std::endl;
            TransitionDensity e(d); e.assign(d.size(),(double)0);
            xout << "nexc="<<nexc<<", d.size()="<<d.size()<<", nsb="<<nsb<<", nsa="<<nsa<<std::endl;
            xout << "h.pairSpace[-1].at(symexc)"<<h.pairSpace[-1][symexc]<<std::endl;
            if (nexc * nsb * nsa != d.size()) throw "nexc";
            for (size_t excd=0; excd<nexc; excd++)
              for (size_t exce=0; exce<nexc; exce++)
                for (size_t ab=0; ab < nsb*nsa; ab++)
                  e[ab+exce*nsb*nsa] += d[ab+excd*nsb*nsa]
                      * (*h.bracket_integrals_ab)[h.pairSpace[0].offset(0,symexc,0)+excd*nexc+exce];
            xout <<"E matrix ab: "<<e<<std::endl;
            e.action(*this);
          xout <<"residual after ab:"<<std::endl<<str(2)<<std::endl;
          }
        }
      }
    }
  }

}
