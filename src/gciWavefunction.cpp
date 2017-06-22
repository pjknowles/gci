#include "gci.h"
#include "gciWavefunction.h"
#include <sstream>
#include <iostream>
#include "gciMolpro.h"
#include "gciStringSet.h"
#include "gciTransitionDensity.h"
#include "Profiler.h"
#include <algorithm>

Wavefunction::Wavefunction(const FCIdump &dump) : State(dump) {
  distributed = false;
  buildStrings();
}

Wavefunction::Wavefunction(std::string filename) : State(filename) {
  distributed = false;
  if (filename!="") buildStrings();
}
Wavefunction::Wavefunction(OrbitalSpace* h, int n, int s, int m2) : State(h,n,s,m2) {
  distributed = false;
  buildStrings();
}

Wavefunction::Wavefunction(const State& state) : State(state) {
  distributed = false;
  buildStrings();
}

Wavefunction::Wavefunction(const Wavefunction &other) : LinearAlgebra::vector<double>(), State(other)
{
  distributed = false;
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
  profiler->start("buildStrings");
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
  profiler->stop("buildStrings");
}

size_t Wavefunction::size() const
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

void Wavefunction::diagonalOperator(const Operator &op)
{
  profiler->start("diagonalOperator");
  auto ha=op.int1(true);
  auto hbb=op.int1(false);
  auto Jaa=op.intJ(true,true);
  auto Jab=op.intJ(true,false);
  auto Jbb=op.intJ(false,false);
  auto Kaa=op.intK(true);
  auto Kbb=op.intK(false);
  size_t offset=0;
  set(op.m_O0);
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
      // static task distribution
      size_t chunk = (nsa-1)/parallel_size+1;
      for (size_t ia=chunk*parallel_rank; ia < chunk*(parallel_rank+1) && ia < nsa; ia++) {
        std::vector<double> on(onb);
        for (size_t i=0; i<nact; i++) {
          for (size_t ib=0; ib < nsb; ib++)
            on[ib+i*nsb] += ona[ia+i*nsa];
        }
        for (size_t i=0; i<nact; i++)
          for (size_t ib=0; ib < nsb; ib++)
            buffer[offset+ia*nsb+ib] += on[ib+i*nsb] * ha[i];
         if (op.m_rank>1){
          for (size_t i=0; i<nact; i++) {
            for (size_t j=0; j<=i; j++) {
              double zz = Jaa(j,i) - (double)0.5 * Kaa(j,i);
              if (i == j) zz *= (double)0.5;
              for (size_t ib=0; ib < nsb; ib++)
                buffer[offset+ia*nsb+ib] += on[ib+i*nsb] * on[ib+j*nsb] * zz;
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
              double zz = -(double)0.5 * Kaa(i,j);
              for (size_t ib=0; ib < nsb; ib++)
                buffer[offset+ia*nsb+ib] += f[ib] * on[ib+i*nsb] * on[ib+j*nsb] * zz;
            }
            double zz = -(double)0.25 * Kaa(i,i);
            for (size_t ib=0; ib < nsb; ib++)
              buffer[offset+ia*nsb+ib] += on[ib+i*nsb] * zz;
          }
        }
      }
      gather_chunks(&buffer[offset],nsa*nsb,chunk*nsb);
    }
    offset +=nsa*nsb;
  }
//  xout << "diagonal elements"<<std::endl; for (size_t i=0; i < buffer.size(); i++) xout <<" "<<buffer[i]; xout <<std::endl;
  profiler->stop("diagonalOperator");
}

void Wavefunction::axpy(double a, const LinearAlgebra::vector<double> *x)
{
  const Wavefunction* xx=dynamic_cast <const Wavefunction*> (x);
//  xout << "Wavefunction::axpy initial=";
//  for (size_t i=0; i<buffer.size(); i++) xout<<" "<<buffer[i]; xout << std::endl;
//  xout << "Wavefunction::axpy x=";
//  for (size_t i=0; i<buffer.size(); i++) xout<<" "<<xx->buffer[i]; xout << std::endl;
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++) buffer[i] += xx->buffer[i]*a;
  else
    for (size_t i=0; i<buffer.size(); i++) buffer[i] += xx->buffer[i]*a;
//  xout << "Wavefunction::axpy result=";
//  for (size_t i=0; i<buffer.size(); i++) xout<<" "<<buffer[i]; xout << std::endl;
}

void Wavefunction::scal(double a)
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++) buffer[i] *= a;
  else
    for (size_t i=0; i<buffer.size(); i++) buffer[i] *= a;
}

Wavefunction& Wavefunction::operator*=(const double &value)
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (std::vector<double>::iterator b=buffer.begin()+parallel_rank*chunk; b != buffer.end() && b < buffer.begin()+(parallel_rank+1)*chunk; b++) *b*=value;
  else
    for (std::vector<double>::iterator b=buffer.begin(); b != buffer.end(); b++) *b*=value;
  return *this;
}

//Wavefunction& Wavefunction::operator = ( const Wavefunction &other)
//{
//    if (this != &other) {
//    if (! compatible(other)) throw std::domain_error("attempt to copy between incompatible Wavefunction objects");
//        *this = other;
//    }
//    return *this;
//}

Wavefunction& Wavefunction::operator+=(const Wavefunction &other)
{
  if (! compatible(other)) throw std::domain_error("attempt to add incompatible Wavefunction objects");
//  xout << "Wavefunction::operator += &this=" << this <<" &other="<<&other <<std::endl;
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i] += other.buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)  buffer[i] += other.buffer[i];
  return *this;
}


Wavefunction& Wavefunction::operator-=(const Wavefunction &other)
{
  if (! compatible(other)) throw std::domain_error("attempt to add incompatible Wavefunction objects");
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i] -= other.buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] -= other.buffer[i];
  return *this;
}

Wavefunction& Wavefunction::operator-=(const double other)
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i] -= other;
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] -= other;
  return *this;
}


Wavefunction& Wavefunction::operator/=(const Wavefunction &other)
{
  if (! compatible(other)) throw std::domain_error("attempt to add incompatible Wavefunction objects");
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i] /= other.buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] /= other.buffer[i];
  return *this;
}


double Wavefunction::update(const Wavefunction &diagonalH, double & eTruncated, double const dEmax)
{
  if (! compatible(diagonalH)) throw std::domain_error("attempt to combine incompatible Wavefunction objects");
//  xout << "Wavefunction::update this="; for (size_t i=0; i<buffer.size(); i++) xout<<buffer[i]<<" ";xout <<std::endl;
//  xout << "Wavefunction::update diagonalH="; for (size_t i=0; i<diagonalH.buffer.size(); i++) xout<<diagonalH.buffer[i]<<" ";xout <<std::endl;
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  size_t imin = distributed ? parallel_rank*chunk : 0;
  size_t imax = distributed ? (parallel_rank+1*chunk) : buffer.size();
  eTruncated=0;
  double ePredicted=0;
  if (imax > buffer.size()) imax = buffer.size();
  for (size_t i=imin; i<imax; i++) {
    double dE=buffer[i]*buffer[i]/diagonalH.buffer[i];
    if (diagonalH.buffer[i] > 1e-5) ePredicted -= dE;
    if (dE > dEmax)
      buffer[i] = -buffer[i] / diagonalH.buffer[i];
    else {
      buffer[i]=(double)0;
      eTruncated += dE;
    }
  }
  return ePredicted;
}


Wavefunction operator+(const Wavefunction &w1, const Wavefunction &w2)
{
  Wavefunction result = w1;
  return result += w2;
}

Wavefunction operator-(const Wavefunction &w1, const Wavefunction &w2)
{
  Wavefunction result = w1;
  return result -= w2;
}

Wavefunction operator/(const Wavefunction &w1, const Wavefunction &w2)
{
  Wavefunction result = w1;
  return result /= w2;
}

Wavefunction operator*(const Wavefunction &w1, const double &value)
{
  Wavefunction result = w1;
  return result *= value;
}

Wavefunction operator*(const double &value, const Wavefunction &w1)
{
  Wavefunction result = w1;
  return result *= value;
}

double Wavefunction::dot(const LinearAlgebra::vector<double> *other) const
{
 return (*this)*(*(dynamic_cast<const Wavefunction*>(other)));
}

void Wavefunction::times(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b)
{
 const Wavefunction* wa = dynamic_cast<const Wavefunction*>(a);
 const Wavefunction* wb = dynamic_cast<const Wavefunction*>(b);
 size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i]=wa->buffer[i]*wb->buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] = wa->buffer[i]*wb->buffer[i];
}

void Wavefunction::divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b, double shift, bool append, bool negative)
{
  const Wavefunction* wa = dynamic_cast<const Wavefunction*>(a);
  const Wavefunction* wb = dynamic_cast<const Wavefunction*>(b);
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (negative) {
      if (append) {
          if (distributed)
            for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
              buffer[i]-=wa->buffer[i]/(wb->buffer[i]+shift);
          else
            for (size_t i=0; i<buffer.size(); i++)
              buffer[i] -= wa->buffer[i]/(wb->buffer[i]+shift);
        } else {
          if (distributed)
            for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
              buffer[i]=-wa->buffer[i]/(wb->buffer[i]+shift);
          else
            for (size_t i=0; i<buffer.size(); i++)
              buffer[i] =- wa->buffer[i]/(wb->buffer[i]+shift);
        }
    } else {
      if (append) {
          if (distributed)
            for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
              buffer[i]+=wa->buffer[i]/(wb->buffer[i]+shift);
          else
            for (size_t i=0; i<buffer.size(); i++)
              buffer[i] += wa->buffer[i]/(wb->buffer[i]+shift);
        } else {
          if (distributed)
            for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
              buffer[i]=wa->buffer[i]/(wb->buffer[i]+shift);
          else
            for (size_t i=0; i<buffer.size(); i++)
              buffer[i] = wa->buffer[i]/(wb->buffer[i]+shift);
        }
    }
}

void Wavefunction::zero()
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i]=0;
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] = 0;
}

double gci::operator *(const Wavefunction &w1, const Wavefunction &w2)
{
  if (! w1.compatible(w2)) throw std::domain_error("attempt to form scalar product between incompatible Wavefunction objects");
  double result=(double)0;
  size_t chunk = (w1.buffer.size()-1)/parallel_size+1;
  if (w1.distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<w1.buffer.size(); i++)
      result += w1.buffer[i]*w2.buffer[i];
  else
    for (size_t i=0; i<w1.buffer.size(); i++)
      result += w1.buffer[i]*w2.buffer[i];
  return result;
}

Wavefunction& gci::Wavefunction::operator -()
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  if (distributed)
    for (size_t i=parallel_rank*chunk; i<(parallel_rank+1)*chunk && i<buffer.size(); i++)
      buffer[i]=-buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i]=-buffer[i];
  return *this;
}

std::string gci::Wavefunction::str(int verbosity, unsigned int columns) const
{
  std::ostringstream s;
  if (verbosity >= 1) {
    s<<std::endl<<"Wavefunction object at address " << this ;
    s<<std::endl<<"values at address " << &buffer <<" and of length " << buffer.size();
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

size_t Wavefunction::minloc(size_t n)
{
  std::vector<size_t> results;
  for (size_t k=0; k<n; k++) {
      auto m=begin(); while(std::count(results.begin(),results.end(),m-begin())!=0) m++;
      for (auto i=begin(); i!=end(); i++) {
        if (*i < *m && std::count(results.begin(),results.end(),i-begin())==0) m=i;
        }
      if (distributed && parallel_size>1)
        throw std::logic_error("Wavefunction::minloc: parallel implementation unfinished"); //FIXME
      results.push_back(m-begin());
    }
  return results.back();
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
  throw std::runtime_error("Wavefunction::determinantAt cannot find");
}

size_t Wavefunction::blockOffset(const unsigned int syma) const
{
  return _blockOffset.at(syma);
}

#ifdef MOLPRO
#include "gciMolpro.h"
using namespace itf;
#endif

void Wavefunction::operatorOnWavefunction(const Operator &h, const Wavefunction &w, bool parallel_stringset)
{
  auto prof = profiler->push("operatorOnWavefunction");
  if (parallel_rank == 0)
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] += h.m_O0 * w.buffer[i];
  else
    for (size_t i=0; i<buffer.size(); i++)
      buffer[i] = (double)0;
//  xout <<"residual after 0-electron:"<<std::endl<<str(2)<<std::endl;

//  xout <<std::endl<<"w in operatorOnWavefunction="<<w.str(2)<<std::endl;
  DivideTasks(99999999,1,1);

  {
  auto p = profiler->push("1-electron");
  size_t offset=0, nsa=0, nsb=0;
  for (unsigned int syma=0; syma<8; syma++) {
    offset += nsa*nsb;
    unsigned int symb = w.symmetry^syma;
    nsa = alphaStrings[syma].size();
    nsb = betaStrings[symb].size();
    if (!NextTask()) continue;
    {
      profiler->start("1-electron TransitionDensity");
      TransitionDensity d(w,
                          w.alphaStrings[syma].begin(),
                          w.alphaStrings[syma].end(),
                          w.betaStrings[symb].begin(),
                          w.betaStrings[symb].end(),
                          1,true, !h.m_uhf);
      profiler->stop("1-electron TransitionDensity");
      profiler->start("1-electron MXM");
      MxmDrvNN(&buffer[offset],&d[0], &(*h.O1(true).data())[0],
          nsa*nsb,w.orbitalSpace->total(0,1),1,true);
      profiler->stop("1-electron MXM",2*nsa*nsb*w.orbitalSpace->total(0,1));
    }
    if (h.m_uhf) {
      profiler->start("1-electron TransitionDensity");
      TransitionDensity d(w,
                          w.alphaStrings[syma].begin(),
                          w.alphaStrings[syma].end(),
                          w.betaStrings[symb].begin(),
                          w.betaStrings[symb].end(),
                          1,false, true);
      profiler->stop("1-electron TransitionDensity");
      profiler->start("1-electron MXM");
      MxmDrvNN(&buffer[offset],&d[0], &(*h.O1(false).data())[0],
          nsa*nsb,w.orbitalSpace->total(0,1),1,true);
      profiler->stop("1-electron MXM",2*nsa*nsb*w.orbitalSpace->total(0,1));
    }
  }

  }
//  xout <<"residual after 1-electron:"<<std::endl<<str(2)<<std::endl;

  { // two-electron contribution, alpha-alpha
    auto p = profiler->push("aa integrals");
    size_t nsbbMax = 64; // temporary static
    for (unsigned int syma=0; syma<8; syma++) {
        profiler->start("StringSet aa");
        StringSet aa(w.alphaStrings,2,0,syma,parallel_stringset);
        profiler->stop("StringSet aa");
        if (aa.size()==0) continue;
        for (unsigned int symb=0; symb<8; symb++) {
            if (!NextTask()) continue;
            auto praa = profiler->push("aa1 loop");
            unsigned int symexc = syma^symb^w.symmetry;
            //size_t nexc = h.pairSpace.find(-1)->second[symexc];
            size_t nexc = h.O2(true,true,false).block_size(symexc);
            size_t nsb = betaStrings[symb].size(); if (nsb==0) continue;
            for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsbbMax > aa.end() ? aa.end() : aa0+nsbbMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
                size_t nsa = aa1-aa0;
                profiler->start("TransitionDensity aa");
                TransitionDensity d(w,aa0,aa1,w.betaStrings[symb].begin(),w.betaStrings[symb].end(),-1,false,false);
                profiler->stop("TransitionDensity aa");
                TransitionDensity e(d);
                { auto pr1 = profiler->push("MXM aa");
                  MxmDrvNN(&e[0],&d[0], &h.O2(true,true,false).block(symexc)[0], nsa*nsb,nexc,nexc,false);
                }
                { auto pr1 = profiler->push("action aa"); e.action(*this); }
              }
          }
      }
  }
//  xout <<"residual after alpha-alpha on process "<<parallel_rank<<" "<<buffer[0]<<std::endl<<str(2)<<std::endl;

  if (true || h.m_uhf) { // two-electron contribution, beta-beta
      auto p = profiler->push("bb integrals");
      size_t nsbbMax = 64; // temporary static
      for (unsigned int symb=0; symb<8; symb++) {
          StringSet bb(w.betaStrings,2,0,symb,parallel_stringset);
          if (bb.size()==0) continue;
          for (unsigned int syma=0; syma<8; syma++) {
              if (!NextTask()) continue;
              unsigned int symexc = symb^syma^w.symmetry;
              size_t nexc = h.O2(false,false,false).block_size(symexc);
              size_t nsa = alphaStrings[syma].size(); if (nsa==0) continue;
              for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
                  size_t nsb = bb1-bb0;
                  profiler->start("TransitionDensity bb");
                  TransitionDensity d(w,w.alphaStrings[syma].begin(),w.alphaStrings[syma].end(),bb0,bb1,-1,false,false);
                  profiler->stop("TransitionDensity bb");
                  TransitionDensity e(d);
                  { auto pr1 = profiler->push("MXM bb");
                    MxmDrvNN(&e[0],&d[0], &h.O2(false,false,false).block(symexc)[0], nsa*nsb,nexc,nexc,false);
                  }
                  {auto pr1 = profiler->push("action bb"); e.action(*this);}
                }
            }
        }
//      xout <<"residual after beta-beta on process "<<parallel_rank<<" "<<buffer[0]<<std::endl<<str(2)<<std::endl;
    }

  { // two-electron contribution, alpha-beta
    auto p = profiler->push("ab integrals");
    size_t nsaaMax = 640; // temporary static
    size_t nsbbMax = 640; // temporary static
    for (unsigned int symb=0; symb<8; symb++) {
        StringSet bb(w.betaStrings,1,0,symb,parallel_stringset);
        if (bb.size()==0) continue;
        for (unsigned int syma=0; syma<8; syma++) {
            StringSet aa(w.alphaStrings,1,0,syma,parallel_stringset);
            if (aa.size()==0) continue;
            unsigned int symexc = symb^syma^w.symmetry;
            size_t nexc = h.O2(true,false,false).block_size(symexc);
            {
              auto pro = profiler->push("StringSet iterator loops");
              for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsaaMax > aa.end() ? aa.end() : aa0+nsaaMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
                  size_t nsa = aa1-aa0;
                  for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
                      size_t nsb = bb1-bb0;
                      if (!NextTask()) continue;
                      profiler->start("TransitionDensity ab");
                      TransitionDensity d(w,aa0,aa1, bb0,bb1,0,false,false);
                      profiler->stop("TransitionDensity ab");
                      TransitionDensity e(d);
                      { auto pr1 = profiler->push("MXM ab");
                        MxmDrvNN(&e[0],&d[0], &h.O2(true,false,false).block(symexc)[0], nsa*nsb,nexc,nexc,false);
                      }
                      { auto pr1 = profiler->push("action ab"); e.action(*this);}
                    }
                }
            }
          }
      }
//      xout <<"residual after alpha-beta on process "<<parallel_rank<<" "<<buffer[0]<<std::endl<<str(2)<<std::endl;
  }

  EndTasks();

  gsum(&buffer[0],buffer.size());
}

void Wavefunction::density(memory::vector<double>& den1)
{
 memory::vector<double> den2(0);
 density(den1,den2,true,false,*this);
}
void Wavefunction::density(memory::vector<double>& den1,memory::vector<double>& den2)
{
 density(den1,den2,true,true,*this);
}
void Wavefunction::density(memory::vector<double>& den1,const Wavefunction& bra)
{
 memory::vector<double> den2(0);
 density(den1,den2,true,false,bra);
}
void Wavefunction::density(memory::vector<double>& den1,memory::vector<double>& den2,const Wavefunction& bra)
{
 density(den1,den2,true,true,bra);
}
void Wavefunction::density(memory::vector<double>& den1, memory::vector<double>& den2, bool d1, bool d2, const Wavefunction &bra)
{
  profiler->start("1-electron");
  size_t offset=0, nsa=0, nsb=0;
  unsigned int symdens=symmetry^bra.symmetry;
//  den1.resize(orbitalSpace->total(symdens,1),(double)0);
  den1.assign(0);
  for (unsigned int syma=0; syma<8; syma++) {
    offset += nsa*nsb;
    unsigned int symb = bra.symmetry^syma;
    nsa = bra.alphaStrings[syma].size();
    nsb = bra.betaStrings[symb].size();
    if (!NextTask()) continue;
    profiler->start("1-electron TransitionDensity");
    TransitionDensity d(*this,
                        bra.alphaStrings[syma].begin(),
                        bra.alphaStrings[syma].end(),
                        bra.betaStrings[symb].begin(),
                        bra.betaStrings[symb].end(),
                        1,true, true);
    profiler->stop("1-electron TransitionDensity");
    profiler->start("1-electron MXM");
      MxmDrvTN(&(den1[0]),&(bra.buffer[offset]),&d[0],
          1,nsa*nsb,1,orbitalSpace->total(symdens,1),true);
      profiler->stop("1-electron MXM",2*nsa*nsb*bra.orbitalSpace->total(0,1));
    }
  profiler->stop("1-electron");

}

gci::Operator Wavefunction::density(int rank, bool uhf, const Wavefunction *bra, std::string description, bool parallel_stringset)
{
  if (bra==nullptr) bra=this;
  auto prof = profiler->push("density");

  std::vector<int> symmetries; for (const auto& s : orbitalSpace->orbital_symmetries) symmetries.push_back(s+1); // only a common orbital space is implemented
  dim_t dim; for (const auto s: *orbitalSpace) dim.push_back(s);
  Operator result(dim,
                  symmetries,
             rank,
             uhf,
             symmetry^bra->symmetry,
             description);
  result.zero();

//  xout <<std::endl<<"w in density="<<w.str(2)<<std::endl;
  DivideTasks(99999999,1,1);

  {
  auto p = profiler->push("1-electron");
  size_t offset=0, nsa=0, nsb=0;
  for (unsigned int syma=0; syma<8; syma++) {
    offset += nsa*nsb;
    unsigned int symb = symmetry^syma;
    nsa = alphaStrings[syma].size();
    nsb = betaStrings[symb].size();
    if (!NextTask()) continue;
    {
      profiler->start("1-electron TransitionDensity");
      TransitionDensity d(*this,
                          alphaStrings[syma].begin(),
                          alphaStrings[syma].end(),
                          betaStrings[symb].begin(),
                          betaStrings[symb].end(),
                          1,true, !result.m_uhf);
      profiler->stop("1-electron TransitionDensity");
      profiler->start("1-electron MXM");
//      MxmDrvNN(&buffer[offset],&d[0], &(*h.O1(true).data())[0],
//          nsa*nsb,w.orbitalSpace->total(0,1),1,true);
      MxmDrvTN(&((*result.O1(true).data())[0]),&d[0],&buffer[offset],orbitalSpace->total(0,1),nsa*nsb,1,1);
      profiler->stop("1-electron MXM",2*nsa*nsb*orbitalSpace->total(0,1));
    }
    if (result.m_uhf) {
      profiler->start("1-electron TransitionDensity");
      TransitionDensity d(*this,
                          alphaStrings[syma].begin(),
                          alphaStrings[syma].end(),
                          betaStrings[symb].begin(),
                          betaStrings[symb].end(),
                          1,false, true);
      profiler->stop("1-electron TransitionDensity");
      profiler->start("1-electron MXM");
      MxmDrvTN(&((*result.O1(false).data())[0]),&d[0],&buffer[offset],orbitalSpace->total(0,1),nsa*nsb,1,1);
      profiler->stop("1-electron MXM",2*nsa*nsb*orbitalSpace->total(0,1));
    }
  }

  }
//  xout <<"residual after 1-electron:"<<std::endl<<str(2)<<std::endl;

  { // two-electron contribution, alpha-alpha
    auto p = profiler->push("aa integrals");
    size_t nsbbMax = 64; // temporary static
    for (unsigned int syma=0; syma<8; syma++) {
        profiler->start("StringSet aa");
        StringSet aa(alphaStrings,2,0,syma,parallel_stringset);
        profiler->stop("StringSet aa");
        if (aa.size()==0) continue;
        for (unsigned int symb=0; symb<8; symb++) {
            if (!NextTask()) continue;
            auto praa = profiler->push("aa1 loop");
            unsigned int symexc = syma^symb^symmetry;
            //size_t nexc = h.pairSpace.find(-1)->second[symexc];
            size_t nexc = result.O2(true,true).block_size(symexc);
            size_t nsb = betaStrings[symb].size(); if (nsb==0) continue;
            for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsbbMax > aa.end() ? aa.end() : aa0+nsbbMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
                size_t nsa = aa1-aa0;
                profiler->start("TransitionDensity aa");
                TransitionDensity d(*this,aa0,aa1,betaStrings[symb].begin(),betaStrings[symb].end(),-1,false,false);
                TransitionDensity e(*bra,aa0,aa1,betaStrings[symb].begin(),betaStrings[symb].end(),-1,false,false);
                profiler->stop("TransitionDensity aa");
                { auto pr1 = profiler->push("MXM aa");
                  MxmDrvTN(&result.O2(true,true).block(symexc)[0], &d[0],&e[0], nexc, nsa*nsb, nsa*nsb, nexc, false); // not yet right for non-symmetric tdm
                }
              }
          }
      }
  }

  if (result.m_uhf) { // two-electron contribution, beta-beta
      auto p = profiler->push("bb integrals");
      size_t nsbbMax = 64; // temporary static
      for (unsigned int symb=0; symb<8; symb++) {
          StringSet bb(betaStrings,2,0,symb,parallel_stringset);
          if (bb.size()==0) continue;
          for (unsigned int syma=0; syma<8; syma++) {
              if (!NextTask()) continue;
              unsigned int symexc = symb^syma^symmetry;
              size_t nexc = result.O2(false,false).block_size(symexc);
              size_t nsa = alphaStrings[syma].size(); if (nsa==0) continue;
              for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
                  size_t nsb = bb1-bb0;
                  profiler->start("TransitionDensity bb");
                  TransitionDensity d(*this,alphaStrings[syma].begin(),alphaStrings[syma].end(),bb0,bb1,-1,false,false);
                  TransitionDensity e(*bra,alphaStrings[syma].begin(),alphaStrings[syma].end(),bb0,bb1,-1,false,false);
                  profiler->stop("TransitionDensity bb");
                  { auto pr1 = profiler->push("MXM bb");
                    MxmDrvTN(&result.O2(false,false).block(symexc)[0], &d[0],&e[0], nexc, nsa*nsb, nsa*nsb, nexc, false); // not yet right for non-symmetric tdm
                  }
                }
            }
        }
    }

  { // two-electron contribution, alpha-beta
    auto p = profiler->push("ab integrals");
    size_t nsaaMax = 640; // temporary static
    size_t nsbbMax = 640; // temporary static
    for (unsigned int symb=0; symb<8; symb++) {
        StringSet bb(betaStrings,1,0,symb,parallel_stringset);
        if (bb.size()==0) continue;
        for (unsigned int syma=0; syma<8; syma++) {
            StringSet aa(alphaStrings,1,0,syma,parallel_stringset);
            if (aa.size()==0) continue;
            unsigned int symexc = symb^syma^symmetry;
            size_t nexc = result.O2(true,false,false).block_size(symexc);
            {
              auto pro = profiler->push("StringSet iterator loops");
              for (StringSet::iterator aa1, aa0=aa.begin(); aa1=aa0+nsaaMax > aa.end() ? aa.end() : aa0+nsaaMax, aa0 <aa.end(); aa0=aa1) { // loop over alpha batches
                  size_t nsa = aa1-aa0;
                  for (StringSet::iterator bb1, bb0=bb.begin(); bb1=bb0+nsbbMax > bb.end() ? bb.end() : bb0+nsbbMax, bb0 <bb.end(); bb0=bb1) { // loop over beta batches
                      size_t nsb = bb1-bb0;
                      if (!NextTask()) continue;
                      profiler->start("TransitionDensity ab");
                      TransitionDensity d(*this,aa0,aa1, bb0,bb1,0,false,false);
                      TransitionDensity e(*bra,aa0,aa1, bb0,bb1,0,false,false);
                      profiler->stop("TransitionDensity ab");
                      { auto pr1 = profiler->push("MXM ab");
                        MxmDrvTN(&result.O2(true,false,false).block(symexc)[0], &d[0],&e[0], nexc, nsa*nsb, nsa*nsb, nexc, false); // not yet right for non-symmetric tdm
                      }
                    }
                }
            }
          }
      }
  }

  EndTasks();

  //FIXME implement gsum
//  gsum(&buffer[0],buffer.size());
  std::cout << "Density:\n"<<result<<std::endl;
  return result;
}

SMat Wavefunction::naturalOrbitals()
{
//    xout <<"naturalOrbitals"<<std::endl;
  SMat natorb(std::vector<std::vector<size_t> >{*orbitalSpace,*orbitalSpace},parityNone,0,"Natural orbitals");
    SMat dens1(std::vector<std::vector<size_t> >{*orbitalSpace,*orbitalSpace},parityEven,0,"Density");
    SMat ee({*orbitalSpace},parityNone,-1,"Occupation numbers");
    density(*dens1.data());
    // at this point, the off-diagonal elements of the density are twice what they should be, so sort that out.
    dens1 *=0.5; for (int k=0; k<dens1.max_symmetry(); k++) for (size_t l=0;l<dens1.dimension(k);l++) dens1.block(k)[(l+2)*(l+1)/2-1]*=2;

    xout << dens1;
    dens1.ev(ee,&natorb);
    xout <<ee<<natorb;
    return natorb;
}

void Wavefunction::putw(File& f, int index)
{
  auto p = profiler->push("Wavefunction::putw");
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  size_t offset = chunk*parallel_rank;
  if (chunk+offset > buffer.size()) chunk = buffer.size()-offset;
  if (offset < buffer.size())
    f.write(&buffer[offset],chunk,index*chunk);
  p+=chunk;
}


void Wavefunction::getw(File& f, int index)
{
  auto p = profiler->push("Wavefunction::getw");
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  size_t offset = chunk*parallel_rank;
  if (chunk+offset > buffer.size()) chunk = buffer.size()-offset;
  if (offset < buffer.size())
    f.read(&buffer[offset],chunk,index*chunk);
  p+=chunk;
}


void Wavefunction::getAll(File& f, int index)
{
  auto p = profiler->push("Wavefunction::getAll");
  getw(f,index);
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  gather_chunks(&buffer[0],buffer.size(),chunk);
  p+=buffer.size();
}

void Wavefunction::gather()
{
  size_t chunk = (buffer.size()-1)/parallel_size+1;
  gather_chunks(&buffer[0],buffer.size(),chunk);

}

#include <cmath>
std::vector<std::size_t> Wavefunction::histogram(const std::vector<double> edges, bool parallel, std::size_t start, std::size_t stop)
{
  if (stop > edges.size()) stop = edges.size();
  std::vector<std::size_t> cumulative(edges.size(),0);
  if (parallel) DivideTasks(edges.size());
  for (size_t i=start; i<stop;i++)
    if (parallel && NextTask())
      for (size_t j=0; j<buffer.size(); j++)
        if (std::fabs(buffer[j]) > edges[i]) cumulative[i]++;
  if (parallel) {
  EndTasks();
  std::vector<double> dcumulative(cumulative.size());
  for (size_t i=0; i<cumulative.size(); i++) dcumulative[i]=(double)cumulative[i];
  gsum(&dcumulative[0],cumulative.size());
  for (size_t i=0; i<cumulative.size(); i++) cumulative[i]=(std::size_t)dcumulative[i];
  }
  return cumulative;
}

double Wavefunction::norm(const double k)
{
  double result = (double)0;
  size_t imin = 0;
  size_t imax = buffer.size();
  if (distributed) {
    size_t chunk = (buffer.size()-1)/parallel_size+1;
    imin=parallel_rank*chunk;
    imax=(parallel_rank+1)*chunk;
  }
  for (size_t i=imin; i<imax; i++)
    if (buffer[i])
      result += pow(fabs(buffer[i]),k);
  return result;
}

Wavefunction& Wavefunction::addAbsPower(const Wavefunction& c, const double k, const double factor)
{

  if (! compatible(c)) throw std::domain_error("attempt to add incompatible Wavefunction objects");
//  xout <<"addAbsPower initial=";
//  for (size_t i=0; i<buffer.size(); i++)
//    xout <<" "<<buffer[i];
//  xout <<std::endl;
  size_t imin = 0;
  size_t imax = buffer.size();
  if (distributed) {
    size_t chunk = (buffer.size()-1)/parallel_size+1;
    imin=parallel_rank*chunk;
    imax=(parallel_rank+1)*chunk;
  }
  if (imax > buffer.size()) imax = buffer.size();
  for (size_t i=imin; i<imax; i++)
    if (k == -1)
      buffer[i] += c.buffer[i] <0 ? -factor : factor;
    else if (c.buffer[i])
      buffer[i] += factor*pow(fabs(c.buffer[i]),k)*c.buffer[i];
//  xout <<"addAbsPower final=";
//  for (size_t i=0; i<buffer.size(); i++)
//    xout <<" "<<buffer[i];
//  xout <<std::endl;
  return *this;
}

std::vector<double>::iterator Wavefunction::begin()
{
  if (distributed)
    return buffer.begin() + parallel_rank*((buffer.size()-1)/parallel_size+1);
  else
    return buffer.begin();
}

std::vector<double>::iterator Wavefunction::end()
{
  if (distributed)
    return buffer.begin() + std::min((size_t)buffer.size(),(size_t)((parallel_rank+1)*((buffer.size()-1)/parallel_size+1)));
  else
    return buffer.end();
}
