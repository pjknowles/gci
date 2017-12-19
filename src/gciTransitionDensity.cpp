#include "gciTransitionDensity.h"
#include <iostream>
#include <sstream>
#include <iterator>

TransitionDensity::TransitionDensity(const Wavefunction &w,
                                     const StringSet::const_iterator &alphaStringsBegin,
                                     const StringSet::const_iterator &alphaStringsEnd,
                                     const StringSet::const_iterator &betaStringsBegin,
                                     const StringSet::const_iterator &betaStringsEnd,
                                     const parity_t parity, const bool doAlpha, const bool doBeta)
  : m_alphaStringsBegin(alphaStringsBegin)
  , m_alphaStringsEnd(alphaStringsEnd)
  , m_betaStringsBegin(betaStringsBegin)
  , m_betaStringsEnd(betaStringsEnd)
  , m_parity(parity)
  , m_nsa(std::distance(m_alphaStringsBegin,m_alphaStringsEnd))
  , m_nsb(std::distance(m_betaStringsBegin,m_betaStringsEnd))
  , m_syma(m_alphaStringsBegin->computed_symmetry())
  , m_symb(m_betaStringsBegin->computed_symmetry())
  , m_symexc(w.symmetry ^ m_syma ^ m_symb)
  , m_deltaAlpha(w.alphaStrings[0].proto.nelec - m_alphaStringsBegin->nelec)
  , m_deltaBeta(w.betaStrings[0].proto.nelec - m_betaStringsBegin->nelec)
  , m_excitations((m_deltaAlpha+m_deltaBeta)%2 ? w.orbitalSpace->operator[](m_symexc) : w.orbitalSpace->total(m_symexc,m_parity))
{
  resize(m_nsa*m_nsb*m_excitations,(double)0);
  if (m_nsa*m_nsb*m_excitations == 0) return;
  auto prof = profiler->push("TransitionDensity");

  if (m_deltaAlpha==0 && m_deltaBeta==0) { // number of electrons preserved, so one-electron excitation
      m_hasAlpha=doAlpha;
      m_hasBeta=doBeta;
      if (doAlpha) {
          auto prof2=profiler->push("TransitionDensity_alpha");
          // alpha excitations
          unsigned int wsymb = m_symb;
          unsigned int wsyma = w.symmetry^wsymb;
          size_t wnsb = w.betaStrings[wsymb].size();
          // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
          size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
          size_t offa=0;
          for (StringSet::const_iterator s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
	      ExcitationSet ee(*s,w.alphaStrings[wsyma],1,1,m_parity);
              prof2 += ee.size()*m_nsb*2;
              for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
                  if (e->phase < 0)
                    for (size_t ib=0; ib<m_nsb; ib++)
                      (*this)[offa+m_nsa*m_nsb*e->orbitalAddress+ib] -= w.buffer[woffset+e->stringIndex*wnsb+ib];
                  else
                    for (size_t ib=0; ib<m_nsb; ib++)
                      (*this)[offa+m_nsa*m_nsb*e->orbitalAddress+ib ]+= w.buffer[woffset+e->stringIndex*wnsb+ib];
                }
              offa += m_nsb;
            }
        }

      if (doBeta) {
          auto prof2=profiler->push("TransitionDensity_beta");
          // beta excitations
          unsigned int wsyma = m_syma;
          unsigned int wsymb = w.symmetry^wsyma;
          size_t wnsb = w.betaStrings[wsymb].size();
          // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
          size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
          size_t offb = 0;
          for (StringSet::const_iterator s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
	      ExcitationSet ee(*s,w.betaStrings[wsymb],1,1,m_parity);
              prof2 += ee.size()*m_nsa*2;
              auto v = &((*this)[offb]);
              for (const auto& e : ee) {
                  if (e.phase < 0)
                    for (size_t ia=0; ia<m_nsa; ia++) {
                        *(v+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb) -= w.buffer[woffset+e.stringIndex+wnsb*ia];
                        //                    (*this)[offb+nsa*nsb*e.orbitalAddress+ia*nsb] -= w.buffer[woffset+e.stringIndex+wnsb*ia];
                      }
                  else
                    for (size_t ia=0; ia<m_nsa; ia++) {
                        *(v+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb) += w.buffer[woffset+e.stringIndex+wnsb*ia];
                        //                    (*this)[offb+nsa*nsb*e.orbitalAddress+ia*nsb] += w.buffer[woffset+e.stringIndex+wnsb*ia];
                      }
                }
              offb ++;
            }
        }

    } else if (m_deltaBeta==0 && m_deltaAlpha>0) { // wavefunction has deltaAlpha more alpha electrons than interacting states
      auto prof2=profiler->push("TransitionDensity_alpha_alpha");
      m_hasAlpha=true; m_hasBeta=false;
      unsigned int wsymb = m_symb;
      unsigned int wsyma = w.symmetry^wsymb;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
      //    xout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

      size_t offa=0;
      for (StringSet::const_iterator s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
          ExcitationSet ee(*s,w.alphaStrings[wsyma],0,m_deltaAlpha,m_parity);
          prof2 += ee.size()*m_nsb*2;
          for (const auto& e : ee) {
              if (e.phase < 0)
                for (size_t ib=0; ib<m_nsb; ib++)
                  (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib]-=
                    w.buffer[woffset+e.stringIndex*wnsb+ib];
              else
                for (size_t ib=0; ib<m_nsb; ib++)
                  (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib]+=
                    w.buffer[woffset+e.stringIndex*wnsb+ib];
            }
          offa += m_nsb;
        }
    } else if (m_deltaAlpha==0 && m_deltaBeta>0) { // wavefunction has deltaBeta more beta electrons than interacting states
      auto prof2=profiler->push("TransitionDensity_beta_beta");
      m_hasAlpha=false; m_hasBeta=true;
      unsigned int wsyma = m_syma;
      unsigned int wsymb = w.symmetry^wsyma;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      for (StringSet::const_iterator s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
          ExcitationSet ee(*s,w.betaStrings[wsymb],0,m_deltaBeta,m_parity);
          prof2 += ee.size()*m_nsa*2;
          for (const auto& e : ee) {
              if (e.phase < 0)
                for (size_t ia=0; ia<m_nsa; ia++)
                  (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb]-=
                    w.buffer[woffset+e.stringIndex+wnsb*ia];
              else
                for (size_t ia=0; ia<m_nsa; ia++)
                  (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb]+=
                    w.buffer[woffset+e.stringIndex+wnsb*ia];
            }
          offb ++;
        }
    } else if (m_deltaAlpha==1 && m_deltaBeta==1) { // wavefunction has 1 more alpha and beta electrons than interacting states
      auto prof2=profiler->push("TransitionDensity_alpha_beta");
      if (m_parity!=parityNone) throw std::logic_error("wrong parity in alpha-beta");
      m_hasAlpha=true; m_hasBeta=true;
      for (unsigned int wsyma=0; wsyma<8; wsyma++) {
          unsigned int wsymb = w.symmetry^wsyma;
          unsigned int symexca = wsyma ^ m_syma;
          const size_t intoff = w.orbitalSpace->offset(m_symexc,symexca,m_parity);
          const size_t wnsb = w.betaStrings[wsymb].size();
          if (wnsb==0) continue;
          const size_t woffset = w.blockOffset(wsyma);
          size_t offb = m_nsa*m_nsb*intoff;
          std::vector<ExcitationSet> eebs;
          for (StringSet::const_iterator sb = m_betaStringsBegin; sb != m_betaStringsEnd; sb++) {
            eebs.emplace_back(*sb,w.betaStrings[wsymb],0,1,m_parity);
            for (auto& eb : eebs.back()) // optimisation for non-standard use in following loop
             const_cast<Excitation&>(eb).orbitalAddress *= m_nsa*m_nsb*(*w.orbitalSpace)[symexca];
          }
          for (StringSet::const_iterator sa = m_alphaStringsBegin; sa != m_alphaStringsEnd; sa++) {
              ExcitationSet eea(*sa,w.alphaStrings[wsyma],0,1,m_parity);
              for (auto& ea : eea) { // optimisation for non-standard use in following loop
               const_cast<Excitation&>(ea).orbitalAddress*=m_nsa*m_nsb;
               const_cast<Excitation&>(ea).stringIndex*=wnsb;
              }
              for (const auto& eeb : eebs) {
                  prof2 += eea.size()*eeb.size()*2;
                      for (const auto& ea : eea) {
                  for (const auto& eb : eeb) {
                          (*this)[offb + eb.orbitalAddress + ea.orbitalAddress] += ea.phase * eb.phase * w.buffer[woffset + eb.stringIndex + ea.stringIndex];
                        }
                    }
                  offb ++;
                }
            }
        }
    } else {
      xout <<"deltaAlpha="<<m_deltaAlpha<<", deltaBeta="<<m_deltaBeta<<std::endl;
      throw std::logic_error("unimplemented case");
    }
  //  size_t populated=0; for (const_iterator x=begin(); x!=end(); x++) if (*x !=(double)0) ++populated; xout <<"TransitionDensity population="<<((double)populated)/((double)size()+1)<<std::endl;

}

void TransitionDensity::action(Wavefunction &w) const
{
  if (m_nsa*m_nsb*m_excitations == 0) return;
  auto prof = profiler->push("TransitionDensity::action");

  if (m_deltaAlpha==0 && m_deltaBeta==0) { // number of electrons preserved, so one-electron excitation

      if (true) {
          // alpha excitations
          unsigned int wsymb = m_symb;
          unsigned int wsyma = w.symmetry^wsymb;
          size_t wnsb = w.betaStrings[wsymb].size();
          // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
          size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
          size_t offa=0;
          for (StringSet::const_iterator s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
              ExcitationSet ee(*s,w.alphaStrings[wsyma],1,1,m_parity);
              prof += ee.size()*m_nsb*2;
              for (const auto& e : ee) {
                  if (e.phase < 0)
                    for (size_t ib=0; ib<m_nsb; ib++)
                      w.buffer[woffset+e.stringIndex*wnsb+ib] -=
                          (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib];
                  else
                    for (size_t ib=0; ib<m_nsb; ib++)
                      w.buffer[woffset+e.stringIndex*wnsb+ib] +=
                          (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib];
                }
              offa += m_nsb;
            }
        }

      if (true) {
          // beta excitations
          unsigned int wsyma = m_syma;
          unsigned int wsymb = w.symmetry^wsyma;
          size_t wnsb = w.betaStrings[wsymb].size();
          // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
          size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
          size_t offb = 0;
          for (StringSet::const_iterator s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
              ExcitationSet ee(*s,w.betaStrings[wsymb],1,1,m_parity);
              prof += ee.size()*m_nsa*2;
              for (const auto& e : ee) {
                  if (e.phase < 0)
                    for (size_t ia=0; ia<m_nsa; ia++)
                      w.buffer[woffset+e.stringIndex+wnsb*ia]-=
                          (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb];
                  else
                    for (size_t ia=0; ia<m_nsa; ia++)
                      w.buffer[woffset+e.stringIndex+wnsb*ia]+=
                          (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb];
                }
              offb ++;
            }
        }

    } else if (m_deltaBeta==0 && m_deltaAlpha > 0) { // wavefunction has deltaAlpha more alpha electrons than interacting states
      unsigned int wsymb = m_symb;
      unsigned int wsyma = w.symmetry^wsymb;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
      size_t offa=0;
      for (StringSet::const_iterator s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
           ExcitationSet ee(*s,w.alphaStrings[wsyma],0,m_deltaAlpha,m_parity);
          //      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
          for (const auto& e : ee) {
              if (e.phase < 0)
                for (size_t ib=0; ib<m_nsb; ib++)
                  w.buffer[woffset+e.stringIndex*wnsb+ib]-=
                      (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib];
              else
                for (size_t ib=0; ib<m_nsb; ib++)
                  w.buffer[woffset+e.stringIndex*wnsb+ib]+=
                      (*this)[offa+m_nsa*m_nsb*e.orbitalAddress+ib];
            }
          offa += m_nsb;
        }
    } else if (m_deltaAlpha==0 && m_deltaBeta>0) { // wavefunction has deltaBeta more beta electrons than interacting states
      unsigned int wsyma = m_syma;
      unsigned int wsymb = w.symmetry^wsyma;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      for (StringSet::const_iterator s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
          ExcitationSet ee(*s,w.betaStrings[wsymb],0,m_deltaBeta,m_parity);
          for (const auto& e : ee) {
              if (e.phase < 0)
                for (size_t ia=0; ia<m_nsa; ia++)
                  w.buffer[woffset+e.stringIndex+wnsb*ia]-=
                      (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb];
              else
                for (size_t ia=0; ia<m_nsa; ia++)
                  w.buffer[woffset+e.stringIndex+wnsb*ia]+=
                      (*this)[offb+m_nsa*m_nsb*e.orbitalAddress+ia*m_nsb];
            }
          offb ++;
        }
    } else if (m_deltaAlpha==1 && m_deltaBeta==1) { // wavefunction has 1 more alpha and beta electrons than interacting states
      if (m_parity) throw std::logic_error("wrong parity in alpha-beta");
      for (unsigned int wsyma=0; wsyma<8; wsyma++) {
          unsigned int wsymb = w.symmetry^wsyma;
          unsigned int symexca = wsyma ^ m_syma;
          size_t intoff = w.orbitalSpace->offset(m_symexc,symexca,m_parity);
          size_t wnsb = w.betaStrings[wsymb].size();
          size_t woffset = w.blockOffset(wsyma);
          size_t offb = 0;
          std::vector<ExcitationSet> eebs;
          for (StringSet::const_iterator sb = m_betaStringsBegin; sb != m_betaStringsEnd; sb++)
            eebs.emplace_back(*sb,w.betaStrings[wsymb],0,1);
          for (StringSet::const_iterator sa = m_alphaStringsBegin; sa != m_alphaStringsEnd; sa++) {
              ExcitationSet eea(*sa,w.alphaStrings[wsyma],0,1);
              for (const auto& eeb : eebs) {
                  prof += eea.size()*eeb.size()*2;
                  for (const auto& eb : eeb) {
                      auto intoffb = offb+m_nsa*m_nsb*(intoff+eb.orbitalAddress*(*w.orbitalSpace)[symexca]);
                      for (const auto& ea : eea) {
                          w.buffer[woffset + eb.stringIndex + wnsb * ea.stringIndex]
                              += ea.phase * eb.phase *
                              (*this)[intoffb + m_nsa*m_nsb*ea.orbitalAddress];
                        }
                    }
                  offb ++;
                }
           }
        }
    } else {
      xout <<"deltaAlpha="<<m_deltaAlpha<<", deltaBeta="<<m_deltaBeta<<std::endl;
      throw std::logic_error("unimplemented case");
    }
}

#include "gciMolpro.h"

SymmetryMatrix::Operator TransitionDensity::density(const Wavefunction &w) const
{
  dim_t dimension; for (auto i=0; i<8; i++) dimension[i]=w[i];
  std::vector<int> symmetries; for (const auto& s : w.orbitalSpace->orbital_symmetries) symmetries.push_back(s+1);
  gci::Operator result(dimension, symmetries, 1, !(m_hasAlpha&&m_hasBeta), m_symexc );

  unsigned int syma = m_alphaStringsBegin->computed_symmetry();
  unsigned int symb = m_betaStringsBegin->computed_symmetry();
  if ((syma^symb) != w.symmetry) throw std::runtime_error("Symmetry mismatch, TransitionDensity::density");
  if (m_nsa*m_nsb != (w.blockOffset(syma+1)-w.blockOffset(syma))) throw std::range_error("Wrong dimensions, TransitionDensity::density"); // present implementation only for complete space in TransitionDensity
  MxmDrvTN(&((*result.O1(true).data())[0]),&this->at(0),&(w.buffer[w.blockOffset(syma)]),m_excitations,m_nsa*m_nsb,1,1);
  return result;
}

std::string TransitionDensity::str(int verbosity, unsigned int columns) const
{
  if (verbosity < 0) return std::string("");
  std::ostringstream s;
  //  xout <<"TransitionDensity::str size()="<<size()<<std::endl;
  if (size()!= m_nsa*m_nsb*m_excitations) {xout << "wrong size in TransitionDensity::str "<<size()<<" "<<m_nsa*m_nsb*m_excitations<<std::endl; throw std::range_error("help");}
  if (size())
    for (size_t ij=0; ij < m_excitations; ij++) {
        s <<std::endl;
        for (size_t ab=0; ab < m_nsa*m_nsb; ab++)
          s << (*this)[ab+ij*m_nsa*m_nsb] <<" ";
      }
  return s.str();
}
