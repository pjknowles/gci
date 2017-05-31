#include "gciTransitionDensity.h"
#include <iostream>
#include <sstream>
#include <iterator>

TransitionDensity::TransitionDensity(const Wavefunction &w,
                                     const StringSet::const_iterator &alphaStringsBegin_,
                                     const StringSet::const_iterator &alphaStringsEnd_,
                                     const StringSet::const_iterator &betaStringsBegin_,
                                     const StringSet::const_iterator &betaStringsEnd_,
                                     const int parity_, const bool doAlpha, const bool doBeta)
{
  profiler->start("TransitionDensity preamble");
  this->parity = parity_;
  this->alphaStringsBegin = alphaStringsBegin_;
  this->alphaStringsEnd = alphaStringsEnd_;
  this->betaStringsBegin = betaStringsBegin_;
  this->betaStringsEnd = betaStringsEnd_;
  // first parse the type of transition
  nsa = std::distance(alphaStringsBegin,alphaStringsEnd);
  nsb = std::distance(betaStringsBegin,betaStringsEnd);
  if (nsa==0 || nsb==0) {profiler->stop();return;}
  unsigned int syma = alphaStringsBegin->computed_symmetry();
  unsigned int symb =  betaStringsBegin->computed_symmetry();
  //  xout << "syma="<<syma<<", nsa="<<nsa<<std::endl;
  //  xout << "symb="<<symb<<", nsb="<<nsb<<std::endl;
  symexc = w.symmetry ^ syma ^ symb;
  excitations = w.orbitalSpace->total(symexc,parity);
  int deltaAlpha = w.alphaStrings[0].proto.nelec - alphaStringsBegin->nelec;
  int deltaBeta = w.betaStrings[0].proto.nelec - betaStringsBegin->nelec;

  //  xout <<"TransitionDensity "<<symexc<<" "<<nsa*nsb*excitations<<std::endl;
  resize(nsa*nsb*excitations,(double)0);
  profiler->stop("TransitionDensity preamble",nsa*nsb*excitations);
  if (nsa*nsb*excitations == 0) return;

  if (deltaAlpha==0 && deltaBeta==0) { // number of electrons preserved, so one-electron excitation

    if (doAlpha) {
      profiler->start("1-electron alpha excitations");
      // alpha excitations
      unsigned int wsymb = symb;
      unsigned int wsyma = w.symmetry^wsymb;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + betaStringsBegin->index(w.betaStrings[wsymb]);
      //    xout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

      size_t offa=0;
      for (StringSet::const_iterator s = alphaStringsBegin; s != alphaStringsEnd; s++) {
        //          xout << "alpha string "<<*s<<std::endl;
        ExcitationSet ee(*s,w.alphaStrings[wsyma],1,1);
        //          xout << "alpha excitations " << ee.str() <<std::endl;
        for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
          //        xout << "alpha excitation " << e->orbitalAddress <<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
          if (e->phase < 0)
            for (size_t ib=0; ib<nsb; ib++)
              (*this)[offa+nsa*nsb*e->orbitalAddress+ib]-=
                w.buffer[woffset+e->stringIndex*wnsb+ib];
          else
            for (size_t ib=0; ib<nsb; ib++) {
              (*this)[offa+nsa*nsb*e->orbitalAddress+ib]+=
                  w.buffer[woffset+e->stringIndex*wnsb+ib];
              //              xout <<"using w.buffer="<<w.buffer[woffset+e->stringIndex*wnsb+ib]<<" "<<woffset<<" "<<e->stringIndex*wnsb<<std::endl;
            }
        }
        offa += nsb;
      }
      //    xout <<"after alpha"<<std::endl;
      //  for (size_t ij=0; ij < excitations; ij++) {
      //    for (size_t ab=0; ab < nsa*nsb; ab++)
      //      xout << (*this)[ab+ij*nsa*nsb] <<" ";
      //    xout <<std::endl;
      //  }
      profiler->stop("1-electron alpha excitations");
    }

    if (doBeta) {
      profiler->start("1-electron beta excitations");
      // beta excitations
      unsigned int wsyma = syma;
      unsigned int wsymb = w.symmetry^wsyma;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + wnsb * alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      //    xout << "beta wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<", woffset="<<woffset<<std::endl;
      for (StringSet::const_iterator s = betaStringsBegin; s != betaStringsEnd; s++) {
        ExcitationSet ee(*s,w.betaStrings[wsymb],1,1);
        for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
          //        xout << "beta excitation " << e->orbitalAddress <<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
          if (e->phase < 0)
            for (size_t ia=0; ia<nsa; ia++)
              (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb]-=
                w.buffer[woffset+e->stringIndex+wnsb*ia];
          else
            for (size_t ia=0; ia<nsa; ia++)
              (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb]+=
                w.buffer[woffset+e->stringIndex+wnsb*ia];
        }
        offb ++;
      }
      profiler->stop("1-electron beta excitations");
    }

  } else if (deltaAlpha==2) { // wavefunction has 2 more alpha electrons than interacting states
    unsigned int wsymb = symb;
    unsigned int wsyma = w.symmetry^wsymb;
    size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    size_t woffset = w.blockOffset(wsyma) + betaStringsBegin->index(w.betaStrings[wsymb]);
    //    xout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

    size_t offa=0;
    for (StringSet::const_iterator s = alphaStringsBegin; s != alphaStringsEnd; s++) {
      //          xout << "alpha string "<<*s<<std::endl;
      ExcitationSet ee(*s,w.alphaStrings[wsyma],0,2);
      //          xout << "alpha excitations " << ee.str() <<std::endl;
      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
        //        xout << "alpha excitation " << e->orbitalAddress <<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
        if (e->phase < 0)
          for (size_t ib=0; ib<nsb; ib++)
            (*this)[offa+nsa*nsb*e->orbitalAddress+ib]-=
              w.buffer[woffset+e->stringIndex*wnsb+ib];
        else
          for (size_t ib=0; ib<nsb; ib++) {
            (*this)[offa+nsa*nsb*e->orbitalAddress+ib]+=
                w.buffer[woffset+e->stringIndex*wnsb+ib];
            //              xout <<"using w.buffer="<<w.buffer[woffset+e->stringIndex*wnsb+ib]<<" "<<woffset<<" "<<e->stringIndex*wnsb<<std::endl;
          }
      }
      offa += nsb;
    }
  } else if (deltaBeta==2) { // wavefunction has 2 more beta electrons than interacting states
    unsigned int wsyma = syma;
    unsigned int wsymb = w.symmetry^wsyma;
    size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    size_t woffset = w.blockOffset(wsyma) + wnsb * alphaStringsBegin->index(w.alphaStrings[wsyma]);
    size_t offb = 0;
    //    xout << "beta wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<", woffset="<<woffset<<std::endl;
    for (StringSet::const_iterator s = betaStringsBegin; s != betaStringsEnd; s++) {
      ExcitationSet ee(*s,w.betaStrings[wsymb],0,2);
      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
        //        xout << "beta excitation " << e->orbitalAddress <<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
        if (e->phase < 0)
          for (size_t ia=0; ia<nsa; ia++)
            (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb]-=
              w.buffer[woffset+e->stringIndex+wnsb*ia];
        else
          for (size_t ia=0; ia<nsa; ia++)
            (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb]+=
              w.buffer[woffset+e->stringIndex+wnsb*ia];
      }
      offb ++;
    }
  } else if (deltaAlpha==1 && deltaBeta==1) { // wavefunction has 1 more alpha and beta electrons than interacting states
    if (parity) throw std::logic_error("wrong parity in alpha-beta");
    for (unsigned int wsyma=0; wsyma<8; wsyma++) {
      unsigned int wsymb = w.symmetry^wsyma;
      unsigned int symexca = wsyma ^ syma;
      size_t intoff = w.orbitalSpace->offset(symexc,symexca,parity);
      size_t wnsb = w.betaStrings[wsymb].size();
      size_t woffset = w.blockOffset(wsyma);
      size_t offb = 0;
      std::vector<ExcitationSet> eebs;
      for (StringSet::const_iterator sb = betaStringsBegin; sb != betaStringsEnd; sb++)
        eebs.push_back(ExcitationSet(*sb,w.betaStrings[wsymb],0,1));
      for (StringSet::const_iterator sa = alphaStringsBegin; sa != alphaStringsEnd; sa++) {
        ExcitationSet eea(*sa,w.alphaStrings[wsyma],0,1);
        for (std::vector<ExcitationSet>::const_iterator eebp=eebs.begin(); eebp!=eebs.end(); eebp++) {
          for (ExcitationSet::const_iterator ea=eea.begin(); ea!=eea.end(); ea++) {
            for (ExcitationSet::const_iterator eb=eebp->begin(); eb!=eebp->end(); eb++) {
              (*this)[offb + nsa*nsb*(intoff +
                                      ea->orbitalAddress +
                                      eb->orbitalAddress*(*w.orbitalSpace)[symexca]) ]
                  += ea->phase * eb->phase *
                  w.buffer[woffset + eb->stringIndex + wnsb * ea->stringIndex];
            }
          }
          offb ++;
        }
      }
    }
  } else {
    xout <<"deltaAlpha="<<deltaAlpha<<", deltaBeta="<<deltaBeta<<std::endl;
    throw std::logic_error("unimplemented case");
  }
//  size_t populated=0; for (const_iterator x=begin(); x!=end(); x++) if (*x !=(double)0) ++populated; xout <<"TransitionDensity population="<<((double)populated)/((double)size()+1)<<std::endl;

}

void TransitionDensity::action(Wavefunction &w)
{
  // first parse the type of transition
  nsa = std::distance(alphaStringsBegin,alphaStringsEnd);
  nsb = std::distance(betaStringsBegin,betaStringsEnd);
  if (nsa==0 || nsb==0) return;
  unsigned int syma = alphaStringsBegin->computed_symmetry();
  unsigned int symb =  betaStringsBegin->computed_symmetry();
  unsigned int symexcw = w.symmetry ^ syma ^ symb;
  excitations = w.orbitalSpace->total(symexcw,parity);
  int deltaAlpha = w.alphaStrings[0].proto.nelec - alphaStringsBegin->nelec;
  int deltaBeta = w.betaStrings[0].proto.nelec - betaStringsBegin->nelec;

  if (nsa*nsb*excitations == 0) return;

  if (deltaAlpha==0 && deltaBeta==0) { // number of electrons preserved, so one-electron excitation

    if (true) {
      // alpha excitations
      unsigned int wsymb = symb;
      unsigned int wsyma = w.symmetry^wsymb;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + betaStringsBegin->index(w.betaStrings[wsymb]);
      //    xout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

      size_t offa=0;
      for (StringSet::const_iterator s = alphaStringsBegin; s != alphaStringsEnd; s++) {
        //          xout << "alpha string "<<*s<<std::endl;
        ExcitationSet ee(*s,w.alphaStrings[wsyma],1,1);
        //          xout << "alpha excitations " << ee.str() <<std::endl;
        for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
          //        xout << "alpha excitation " << e->orbitalAddress <<" "<<e->phase <<" "<<e->stringIndex<<std::endl;
          if (e->phase < 0)
            for (size_t ib=0; ib<nsb; ib++)
              w.buffer[woffset+e->stringIndex*wnsb+ib] -=
                  (*this)[offa+nsa*nsb*e->orbitalAddress+ib];
          else
            for (size_t ib=0; ib<nsb; ib++)
              w.buffer[woffset+e->stringIndex*wnsb+ib] +=
                  (*this)[offa+nsa*nsb*e->orbitalAddress+ib];
        }
        offa += nsb;
      }
    }

    if (true) {
      // beta excitations
      unsigned int wsyma = syma;
      unsigned int wsymb = w.symmetry^wsyma;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + wnsb * alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      //    xout << "beta wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<", woffset="<<woffset<<std::endl;
      for (StringSet::const_iterator s = betaStringsBegin; s != betaStringsEnd; s++) {
        ExcitationSet ee(*s,w.betaStrings[wsymb],1,1);
        for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
          if (e->phase < 0)
            for (size_t ia=0; ia<nsa; ia++)
              w.buffer[woffset+e->stringIndex+wnsb*ia]-=
                  (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb];
          else
            for (size_t ia=0; ia<nsa; ia++)
              w.buffer[woffset+e->stringIndex+wnsb*ia]+=
                  (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb];
        }
        offb ++;
      }
    }

  } else if (deltaAlpha==2) { // wavefunction has 2 more alpha electrons than interacting states
    unsigned int wsymb = symb;
    unsigned int wsyma = w.symmetry^wsymb;
    size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    size_t woffset = w.blockOffset(wsyma) + betaStringsBegin->index(w.betaStrings[wsymb]);
    //    xout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

    size_t offa=0;
    for (StringSet::const_iterator s = alphaStringsBegin; s != alphaStringsEnd; s++) {
      ExcitationSet ee(*s,w.alphaStrings[wsyma],0,2);
      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
        if (e->phase < 0)
          for (size_t ib=0; ib<nsb; ib++)
            w.buffer[woffset+e->stringIndex*wnsb+ib]-=
                (*this)[offa+nsa*nsb*e->orbitalAddress+ib];
        else
          for (size_t ib=0; ib<nsb; ib++)
            w.buffer[woffset+e->stringIndex*wnsb+ib]+=
                (*this)[offa+nsa*nsb*e->orbitalAddress+ib];
      }
      offa += nsb;
    }
  } else if (deltaBeta==2) { // wavefunction has 2 more beta electrons than interacting states
    unsigned int wsyma = syma;
    unsigned int wsymb = w.symmetry^wsyma;
    size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    size_t woffset = w.blockOffset(wsyma) + wnsb * alphaStringsBegin->index(w.alphaStrings[wsyma]);
    size_t offb = 0;
    for (StringSet::const_iterator s = betaStringsBegin; s != betaStringsEnd; s++) {
      ExcitationSet ee(*s,w.betaStrings[wsymb],0,2);
      for (ExcitationSet::const_iterator e=ee.begin(); e!=ee.end(); e++) {
        if (e->phase < 0)
          for (size_t ia=0; ia<nsa; ia++)
            w.buffer[woffset+e->stringIndex+wnsb*ia]-=
                (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb];
        else
          for (size_t ia=0; ia<nsa; ia++)
            w.buffer[woffset+e->stringIndex+wnsb*ia]+=
                (*this)[offb+nsa*nsb*e->orbitalAddress+ia*nsb];
      }
      offb ++;
    }
  } else if (deltaAlpha==1 && deltaBeta==1) { // wavefunction has 1 more alpha and beta electrons than interacting states
    if (parity) throw std::logic_error("wrong parity in alpha-beta");
    for (unsigned int wsyma=0; wsyma<8; wsyma++) {
      unsigned int wsymb = w.symmetry^wsyma;
      unsigned int symexca = wsyma ^ syma;
      size_t intoff = w.orbitalSpace->offset(symexcw,symexca,parity);
      size_t wnsb = w.betaStrings[wsymb].size();
      size_t woffset = w.blockOffset(wsyma);
      size_t offb = 0;
      std::vector<ExcitationSet> eebs;
      for (StringSet::const_iterator sb = betaStringsBegin; sb != betaStringsEnd; sb++)
        eebs.push_back(ExcitationSet(*sb,w.betaStrings[wsymb],0,1));
      for (StringSet::const_iterator sa = alphaStringsBegin; sa != alphaStringsEnd; sa++) {
        ExcitationSet eea(*sa,w.alphaStrings[wsyma],0,1);
        for (std::vector<ExcitationSet>::const_iterator eebp=eebs.begin(); eebp!=eebs.end(); eebp++) {
          for (ExcitationSet::const_iterator ea=eea.begin(); ea!=eea.end(); ea++) {
            for (ExcitationSet::const_iterator eb=eebp->begin(); eb!=eebp->end(); eb++) {
              w.buffer[woffset + eb->stringIndex + wnsb * ea->stringIndex]
                  += ea->phase * eb->phase *
                  (*this)[offb + nsa*nsb*(intoff +
                                          ea->orbitalAddress +
                                          eb->orbitalAddress*(*w.orbitalSpace)[symexca]) ];
            }
          }
          offb ++;
        }
      }
    }
  } else {
    xout <<"deltaAlpha="<<deltaAlpha<<", deltaBeta="<<deltaBeta<<std::endl;
    throw std::logic_error("unimplemented case");
  }
}

#include "gciMolpro.h"
std::vector<double> TransitionDensity::density(Wavefunction &w)
{
  unsigned int syma = alphaStringsBegin->computed_symmetry();
  unsigned int symb = betaStringsBegin->computed_symmetry();
  if ((syma^symb) != w.symmetry) throw std::runtime_error("Symmetry mismatch, TransitionDensity::density");
  if (nsa*nsb != (w.blockOffset(syma+1)-w.blockOffset(syma))) throw std::range_error("Wrong dimensions, TransitionDensity::density"); // present implementation only for complete space in TransitionDensity
  std::vector<double> result(excitations,(double)0);
  MxmDrvTN(&result[0],&this->at(0),&(w.buffer[w.blockOffset(syma)]),excitations,nsa*nsb,1,1);
  return result;
}

std::string TransitionDensity::str(int verbosity, unsigned int columns) const
{
  if (verbosity < 0) return std::string("");
  std::ostringstream s;
  //  xout <<"TransitionDensity::str size()="<<size()<<std::endl;
  if (size()!= nsa*nsb*excitations) {xout << "wrong size in TransitionDensity::str "<<size()<<" "<<nsa*nsb*excitations<<std::endl; throw std::range_error("help");}
  if (size())
    for (size_t ij=0; ij < excitations; ij++) {
      s <<std::endl;
      for (size_t ab=0; ab < nsa*nsb; ab++)
        s << (*this)[ab+ij*nsa*nsb] <<" ";
    }
  return s.str();
}
