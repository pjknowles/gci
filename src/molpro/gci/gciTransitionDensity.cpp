#include "gciTransitionDensity.h"
#include <iostream>
#include <iterator>
#include <sstream>
#if defined __has_include
#if __has_include(<mkl.h>)
#include <mkl.h>
#else
#include <cblas.h>
#endif
#else
#include <cblas.h>
#endif

// using TransitionDensity = gci::TransitionDensity;
namespace molpro {
namespace gci {

TransitionDensity::TransitionDensity(const TransitionDensity &source, bool copy)
    : molpro::array<double>(source.m_nsa * source.m_nsb * source.m_excitations),
      m_alphaStringsBegin(source.m_alphaStringsBegin), m_alphaStringsEnd(source.m_alphaStringsEnd),
      m_betaStringsBegin(source.m_betaStringsBegin), m_betaStringsEnd(source.m_betaStringsEnd),
      m_parity(source.m_parity), m_nsa(source.m_nsa), m_nsb((source.m_nsb)), m_syma(source.m_syma),
      m_symb(source.m_symb), m_symexc(source.m_symexc), m_deltaAlpha(source.m_deltaAlpha),
      m_deltaBeta(source.m_deltaBeta), m_excitations(source.m_excitations) {
  if (copy)
    std::copy(source.begin(), source.end(), begin());
}
TransitionDensity::TransitionDensity(const Wavefunction &w, const StringSet::const_iterator &alphaStringsBegin,
                                     const StringSet::const_iterator &alphaStringsEnd,
                                     const StringSet::const_iterator &betaStringsBegin,
                                     const StringSet::const_iterator &betaStringsEnd, const molpro::parity_t parity,
                                     const bool doAlpha, const bool doBeta)
    : molpro::array<double>(std::distance(alphaStringsBegin, alphaStringsEnd) *
                            std::distance(betaStringsBegin, betaStringsEnd) *
                            ((w.alphaStrings[0].proto.nelec - alphaStringsBegin->nelec + w.betaStrings[0].proto.nelec -
                              betaStringsBegin->nelec) %
                                     2
                                 ? w.orbitalSpace->operator[](w.symmetry ^ alphaStringsBegin->computed_symmetry() ^
                                                              betaStringsBegin->computed_symmetry())
                                 : w.orbitalSpace->total(w.symmetry ^ alphaStringsBegin->computed_symmetry() ^
                                                             betaStringsBegin->computed_symmetry(),
                                                         parity))),
      m_alphaStringsBegin(alphaStringsBegin), m_alphaStringsEnd(alphaStringsEnd), m_betaStringsBegin(betaStringsBegin),
      m_betaStringsEnd(betaStringsEnd), m_parity(parity), m_nsa(std::distance(m_alphaStringsBegin, m_alphaStringsEnd)),
      m_nsb(std::distance(m_betaStringsBegin, m_betaStringsEnd)), m_syma(m_alphaStringsBegin->computed_symmetry()),
      m_symb(m_betaStringsBegin->computed_symmetry()), m_symexc(w.symmetry ^ m_syma ^ m_symb),
      m_deltaAlpha(w.alphaStrings[0].proto.nelec - m_alphaStringsBegin->nelec),
      m_deltaBeta(w.betaStrings[0].proto.nelec - m_betaStringsBegin->nelec),
      m_excitations((m_deltaAlpha + m_deltaBeta) % 2 ? w.orbitalSpace->operator[](m_symexc)
                                                     : w.orbitalSpace->total(m_symexc, m_parity)) {
  //  cout << "time on entering TransitionDensity " << std::chrono::steady_clock::now().time_since_epoch().count() <<"
  //  "<< gci::profiler->getResources().wall << std::endl; cout << "size(): " << memory::array<double>::size() <<
  //  std::endl; cout << m_nsa << " " << m_nsb << " " << m_excitations << " " << m_nsa * m_nsb * m_excitations <<" =
  //  "<<size() << std::endl;
  assert(molpro::array<double>::size() == m_nsa * m_nsb * m_excitations);
  if (empty())
    return;
  auto prof = profiler->push("TransitionDensity");
  {
    auto prof = profiler->push("initialise");
    prof += size();
    assign(0); // TODO needed?
  }

  if (m_deltaAlpha == 0 && m_deltaBeta == 0) { // number of electrons preserved, so one-electron excitation
    m_hasAlpha = doAlpha;
    m_hasBeta = doBeta;
    if (doAlpha) {
      auto prof2 = profiler->push("TransitionDensity_alpha");
      // alpha excitations
      const unsigned int wsymb = m_symb;
      const unsigned int wsyma = w.symmetry ^ wsymb;
      const size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      const size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
      size_t offa = 0;
      for (auto s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
        ExcitationSet ee(*s, w.alphaStrings[wsyma], 1, 1, m_parity);
        if (w.m_sparse) {
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ib = 0; ib < m_nsb; ib++)
                if (w.buffer_sparse.count(woffset + e.stringIndex * wnsb + ib))
                  (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] -=
                      w.buffer_sparse.at(woffset + e.stringIndex * wnsb + ib);
            } else {
              for (size_t ib = 0; ib < m_nsb; ib++)
                if (w.buffer_sparse.count(woffset + e.stringIndex * wnsb + ib))
                  (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] +=
                      w.buffer_sparse.at(woffset + e.stringIndex * wnsb + ib);
            }
          }

        } else {
          prof2 += ee.size() * m_nsb * 2;
          for (const auto &e : ee) {
            if (e.phase < 0)
              for (size_t ib = 0; ib < m_nsb; ib++)
                (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] -= w.buffer[woffset + e.stringIndex * wnsb + ib];
            else
              for (size_t ib = 0; ib < m_nsb; ib++)
                (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] += w.buffer[woffset + e.stringIndex * wnsb + ib];
          }
        }
        offa += m_nsb;
      }
    }

    if (doBeta) {
      auto prof2 = profiler->push("TransitionDensity_beta");
      // beta excitations
      const unsigned int wsyma = m_syma;
      const unsigned int wsymb = w.symmetry ^ wsyma;
      const size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      const size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      for (auto s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
        ExcitationSet ee(*s, w.betaStrings[wsymb], 1, 1, m_parity);
        if (w.m_sparse) {
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ia = 0; ia < m_nsa; ia++) {
                if (w.buffer_sparse.count(woffset + e.stringIndex + wnsb * ia))
                  (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] -=
                      w.buffer_sparse.at(woffset + e.stringIndex + wnsb * ia);
              }
            } else {
              for (size_t ia = 0; ia < m_nsa; ia++) {
                if (w.buffer_sparse.count(woffset + e.stringIndex + wnsb * ia))
                  (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] +=
                      w.buffer_sparse.at(woffset + e.stringIndex + wnsb * ia);
              }
            }
          }
        } else {
          prof2 += ee.size() * m_nsa * 2;
          auto v = &((*this)[offb]);
          for (const auto &e : ee) {
            if (e.phase < 0)
              for (size_t ia = 0; ia < m_nsa; ia++) {
                *(v + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb) -= w.buffer[woffset + e.stringIndex + wnsb * ia];
                //                    (*this)[offb+nsa*nsb*e.orbitalAddress+ia*nsb] -=
                //                    w.buffer[woffset+e.stringIndex+wnsb*ia];
              }
            else
              for (size_t ia = 0; ia < m_nsa; ia++) {
                *(v + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb) += w.buffer[woffset + e.stringIndex + wnsb * ia];
                //                    (*this)[offb+nsa*nsb*e.orbitalAddress+ia*nsb] +=
                //                    w.buffer[woffset+e.stringIndex+wnsb*ia];
              }
          }
        }
        offb++;
      }
    }

  } else if (m_deltaBeta == 0 &&
             m_deltaAlpha > 0) { // wavefunction has deltaAlpha more alpha electrons than interacting states
    auto prof2 = profiler->push("TransitionDensity_alpha_alpha");
    m_hasAlpha = true;
    m_hasBeta = false;
    const unsigned int wsymb = m_symb;
    const unsigned int wsyma = w.symmetry ^ wsymb;
    const size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    const size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
    //    cout << "alpha wsyma="<<wsyma<<", wsymb="<<wsymb<<", woffset="<<woffset<<std::endl;

    size_t offa = 0;
    for (auto s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
      ExcitationSet ee(*s, w.alphaStrings[wsyma], 0, m_deltaAlpha, m_parity);
      if (w.m_sparse) {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ib = 0; ib < m_nsb; ib++)
              if (w.buffer_sparse.count(woffset + e.stringIndex * wnsb + ib))
                (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] -=
                    w.buffer_sparse.at(woffset + e.stringIndex * wnsb + ib);
          } else {
            for (size_t ib = 0; ib < m_nsb; ib++)
              if (w.buffer_sparse.count(woffset + e.stringIndex * wnsb + ib))
                (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] +=
                    w.buffer_sparse.at(woffset + e.stringIndex * wnsb + ib);
          }
        }
      } else {
        prof2 += ee.size() * m_nsb * 2;
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ib = 0; ib < m_nsb; ib++)
              (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] -= w.buffer[woffset + e.stringIndex * wnsb + ib];
          } else {
            for (size_t ib = 0; ib < m_nsb; ib++)
              (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] += w.buffer[woffset + e.stringIndex * wnsb + ib];
          }
        }
      }
      offa += m_nsb;
    }
  } else if (m_deltaAlpha == 0 &&
             m_deltaBeta > 0) { // wavefunction has deltaBeta more beta electrons than interacting states
    auto prof2 = profiler->push("TransitionDensity_beta_beta");
    m_hasAlpha = false;
    m_hasBeta = true;
    const unsigned int wsyma = m_syma;
    const unsigned int wsymb = w.symmetry ^ wsyma;
    const size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    const size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
    size_t offb = 0;
    for (auto s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
      ExcitationSet ee(*s, w.betaStrings[wsymb], 0, m_deltaBeta, m_parity);
      if (w.m_sparse) {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ia = 0; ia < m_nsa; ia++)
              if (w.buffer_sparse.count(woffset + e.stringIndex + wnsb * ia))
                (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] -=
                    w.buffer_sparse.at(woffset + e.stringIndex + wnsb * ia);
          } else {
            for (size_t ia = 0; ia < m_nsa; ia++)
              if (w.buffer_sparse.count(woffset + e.stringIndex + wnsb * ia))
                (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] +=
                    w.buffer_sparse.at(woffset + e.stringIndex + wnsb * ia);
          }
        }
      } else {
        prof2 += ee.size() * m_nsa * 2;
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ia = 0; ia < m_nsa; ia++)
              (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] -=
                  w.buffer[woffset + e.stringIndex + wnsb * ia];
          } else {
            for (size_t ia = 0; ia < m_nsa; ia++)
              (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] +=
                  w.buffer[woffset + e.stringIndex + wnsb * ia];
          }
        }
      }
      offb++;
    }
  } else if (m_deltaAlpha == 1 &&
             m_deltaBeta == 1) { // wavefunction has 1 more alpha and beta electrons than interacting states
    auto prof2 = profiler->push("TransitionDensity_alpha_beta");
    if (m_parity != molpro::parityNone)
      throw std::logic_error("wrong parity in alpha-beta");
    m_hasAlpha = true;
    m_hasBeta = true;
    for (unsigned int wsyma = 0; wsyma < 8; wsyma++) {
      unsigned int wsymb = w.symmetry ^ wsyma;
      unsigned int symexca = wsyma ^ m_syma;
      const size_t wnsb = w.betaStrings[wsymb].size();
      if (wnsb == 0)
        continue;
      const size_t woffset = w.blockOffset(wsyma);
      size_t offb = m_nsa * m_nsb * w.orbitalSpace->offset(m_symexc, symexca, m_parity);
      std::vector<ExcitationSet> eebs;
      for (auto sb = m_betaStringsBegin; sb != m_betaStringsEnd; sb++) {
        eebs.emplace_back(*sb, w.betaStrings[wsymb], 0, 1, m_parity);
        for (auto &eb : eebs.back()) // optimisation for non-standard use in following loop
          const_cast<Excitation &>(eb).orbitalAddress *= m_nsa * m_nsb * (*w.orbitalSpace)[symexca];
      }
      for (auto sa = m_alphaStringsBegin; sa != m_alphaStringsEnd; sa++) {
        ExcitationSet eea(*sa, w.alphaStrings[wsyma], 0, 1, m_parity);
        for (auto &ea : eea) { // optimisation for non-standard use in following loop
          const_cast<Excitation &>(ea).orbitalAddress *= m_nsa * m_nsb;
          const_cast<Excitation &>(ea).stringIndex *= wnsb;
        }
        if (w.m_sparse) {
          for (const auto &eeb : eebs) {
            for (const auto &eb : eeb) {
              if (eb.phase > 0) {
                for (const auto &ea : eea)
                  if (w.buffer_sparse.count(woffset + eb.stringIndex + ea.stringIndex))
                    (*this)[offb + eb.orbitalAddress + ea.orbitalAddress] +=
                        ea.phase * w.buffer_sparse.at(woffset + eb.stringIndex + ea.stringIndex);
              } else {
                for (const auto &ea : eea)
                  if (w.buffer_sparse.count(woffset + eb.stringIndex + ea.stringIndex))
                    (*this)[offb + eb.orbitalAddress + ea.orbitalAddress] -=
                        ea.phase * w.buffer_sparse.at(woffset + eb.stringIndex + ea.stringIndex);
              }
            }
            offb++;
          }
        } else {
          for (const auto &eeb : eebs) {
            prof2 += eea.size() * eeb.size() * 2;
            for (const auto &eb : eeb) {
              if (eb.phase > 0) {
                for (const auto &ea : eea)
                  (*this)[offb + eb.orbitalAddress + ea.orbitalAddress] +=
                      ea.phase * w.buffer[woffset + eb.stringIndex + ea.stringIndex];
              } else {
                for (const auto &ea : eea)
                  (*this)[offb + eb.orbitalAddress + ea.orbitalAddress] -=
                      ea.phase * w.buffer[woffset + eb.stringIndex + ea.stringIndex];
              }
            }
            offb++;
          }
        }
      }
    }
  } else {
    cout << "deltaAlpha=" << m_deltaAlpha << ", deltaBeta=" << m_deltaBeta << std::endl;
    throw std::logic_error("unimplemented case");
  }
  //  size_t populated=0; for (const_iterator x=begin(); x!=end(); x++) if (*x !=(double)0) ++populated; cout
  //  <<"TransitionDensity population="<<((double)populated)/((double)size()+1)<<std::endl;
  //  cout << "time on leaving TransitionDensity " << std::chrono::steady_clock::now().time_since_epoch().count() <<"
  //  "<< gci::profiler->getResources().wall << std::endl;
}

void TransitionDensity::action(Wavefunction &w) const {
  if (m_nsa * m_nsb * m_excitations == 0)
    return;
  auto prof = profiler->push("TransitionDensity::action");

  if (m_deltaAlpha == 0 && m_deltaBeta == 0) { // number of electrons preserved, so one-electron excitation

    if (true) {
      // alpha excitations
      unsigned int wsymb = m_symb;
      unsigned int wsyma = w.symmetry ^ wsymb;
      size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
      size_t offa = 0;
      for (auto s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
        ExcitationSet ee(*s, w.alphaStrings[wsyma], 1, 1, m_parity);
        if (w.m_sparse) {
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ib = 0; ib < m_nsb; ib++)
                if ((*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] != 0)
                  w.buffer_sparse[woffset + e.stringIndex * wnsb + ib] -=
                      (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
            } else {
              for (size_t ib = 0; ib < m_nsb; ib++)
                if ((*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] != 0)
                  w.buffer_sparse[woffset + e.stringIndex * wnsb + ib] +=
                      (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
            }
          }
        } else {
          prof += ee.size() * m_nsb * 2;
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ib = 0; ib < m_nsb; ib++)
                w.buffer[woffset + e.stringIndex * wnsb + ib] -= (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
            } else {
              for (size_t ib = 0; ib < m_nsb; ib++)
                w.buffer[woffset + e.stringIndex * wnsb + ib] += (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
            }
          }
        }
        offa += m_nsb;
      }
    }

    if (true) {
      // beta excitations
      const unsigned int wsyma = m_syma;
      const unsigned int wsymb = w.symmetry ^ wsyma;
      const size_t wnsb = w.betaStrings[wsymb].size();
      // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
      const size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
      size_t offb = 0;
      for (auto s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
        ExcitationSet ee(*s, w.betaStrings[wsymb], 1, 1, m_parity);
        if (w.m_sparse) {
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ia = 0; ia < m_nsa; ia++)
                if ((*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] != 0)
                  w.buffer_sparse[woffset + e.stringIndex + wnsb * ia] -=
                      (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
            } else {
              for (size_t ia = 0; ia < m_nsa; ia++)
                if ((*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] != 0)
                  w.buffer_sparse[woffset + e.stringIndex + wnsb * ia] +=
                      (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
            }
          }
        } else {
          prof += ee.size() * m_nsa * 2;
          for (const auto &e : ee) {
            if (e.phase < 0) {
              for (size_t ia = 0; ia < m_nsa; ia++)
                w.buffer[woffset + e.stringIndex + wnsb * ia] -=
                    (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
            } else {
              for (size_t ia = 0; ia < m_nsa; ia++)
                w.buffer[woffset + e.stringIndex + wnsb * ia] +=
                    (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
            }
          }
        }
        offb++;
      }
    }

  } else if (m_deltaBeta == 0 &&
             m_deltaAlpha > 0) { // wavefunction has deltaAlpha more alpha electrons than interacting states
    const unsigned int wsymb = m_symb;
    const unsigned int wsyma = w.symmetry ^ wsymb;
    const size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    const size_t woffset = w.blockOffset(wsyma) + m_betaStringsBegin->index(w.betaStrings[wsymb]);
    size_t offa = 0;
    for (auto s = m_alphaStringsBegin; s != m_alphaStringsEnd; s++) {
      ExcitationSet ee(*s, w.alphaStrings[wsyma], 0, m_deltaAlpha, m_parity);
      if (w.m_sparse) {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ib = 0; ib < m_nsb; ib++)
              if ((*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] != 0)
                w.buffer_sparse[woffset + e.stringIndex * wnsb + ib] -=
                    (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
          } else {
            for (size_t ib = 0; ib < m_nsb; ib++)
              if ((*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib] != 0)
                w.buffer_sparse[woffset + e.stringIndex * wnsb + ib] +=
                    (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
          }
        }
      } else {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ib = 0; ib < m_nsb; ib++)
              w.buffer[woffset + e.stringIndex * wnsb + ib] -= (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
          } else {
            for (size_t ib = 0; ib < m_nsb; ib++)
              w.buffer[woffset + e.stringIndex * wnsb + ib] += (*this)[offa + m_nsa * m_nsb * e.orbitalAddress + ib];
          }
        }
      }
      offa += m_nsb;
    }
  } else if (m_deltaAlpha == 0 &&
             m_deltaBeta > 0) { // wavefunction has deltaBeta more beta electrons than interacting states
    const unsigned int wsyma = m_syma;
    const unsigned int wsymb = w.symmetry ^ wsyma;
    const size_t wnsb = w.betaStrings[wsymb].size();
    // assumes that alphaStrings, betaStrings are contiguous ordered subsets of wavefunction strings
    const size_t woffset = w.blockOffset(wsyma) + wnsb * m_alphaStringsBegin->index(w.alphaStrings[wsyma]);
    size_t offb = 0;
    for (auto s = m_betaStringsBegin; s != m_betaStringsEnd; s++) {
      ExcitationSet ee(*s, w.betaStrings[wsymb], 0, m_deltaBeta, m_parity);
      if (w.m_sparse) {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ia = 0; ia < m_nsa; ia++)
              if ((*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] != 0)
                w.buffer_sparse[woffset + e.stringIndex + wnsb * ia] -=
                    (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
          } else {
            for (size_t ia = 0; ia < m_nsa; ia++)
              if ((*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb] != 0)
                w.buffer_sparse[woffset + e.stringIndex + wnsb * ia] +=
                    (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
          }
        }
      } else {
        for (const auto &e : ee) {
          if (e.phase < 0) {
            for (size_t ia = 0; ia < m_nsa; ia++)
              w.buffer[woffset + e.stringIndex + wnsb * ia] -=
                  (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
          } else {
            for (size_t ia = 0; ia < m_nsa; ia++)
              w.buffer[woffset + e.stringIndex + wnsb * ia] +=
                  (*this)[offb + m_nsa * m_nsb * e.orbitalAddress + ia * m_nsb];
          }
        }
      }
      offb++;
    }
  } else if (m_deltaAlpha == 1 &&
             m_deltaBeta == 1) { // wavefunction has 1 more alpha and beta electrons than interacting states
    if (m_parity != molpro::parityNone)
      throw std::logic_error("wrong parity in alpha-beta");
    for (unsigned int wsyma = 0; wsyma < 8; wsyma++) {
      unsigned int wsymb = w.symmetry ^ wsyma;
      unsigned int symexca = wsyma ^ m_syma;
      const size_t wnsb = w.betaStrings[wsymb].size();
      const size_t woffset = w.blockOffset(wsyma);
      size_t offb = m_nsa * m_nsb * w.orbitalSpace->offset(m_symexc, symexca, m_parity);
      std::vector<ExcitationSet> eebs;
      for (auto sb = m_betaStringsBegin; sb != m_betaStringsEnd; sb++) {
        eebs.emplace_back(*sb, w.betaStrings[wsymb], 0, 1);
        for (auto &eb : eebs.back()) // optimisation for non-standard use in following loop
          const_cast<Excitation &>(eb).orbitalAddress *= m_nsa * m_nsb * (*w.orbitalSpace)[symexca];
      }
      for (auto sa = m_alphaStringsBegin; sa != m_alphaStringsEnd; sa++) {
        ExcitationSet eea(*sa, w.alphaStrings[wsyma], 0, 1);
        for (auto &ea : eea) { // optimisation for non-standard use in following loop
          const_cast<Excitation &>(ea).orbitalAddress *= m_nsa * m_nsb;
          const_cast<Excitation &>(ea).stringIndex *= wnsb;
        }
        if (w.m_sparse) {
          for (const auto &eeb : eebs) {
            prof += eea.size() * eeb.size() * 2;
            for (const auto &eb : eeb) {
              if (eb.phase > 0)
                for (const auto &ea : eea) {
                  if ((*this)[offb + eb.orbitalAddress + ea.orbitalAddress] != 0)
                    w.buffer_sparse[woffset + eb.stringIndex + ea.stringIndex] +=
                        ea.phase * (*this)[offb + eb.orbitalAddress + ea.orbitalAddress];
                }
              else
                for (const auto &ea : eea) {
                  if ((*this)[offb + eb.orbitalAddress + ea.orbitalAddress] != 0)
                    w.buffer_sparse[woffset + eb.stringIndex + ea.stringIndex] -=
                        ea.phase * (*this)[offb + eb.orbitalAddress + ea.orbitalAddress];
                }
            }
            offb++;
          }
        } else {
          for (const auto &eeb : eebs) {
            prof += eea.size() * eeb.size() * 2;
            for (const auto &eb : eeb) {
              if (eb.phase > 0)
                for (const auto &ea : eea) {
                  w.buffer[woffset + eb.stringIndex + ea.stringIndex] +=
                      ea.phase * (*this)[offb + eb.orbitalAddress + ea.orbitalAddress];
                }
              else
                for (const auto &ea : eea) {
                  w.buffer[woffset + eb.stringIndex + ea.stringIndex] -=
                      ea.phase * (*this)[offb + eb.orbitalAddress + ea.orbitalAddress];
                }
            }
            offb++;
          }
        }
      }
    }
  } else {
    cout << "deltaAlpha=" << m_deltaAlpha << ", deltaBeta=" << m_deltaBeta << std::endl;
    throw std::logic_error("unimplemented case");
  }
}

molpro::Operator TransitionDensity::density(const Wavefunction &w) const {
  molpro::dim_t dimension(8);
  for (unsigned int i = 0; i < 8; i++) {
    dimension[i] = 0;
    for (const auto &s : w.orbitalSpace->orbital_symmetries)
      if (s == i)
        ++dimension[i];
  }
  molpro::Operator result(dimension, 1, !(m_hasAlpha && m_hasBeta), m_symexc);

  unsigned int syma = m_alphaStringsBegin->computed_symmetry();
  unsigned int symb = m_betaStringsBegin->computed_symmetry();
  if ((syma ^ symb) != w.symmetry)
    throw std::runtime_error("Symmetry mismatch, TransitionDensity::density");
  if (m_nsa * m_nsb != (w.blockOffset(syma + 1) - w.blockOffset(syma)))
    throw std::range_error("Wrong dimensions, TransitionDensity::density"); // present implementation only for complete
                                                                            // space in TransitionDensity
  double *out = &((*result.O1(true).data())[0]);
  const double *a = &at(0);
  const double *b = &(w.buffer[w.blockOffset(syma)]);
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m_excitations, 1, m_nsa * m_nsb, 1, a, 1, b, m_nsa * m_nsb,
              false ? 1 : 0, out, m_excitations);
  return result;
}

std::string TransitionDensity::str(int verbosity, unsigned int columns) const {
  if (verbosity < 0)
    return std::string("");
  std::ostringstream s;
  //  cout <<"TransitionDensity::str size()="<<size()<<std::endl;
  if (size() != m_nsa * m_nsb * m_excitations) {
    cout << "wrong size in TransitionDensity::str " << size() << " " << m_nsa * m_nsb * m_excitations << std::endl;
    throw std::range_error("help");
  }
  if (!empty())
    for (size_t ij = 0; ij < m_excitations; ij++) {
      s << std::endl;
      for (size_t ab = 0; ab < m_nsa * m_nsb; ab++)
        s << (*this)[ab + ij * m_nsa * m_nsb] << " ";
    }
  return s.str();
}
} // namespace gci
} // namespace molpro
