#ifndef GCIORBITALS_H
#define GCIORBITALS_H
#include "SMat.h"
#include "gciOrbitalSpace.h"
#include "gciPrintable.h"
#include <memory>
#include <iostream>
#include <numeric>


class Orbitals : public Printable
{
public:
  Orbitals(const OrbitalSpace& orbitalSpace)
    : m_orbitals(dims_t{orbitalSpace,orbitalSpace},parityNone,0,"orbitals")
    , m_energies(dims_t{orbitalSpace},parityNone,-1,"energies")
    , m_occupations(dims_t{orbitalSpace},parityNone,-1,"occupations")
  {  }
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const {
    std::ostringstream s;
    if (std::accumulate(m_orbitals.data()->begin(),m_orbitals.data()->end(),0)!=0)
      s << m_orbitals.str(m_orbitals.m_description,verbosity)<<std::endl;
    if (std::accumulate(m_energies.data()->begin(),m_energies.data()->end(),0)!=0)
      s << m_energies.str(m_energies.m_description,verbosity)<<std::endl;
    if (std::accumulate(m_occupations.data()->begin(),m_occupations.data()->end(),0)!=0)
      s << m_occupations.str(m_occupations.m_description,verbosity)<<std::endl;
    return s.str();
  }
private:
  Orbitals();
public:
  SMat m_orbitals;
  SMat m_energies;
  SMat m_occupations;
};

#endif // GCIORBITALS_H
