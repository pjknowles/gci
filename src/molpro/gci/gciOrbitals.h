#ifndef GCIORBITALS_H
#define GCIORBITALS_H
#include "molpro/gci/gciOrbitalSpace.h"
#include "molpro/gci/gciPrintable.h"
#include <iostream>
#include <memory>
#include <molpro/symmetry_matrix/SMat.h>
#include <numeric>

namespace molpro {
namespace gci {
class Orbitals : public Printable {
public:
  Orbitals(const OrbitalSpace &orbitalSpace)
      : m_orbitals(molpro::dims_t{orbitalSpace, orbitalSpace}, molpro::parityNone, 0, "orbitals"),
        m_energies(molpro::dims_t{orbitalSpace}, molpro::parityNone, -1, "energies"),
        m_occupations(molpro::dims_t{orbitalSpace}, molpro::parityNone, -1, "occupations") {}
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const {
    std::ostringstream s;
    if (std::accumulate(m_orbitals.data()->begin(), m_orbitals.data()->end(), 0) != 0)
      s << m_orbitals.str(m_orbitals.m_description, verbosity) << std::endl;
    if (std::accumulate(m_energies.data()->begin(), m_energies.data()->end(), 0) != 0)
      s << m_energies.str(m_energies.m_description, verbosity) << std::endl;
    if (std::accumulate(m_occupations.data()->begin(), m_occupations.data()->end(), 0) != 0)
      s << m_occupations.str(m_occupations.m_description, verbosity) << std::endl;
    return s.str();
  }

private:
  Orbitals();

public:
  molpro::SMat m_orbitals;
  molpro::SMat m_energies;
  molpro::SMat m_occupations;
};
} // namespace gci
} // namespace molpro

#endif // GCIORBITALS_H
