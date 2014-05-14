#ifndef GCITRANSITIONDENSITY_H
#define GCITRANSITIONDENSITY_H
#include <vector>
#include <gciStringSet.h>
#include <gciOrbitalSpace.h>
#include <gciWavefunction.h>
#include <gciPrintable.h>

namespace gci {

/*!
 * \brief Class to hold transition density matrix,
 * defined by an array of ExcitationSet objects
 */
class TransitionDensity : public std::vector<double>, public Printable
{
public:
  TransitionDensity(const Wavefunction& w, const StringSet& alphaStrings, const StringSet& betaStrings, int parity);
  std::string str(int verbosity) const;
private:
  size_t nsa; ///< number of alpha strings
  size_t nsb; ///< number of beta strings
  size_t excitations; ///< number of excitations
};
}

using namespace gci;

#endif // GCITRANSITIONDENSITY_H
