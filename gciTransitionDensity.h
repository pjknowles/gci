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
  /*!
   * \brief Construct a TransitionDensity from a wavefunction to a subset of space defined by sets of alpha and beta String objects
   * \param w
   * \param alphaStrings
   * \param betaStrings
   * \param parity
   * \param doAlpha whether to process alpha excitations
   * \param doBeta whether to process beta excitations
   */
  TransitionDensity(const Wavefunction& w, const StringSet &alphaStrings, const StringSet &betaStrings, int parity, const bool doAlpha=true, const bool doBeta=true);
  /*!
   * \brief Collapse onto a configuration-space residual
   * g(I) += E(K,exc) <I|exc|K>
   * \param g a Wavefunction object that will receive contributions
   */
  void action(Wavefunction& g);
  std::string str(int verbosity) const;
private:
  size_t nsa; ///< number of alpha strings
  size_t nsb; ///< number of beta strings
  size_t excitations; ///< number of excitations
  const StringSet* alphaStrings;
  const StringSet* betaStrings;
  int parity;
};
}

using namespace gci;

#endif // GCITRANSITIONDENSITY_H
