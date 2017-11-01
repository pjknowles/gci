#ifndef GCITRANSITIONDENSITY_H
#define GCITRANSITIONDENSITY_H
#include "memory.h"
#include <vector>
#include "gciStringSet.h"
#include "gciOrbitalSpace.h"
#include "Operator.h"
#include "gciWavefunction.h"
#include "gciPrintable.h"

namespace gci {

/*!
 * \brief Class to hold transition density matrix,
 * defined by an array of ExcitationSet objects
 */
class TransitionDensity : public memory::vector<double>, public Printable
{
public:
  /*!
   * \brief Construct a TransitionDensity from a wavefunction to a subset of space defined by sets of alpha and beta String objects
   * \param w
   * \param alphaStringsBegin
   * \param alphaStringsEnd
   * \param betaStringsBegin
   * \param betaStringsEnd
   * \param parity
   * \param doAlpha whether to process alpha excitations
   * \param doBeta whether to process beta excitations
   */
  TransitionDensity(const Wavefunction& w,
                    const StringSet::const_iterator &alphaStringsBegin,
                    const StringSet::const_iterator &alphaStringsEnd,
                    const StringSet::const_iterator &betaStringsBegin,
                    const StringSet::const_iterator &betaStringsEnd,
                    const int parity, const bool doAlpha=true, const bool doBeta=true);
  /*!
   * \brief Collapse onto a configuration-space residual
   * w(I) += E(K,exc) <I|exc|K>
   * \param w a Wavefunction object that will receive contributions
   */
  void action(Wavefunction& w);
  /*!
   * \brief Contract TransitionDensity with a bra state to form a 1-particle density matrix
   * \param w the bra state
   * \return  the 1-particle density matrix
   */
  std::vector<double> density(Wavefunction& w);
  SymmetryMatrix::Operator density(const Wavefunction &w);

  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;
private:
  size_t nsa; ///< number of alpha strings
  size_t nsb; ///< number of beta strings
  unsigned int symexc; ///< symmetry of excitations
  size_t excitations; ///< number of excitations
  const StringSet::const_iterator alphaStringsBegin;
  const StringSet::const_iterator alphaStringsEnd;
  const StringSet::const_iterator betaStringsBegin;
  const StringSet::const_iterator betaStringsEnd;
  const int parity;
  bool m_hasAlpha;
  bool m_hasBeta;
};
}

using namespace gci;

#endif // GCITRANSITIONDENSITY_H
