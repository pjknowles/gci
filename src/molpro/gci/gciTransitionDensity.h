#ifndef GCITRANSITIONDENSITY_H
#define GCITRANSITIONDENSITY_H
#include <molpro/memory.h>
#include <molpro/Operator.h>
#include <molpro/SMat.h>
#include <vector>

#include "molpro/gci/gciStringSet.h"
#include "molpro/gci/gciOrbitalSpace.h"
#include "molpro/gci/gciWavefunction.h"
#include "molpro/gci/gciPrintable.h"

namespace gci {

/*!
 * \brief Class to hold transition density matrix,
 * defined by an array of ExcitationSet objects
 */
template <typename T, typename A = std::allocator<T>>
class ctor_allocator : public A
{
  using a_t = std::allocator_traits<A>;
 public:
  using A::A; // Inherit constructors from A

  template <typename U> struct rebind
  {
    using other =
    ctor_allocator
        <  U, typename a_t::template rebind_alloc<U>  >;
  };

  template <typename U>
  void construct(U* ptr)
  noexcept(std::is_nothrow_default_constructible<U>::value)
  {
    ::new(static_cast<void*>(ptr)) U;
  }

  template <typename U, typename...Args>
  void construct(U* ptr, Args&&... args)
  {
    a_t::construct(static_cast<A&>(*this),
                   ptr, std::forward<Args>(args)...);
  }
};
//class TransitionDensity : public memory::vector<double, memory::allocator<double> >, public Printable {
 class TransitionDensity : public molpro::array<double>, public Printable {
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
  TransitionDensity(const Wavefunction &w,
                    const StringSet::const_iterator &alphaStringsBegin,
                    const StringSet::const_iterator &alphaStringsEnd,
                    const StringSet::const_iterator &betaStringsBegin,
                    const StringSet::const_iterator &betaStringsEnd,
                    molpro::parity_t parity, bool doAlpha = true, bool doBeta = true);

  /*!
   * @brief Copy constructor, with the option not to copy data
   * @param source
   * @param copy If false, do not contents of source
   */
  TransitionDensity(const TransitionDensity& source, bool copy=true);

  /*!
   * \brief Collapse onto a configuration-space residual
   * w(I) += E(K,exc) <I|exc|K>
   * \param w a Wavefunction object that will receive contributions
   */
  void action(Wavefunction &w) const;
  /*!
   * \brief Contract TransitionDensity with a bra state to form a 1-particle density matrix
   * \param w the bra state
   * \return  the 1-particle density matrix
   */
  molpro::Operator density(const Wavefunction &w) const;

  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const override;
 private:
  const StringSet::const_iterator m_alphaStringsBegin;
  const StringSet::const_iterator m_alphaStringsEnd;
  const StringSet::const_iterator m_betaStringsBegin;
  const StringSet::const_iterator m_betaStringsEnd;
  const molpro::parity_t m_parity;
  const size_t m_nsa; ///< number of alpha strings
  const size_t m_nsb; ///< number of beta strings
  const unsigned int m_syma;
  const unsigned int m_symb;
  const unsigned int m_symexc; ///< symmetry of excitations
  const int m_deltaAlpha;
  const int m_deltaBeta;
  const size_t m_excitations; ///< number of excitations
  bool m_hasAlpha;
  bool m_hasBeta;
};
}

#endif // GCITRANSITIONDENSITY_H
