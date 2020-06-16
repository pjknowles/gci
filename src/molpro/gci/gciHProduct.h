#ifndef GCI_GCIHPRODUCT_H
#define GCI_GCIHPRODUCT_H

#include <ostream>

#include "molpro/gci/gciVibExcitation.h"
#include "molpro/gci/gciVibSpace.h"

namespace molpro {
namespace gci {

/*!
 * @brief Hartree product is a single vibrational basis represented as excitations out of the ground state.
 *
 * Excitations are stored as {{modeIndex, modalIndex}, ...}:
 *      - modeIndex  -- index of excited mode
 *      - modalIndex -- index of the excitated modal
 *
 * Order:
 */
class HProduct {
public:
  using modal_t = std::array<int, 2>;
  using product_t = std::vector<modal_t>;

protected:
  product_t m_prod; //!< excitations representing the Hartree product
public:
  /*!
   * @brief Constructs an empty Hartree product representing the ground state.
   * @warning This is the preferred way to construct an empty product.
   */
  HProduct() : m_prod(product_t{}) {}
  HProduct(product_t phi);

  product_t::iterator begin() { return m_prod.begin(); }
  product_t::iterator end() { return m_prod.end(); }
  product_t::const_iterator begin() const { return m_prod.begin(); }
  product_t::const_iterator end() const { return m_prod.end(); }
  product_t::const_iterator cbegin() const { return m_prod.cbegin(); }
  product_t::const_iterator cend() const { return m_prod.cend(); }

  bool operator==(const HProduct &other) const { return m_prod == other.m_prod; }
  bool operator!=(const HProduct &other) const { return !(*this == other); }

  //! Checks if the product is empty. This implies ground state.
  bool empty() const { return m_prod.empty(); }

  //! Number of excited modes
  auto excLvl() const { return m_prod.size(); }

  //! Excitation level of a mode
  int excLvl(int iMode) const;

  const auto &operator[](unsigned long i) const { return m_prod[i]; }

  /*!
   * @brief Applies the vibrational excitation operator
   * @note If annihilation operator does not match the occupied modal, that modal is set to -1 and will be declared
   * outside of vibrational space when checked.
   */
  HProduct excite(const VibExcitation &exc) const;

  /*!
   * @brief Excites a mode by one level.
   * @note No check is made whether the resultant product is still within some desired Fock space
   * @param iMode Index of the mode to be excited
   * @return Excited HProduct
   */
  void raise(int iMode) { changeModal(iMode, +1); }

  /*!
   * @brief Loweres a mode by one level.
   * @note No check is made whether the resultant product is still within some desired Fock space
   * @param iMode Index of the mode to be excited
   * @return Excited HProduct
   */
  void lower(int iMode) { changeModal(iMode, -1); }

  /*!
   * @brief Changes occupied modal by `diff` excitations.
   * @note No check is made over legitemacy of the resulted Product. An explicit check whether result is withinSpace
   * has to be made.
   * @param iMode Mode to be changed
   * @param diff Relative number of excitations(de-excitations positive/negative values
   * @return Excited Hartree product
   */
  void changeModal(int iMode, int diff);

  /*!
   * @brief Checks that the product is within a specified space
   * @param vibSpace Vibrational Fock space
   * @return true if it is, false if it's not
   */
  bool withinSpace(const VibSpace &vibSpace) const;

  //! @brief List of excited modes in `this` product
  std::vector<int> excitedModes() const;

protected:
  /*!
   * @brief Orders `this` Hartree product
   *      - only one basis per mode (no modeIndex occurs twice)
   *      - modeIndex in increasing order
   *      - modalIndex != 0, since this is the assumed value for the ground state.
   * @param prod
   */
  void order();

  /*!
   * @brief Checks that `this` product is legitimate.
   * @warning Assumes that the product is ordered
   */
  void check() const;
};

std::ostream &operator<<(std::ostream &os, HProduct const &obj);

} //  namespace gci
} // namespace molpro

#endif // GCI_GCIHPRODUCT_H
