#ifndef GCI_GCIHPRODUCT_H
#define GCI_GCIHPRODUCT_H

#include <vector>

#include "gciVibSpace.h"

namespace gci {

/*!
 * @brief Hartree product is a single vibrational basis represented as occupation numbers for excited states.
 *
 * Excitations are stored as {{modeIndex, modalIndex}, ...}:
 *      - modeIndex  -- mode outside of its ground state
 *      - modalIndex -- index for the occupied excitated state basis function
 *
 * Order:
 *      - only one basis per mode (no modeIndex occurs twice)
 *      - modeIndex in increasing order
 *      - modalIndex != 0, since this is the assumed value for the ground state.
 */
class HProduct {
public:
    using t_Modal = std::vector<int>;
    using t_Product = std::vector<t_Modal>;
    HProduct(t_Product phi);

    t_Product::const_iterator begin() const {return m_prod.cbegin();}
    t_Product::const_iterator end() const {return m_prod.cend();}

    bool operator==(const HProduct &other) const {return m_prod == other.m_prod;}
    bool operator!=(const HProduct &other) const {return !(*this == other);}

    //! Checks if the product is empty. This implies ground state.
    bool empty() const {return m_prod.empty();}

    //! Number of excited modes
    auto excLvl() const {return m_prod.size();}

    //! Excitation level of a mode
    int excLvl(int iMode) const;

    const auto &operator[](unsigned long i) const {return m_prod[i];}

    /*!
     * @brief Excites a mode by one level.
     * @note No check is made whether the resultant product is still within some desired Fock space
     * @param iMode Index of the mode to be excited
     * @return Excited HProduct
     */
    void raise(int iMode) {changeModal(iMode, +1);}

    /*!
     * @brief Loweres a mode by one level.
     * @note No check is made whether the resultant product is still within some desired Fock space
     * @param iMode Index of the mode to be excited
     * @return Excited HProduct
     */
    void lower(int iMode) {changeModal(iMode, -1);}

    /*!
     * @brief Changes occupied modal by `diff` excitations.
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
    bool withinSpace(const VibSpace &vibSpace);

    /*!
     * @brief List of excited modes in the current product
     */
    std::vector<int> excitedModes() const;
protected:
    t_Product m_prod; //!< excitations representing the Hartree product

    //! Enforces ordering of the HartreeProduct..
    void reorder(t_Product &prod) const;

    //! Checks that the product makes sense. @warning Assumes that product is ordered
    void check(t_Product &prod) const;

};
}//  namespace gci

#include <vector>

#endif //GCI_GCIHPRODUCT_H
