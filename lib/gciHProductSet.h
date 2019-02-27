#ifndef GCI_GCIHPRODUCTSET_H
#define GCI_GCIHPRODUCTSET_H

#include <vector>

#include "gciHProduct.h"
#include "gciVibSpace.h"
#include "gciMixedOperator.h"

namespace gci {

/*!
 * @brief Collection of Hartree products.
 */
class HProductSet {
public:
    using const_iterator = std::vector<HProduct>::const_iterator;
    using value_type = HProduct;
    //! Constructs the complete vibrational Fock space
    HProductSet(const VibSpace &vibSpace)
            : m_vibSpace(vibSpace), m_basis(), m_vibDim(0), m_excLvlDim(vibSpace.excLvl + 1, 0),
              m_connectedSet(false) {setVibDim();};

    /*!
     * @brief Generates a set of products connected to the `bra` through an operator. This is a subset of the full space.
     *
     * In ascending order with respect to the excitation level and mode index.
     * @param vibSpace Definition of the full vibrational Fock space
     * @param bra Reference product
     * @param vibOp Operator connecting reference to products in this set
     */
    HProductSet(const HProduct &bra, const VibSpace &vibSpace, const VibOp &vibOp);

    const_iterator begin() const {return m_basis.cbegin();}
    const_iterator end() const {return m_basis.cend();}
    const auto &operator[](unsigned long i) const {return m_basis[i];}

    /*!
     * @brief Generates the full vibrational basis
     * @warning Cannot be called `if (m_connectedSet)`
     *
     * Ordering : {GS, S, D, T, Q, ...}
     * with GS -- ground state, S -- single excitations, D -- double excitations etc.
     * S are ordered in ascending mode index
     * D, T, Q etc are stored in a lower triangular format
     * E.g. Loop structure for D,
     *      for ( i = 0; i < nMode; ++i)
     *          for (j = 0; j < i; ++j)
     */
    void generateFullSpace();

    /*!
     * @brief Returns index to a basis function within the full space
     * @warning Cannot be called `if (connectedSet)`
     * @param phi Vibrational basis function
     */
    size_t index(const HProduct &phi) const;

    //! Total size of the vibrational space
    auto vibDim() const {return m_vibDim;}

    //! Size of the basis stored in `this`
    auto size() const {return m_basis.size();}

    //! Total size of the vibrational space at each excitation level, starting from zero (i.e. at 2 = CI0 + CI1 + CI2)
    auto excLvlDim() const {return m_excLvlDim;}

    //! Copy of the vibrational space defining this set of products
    VibSpace vibSpace() const {return m_vibSpace;}

    //! Flags if this product was generated as a connected set
    bool connectedSet() const {return m_connectedSet;}

protected:
    VibSpace m_vibSpace; //!< Definition of the *full* vibrational space of which `this` may be a subset
    std::vector<HProduct> m_basis; //!< vibrational basis as a set of Hartree products
    size_t m_vibDim; //!< Dimension of the vibrational space
    //! Number of vibrational basis functions up to each excitation level, counting from 0 as the GS.
    std::vector<size_t> m_excLvlDim;
    //! Flags `this` is as a subset of the Fock space, connected to another product through an operator.
    //! Some member functions won't work.
    bool m_connectedSet;

    /*!
     * @brief Sets the dimensionality of the vibrational space
     * @warning Cannot be called `if (m_connectedSet)`
     */
    void setVibDim();

    /*!
     * @brief @copybrief HProductSet(const VibSpace&, const HProduct&, const &VibOp)
     * @param vibOp Operator that generates same coupling as Q_A
     * @param bra Reference product
     */
    void generateQcoupledSpace(const HProduct &bra, const VibOp &vibOp);

    /*!
     * @brief @copybrief HProductSet(const VibSpace&, const HProduct&, const &VibOp)
     * @param vibOp Operator that generates same coupling as Q_A*Q_B or Q_A*Q_A
     * @param bra Reference product
     */
    void generateQsqCoupledSpace(const HProduct &bra, const VibOp &vibOp);
};


} //  namespace gci

#endif //GCI_GCIHPRODUCTSET_H
