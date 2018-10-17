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
    HProductSet(const VibSpace &vibSpace)
            : m_vibSpace(vibSpace), m_basis(), m_vibDim(0), m_vibExcLvlDim() {setVibDim();};

    /*!
     * @brief Generates a set of products connected to the `bra` through an operator. This is a subset of the full space.
     * @param vibSpace Definition of the full vibrationa lFock space
     * @param bra Reference product
     * @param vibOp Operator connecting reference to products in this set
     */
    HProductSet(const VibSpace &vibSpace, const HProduct &bra, MixedOperator::VibOp vibOp);

    /*!
     * @brief Sets the dimensionality of the vibrational space
     */
    void setVibDim();

    /*!
     * @brief Generates the full vibrational basis
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
     * @param phi Vibrational basis function
     */
    size_t index(const HProduct &phi);

    size_t vibDim() const {return m_vibDim;}
    VibSpace vibSpace() const {return m_vibSpace;}
protected:
    VibSpace m_vibSpace; //!< Definition of the *full* vibrational space of which `this` may be a subset
    std::vector<HProduct> m_basis; //!< vibrational basis as a set of Hartree products
    size_t m_vibDim; //!< Dimension of the vibrational space
    std::vector<size_t> m_vibExcLvlDim; //!< Number of vibrational basis functions up to each excitation level,
    //!< counting from 0 as the GS.
};

HProductSet connectedSpace(const HProduct &bra);


} //  namespace gci

#endif //GCI_GCIHPRODUCTSET_H
