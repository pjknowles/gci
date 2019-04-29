#ifndef GCI_GCIMIXEDTENSOR_H
#define GCI_GCIMIXEDTENSOR_H

#include <Operator.h>

namespace gci {

/*!
 * @brief Vibrational excitation opeartor string that can operate on a HProduct.
 * @note excitation strings are ordered with increasing mode index
 */
class VibExcitation {
public:
    VibExcitation(std::vector<int> modes, std::vector<std::pair<int, int>> excs) : m_modes(std::move(modes)),
                                                                                   m_excitations(std::move(excs)),
                                                                                   m_mc_lvl(m_modes.size()) { }
    int mc_lvl() const {return m_mc_lvl;}
    const std::vector<int> &modes() const {return m_modes;}
    const std::vector<std::pair<int, int>> &excitations() const {return m_excitations;}
protected:
    std::vector<int> m_modes; //!< mode indices
    std::vector<std::pair<int, int>> m_excitations; //!< modal indices of excitation E_{ij}, first int annihilates, second creates
    int m_mc_lvl; //!< number of modes involved in the excitation
};

/*!
 * @brief Integrals over electronic and vibrational basis, stored as a collection of SMat's
 * with corresponding VibExcitation's.
 */
class MixedTensor {
public:
    using tensor_el_t = std::pair<SymmetryMatrix::Operator, VibExcitation>;
    using tensor_t = std::map<size_t, tensor_el_t>;
    enum class parity_t : int {
        none = 0, //!< no symmetry
        even = 1, //!< symmetric
        odd = -1 //!< antisymmetric
    };

    tensor_t tensor; //!< Mixed tensor
    std::string name; //!< Name of the tensor

    explicit MixedTensor(int nMode, int nModal, parity_t herm = parity_t::even,
                         parity_t exch = parity_t::even, std::string name_ = "") :
            m_nMode(nMode), m_nModal(nModal), m_hermiticity(herm), m_exchange(exch), name(std::move(name_)) { }

    /*!
     * @brief Appends a new term to the tensor.
     * @note If a parity equivalent tensor element exists, than nothing is done.
     * @param op Electronic operator
     * @param vibExc Vibrational excitation string
     */
    void append(SymmetryMatrix::Operator &&op, const VibExcitation &vibExc);

    /*!
     * @brief Access tensor element corresponding to a vibrational excitation string
     * @note Throws std::out_of_range error if element does not exist
     * @param exc Vibrational excitation string
     * @return tensor element
     */
    SymmetryMatrix::Operator &at(const VibExcitation &exc);


    /*!
     * @brief Unique hash value for this excitation string for a given parity
     * @note Exchange of excitation strings for different modes is alway symmetric.
     * @param parity Symmetry for exchange of modals in excitation string
     * @return
     */
    static size_t hash(const VibExcitation &exc, int nMode, int nModal, parity_t hermiticity, parity_t exchange);
protected:
    parity_t m_hermiticity; //!< conjugation symmetry E_{ij} -> E_{ji}
    parity_t m_exchange;//!< symmetry under echange of mode indices E_{ij}^A E_{kl}^B -> E_{kl}^B E_{ij}^A
    int m_nMode; //!< number of modes
    int m_nModal; //!< number of modals per mode (assumed to be the same for each mode)

    size_t hash(const VibExcitation &exc) {return hash(exc, m_nMode, m_nModal, m_hermiticity, m_exchange);}

    /*!
     * @brief Hash value for 1 mode operator with conjugation symmetry
     * Order:
     *      1	0	0	0
     *      2	3	0	0
     *      4	5	6 	0
     *      7	8	9	10
     */
    static size_t hash_mc1_sym(const VibExcitation &exc, int nModal);

    /*!
     * @brief Hash value for 1 mode operator with no conjugation symmetry
     * @note this ordering has no value anymore. I thought I could get around needing nModal, but couldn't.
     * Order:
     *      1	4	9	16
     *      2	3	8	15
     *      5	6	7	14
     *      10	11	12	13
     */
    static size_t hash_mc1_nosym(const VibExcitation &exc, int nModal);

    /*!
     * @brief Hash value for 1 mode operator with no conjugation symmetry
     * @note this ordering has no value anymore. I thought I could get around needing nModal, but couldn't.
     * Order:
     *      1	4	9	16
     *      2	3	8	15
     *      5	6	7	14
     *      10	11	12	13
     */
    static size_t hash_mc1_nosym_old(const VibExcitation &exc, int nModal);
};

} // namespace gci

#endif //GCI_GCIMIXEDTENSOR_H
