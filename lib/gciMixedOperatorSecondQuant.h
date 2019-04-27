#ifndef GCI_GCIMIXEDOPERATORSECONDQUANT_H
#define GCI_GCIMIXEDOPERATORSECONDQUANT_H

#include <vector>
#include <SMat.h>

namespace gci {

//using parity_t = SymmetryMatrix::parity_t;
enum class parity_t : int {
    none = 0, //!< no symmetry
    even = 1, //!< symmetric
    odd = -1 //!< antisymmetric
};

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
    explicit MixedTensor(int nMode, int nModal, parity_t hermiticity = parity_t::even,
                         parity_t exchange = parity_t::even) :
            m_nMode(nMode), m_nModal(nModal), m_hermiticity(hermiticity), m_exchange(exchange) { }

    std::pair<parity_t, parity_t> parity() {return {m_hermiticity, m_exchange};}

    /*!
     * @brief Appends a new term to the tensor.
     * @note If a parity equivalent tensor element exists, than nothing is done.
     * @param op Electronic operator
     * @param vibExc Vibrational excitation string
     */
    void append(SymmetryMatrix::SMat &&op, const VibExcitation &vibExc) {
        auto h = hash(vibExc);
        m_tensor.insert({h, op});
        m_excStrings.insert({h, vibExc});
    }

    /*!
     * @brief Unique hash value for this excitation string for a given parity
     * @note Exchange of excitation strings for different modes is alway symmetric.
     * @param parity Symmetry for exchange of modals in excitation string
     * @return
     */
    static size_t hash(const VibExcitation &exc, int nMode, int nModal, parity_t hermiticity, parity_t exchange) {
        if (exc.mc_lvl() == 1) {
            if (hermiticity == parity_t::even || hermiticity == parity_t::odd) return hash_mc1_sym(exc, nMode, nModal);
            else return hash_mc1_nosym(exc, nMode, nModal);
        } else {
            throw std::logic_error("Hash number for MC > 1 not implemented yet");
        }
    }

protected:
    std::map<size_t, SymmetryMatrix::SMat> m_tensor; //! Electronic component of the tensor for each set of vibrational excitations.
    std::map<size_t, VibExcitation> m_excStrings; //! Vibrational exciation strings.
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
    static size_t hash_mc1_sym(const VibExcitation &exc, int nMode, int nModal) {
        size_t A = exc.modes()[0];
        size_t n = std::max(exc.excitations()[0].first, exc.excitations()[0].second);
        size_t m = std::min(exc.excitations()[0].first, exc.excitations()[0].second);
        size_t h_per_mode = nModal * (nModal + 1) / 2;
        return (A - 1) * h_per_mode + n * (n - 1) / 2 + m;
    }

    /*!
     * @brief Hash value for 1 mode operator with no conjugation symmetry
     * Order:
     *      1	4	9	16
     *      2	3	8	15
     *      5	6	7	14
     *      10	11	12	13
     */
    static size_t hash_mc1_nosym(const VibExcitation &exc, int nMode, int nModal) {
        size_t A = exc.modes()[0];
        size_t i = exc.excitations()[0].first;
        size_t j = exc.excitations()[0].second;
        size_t h_per_mode = nModal * nModal;
        if (i >= j) return (A - 1) * h_per_mode + (i - 1) * (i - 1) + j;
        else return (A - 1) * h_per_mode + (j - 1) * (j - 1) + 2 * j - i;
    }
};

/*!
 * @brief Mixed electron-vibration operator in a fully second quantized form.
 *
 * sqH_el = pd q h_pq + 2e- + pd q E_ij^A h_{pq,ij}^A
 *
 * @note Only 1MC operators are currently implemented
 *
 *
 */
class MixedOperatorSecondQuant {
public:
protected:
    SymmetryMatrix::SMat Hel;
    SymmetryMatrix::SMat Hvib;
    MixedTensor mixedHam;
};
} // namespace gci

#endif //GCI_GCIMIXEDOPERATORSECONDQUANT_H
