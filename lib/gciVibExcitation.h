#ifndef GCI_GCIVIBEXCITATION_H
#define GCI_GCIVIBEXCITATION_H

#include <vector>
#include <array>

namespace gci {

/*!
 * @brief Vibrational excitation opeartor string that can operate on a HProduct.
 * @note excitation strings are ordered with increasing mode index
 */
class VibExcitation {
public:
    using exc_t = std::array<int, 3>;
    std::vector<exc_t> excitations; //!< mode indices
    int mc_lvl; //!< number of modes involved in the excitation

    explicit VibExcitation(std::vector<exc_t> excs) : excitations(std::move(excs)), mc_lvl(excitations.size()) { }

    /*!
     * @brief Applies complex conjugation to the excitation operator.
     * Swaps modal labels.
     */
    void conjugate() {
        for (auto &exc : excitations) std::swap(exc[1], exc[2]);
    }

};

} //  namespace gci

#endif //GCI_GCIVIBEXCITATION_H
