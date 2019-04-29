#ifndef GCI_GCIMIXEDOPERATORSECONDQUANT_H
#define GCI_GCIMIXEDOPERATORSECONDQUANT_H

#include "gciVibOperator.h"

#include <FCIdump.h>

namespace gci {

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
    SymmetryMatrix::Operator Hel; //!< Purely electronic term
    std::shared_ptr<VibOperator<double>> Hvib;//!< Purely vibrational term
    std::vector<VibOperator<SymmetryMatrix::Operator >> mixedHam;//!< Mixed electronic-vibrational terms
    int nMode; //!< Number of vibrational modes
    int nModal; //!< Number of modals per mode (for now assumed the same for each mode)

    explicit MixedOperatorSecondQuant(const FCIdump &fcidump);
    MixedOperatorSecondQuant() : Hel({0, 0}, 0, false, 0, true, "dummy"), mixedHam({}), nMode(0), nModal(0) { }
};
} // namespace gci

#endif //GCI_GCIMIXEDOPERATORSECONDQUANT_H
