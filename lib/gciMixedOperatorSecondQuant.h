#ifndef GCI_GCIMIXEDOPERATORSECONDQUANT_H
#define GCI_GCIMIXEDOPERATORSECONDQUANT_H

#include "gciMixedTensor.h"

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
protected:
    SymmetryMatrix::SMat Hel;
    SymmetryMatrix::SMat Hvib;
    std::vector<MixedTensor> mixedHam;
};
} // namespace gci

#endif //GCI_GCIMIXEDOPERATORSECONDQUANT_H
