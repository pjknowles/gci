#include "gciMixedOperator.h"

#include <stdexcept>
#include <valarray>

namespace gci {

double MixedOperator::O_Hvib(const HProduct &bra, const HProduct &ket) const {
    if (bra != ket) throw std::logic_error("HO operator is diagonal. Always 0.");
    double O_Hvib = 0.0;
    for (const auto &el : bra) {
        auto iMode = el[0];
        auto iModal = el[1];
        O_Hvib += 0.5 * freq[iMode] * iModal;
    }
    return O_Hvib;
}

double MixedOperator::expectVal(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) {
    switch (vibOp.type) {
        case VibOpType::HO: return O_Hvib(bra, ket);
        case VibOpType::Q: return O_Q(bra, ket, vibOp);
        case VibOpType::dQ: return O_dQ(bra, ket, vibOp);
        case VibOpType::Qsq: return O_Qsq(bra, ket, vibOp);
    }
}

double MixedOperator::O_Q(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const {
    if (vibOp.type != VibOpType::Q) throw std::logic_error("Wrong operator type");
    auto func = [](double w, const int n, const int diff) {
        return std::sqrt(n / (2.0 * w));
    };
    return QtypeOperator(bra, ket, func, vibOp.mode[0]);
}

double MixedOperator::O_dQ(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const {
    if (vibOp.type != VibOpType::dQ) throw std::logic_error("Wrong operator type");
    auto func = [](double w, int n, int diff) {
        return (double) (diff > 0 ? -1 : 1) * std::sqrt(0.5 * w * n);
    };
    return QtypeOperator(bra, ket, func, vibOp.mode[0]);
}

double MixedOperator::O_Qsq(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const {
    if (vibOp.type != VibOpType::Qsq) throw std::logic_error("Wrong operator type");
    //TODO Just need to finish Qsq and than I can finish MixedWavefunction
}

double MixedOperator::QtypeOperator(const HProduct &bra, const HProduct &ket,
                                    const std::function<double(double, int, int)> &func, const int targetMode) const {
    std::valarray<int> braOcc{0, nMode};
    std::valarray<int> ketOcc{0, nMode};
    for (const auto &modal : bra) braOcc[modal[0]] = modal[1];
    for (const auto &modal : ket) ketOcc[modal[0]] = modal[1];
    int diff = std::abs(braOcc - ketOcc).sum();
    if (diff != 1) throw std::logic_error("Bra and ket are separated by more than 1 excitation. Always 0.");
    diff = braOcc[targetMode] - ketOcc[targetMode];
    int n = diff > 0 ? braOcc[targetMode] : ketOcc[targetMode];
    return func(targetMode, n, diff);
}

}  // namespace gci
