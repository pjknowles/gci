#include "gciMixedOperator.h"

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

double MixedOperator::O_Q(const HProduct &bra, const HProduct &ket) const {
    auto func = [&](const int mode, const int n, const int diff) {return std::sqrt(n / (2.0 * freq[mode]));};
    return QtypeOperator(bra, ket, func);
}

double MixedOperator::O_dQ(const HProduct &bra, const HProduct &ket) const {
    auto sign = [](const int n) {return n > 0 ? 1 : -1;};
    auto func = [&](const int mode, const int n, const int diff) {
        return sign(diff) * std::sqrt(0.5 * freq[mode] * n);
    };
    return QtypeOperator(bra, ket, func);
}

double MixedOperator::expectVal(const HProduct &bra, const HProduct &ket, const VibOp vibOp) {
    switch (vibOp) {
        case VibOp::HO: return O_Hvib(bra, ket);
        case VibOp::Q: return O_Q(bra, ket);
        case VibOp::dQ: return O_dQ(bra, ket);
        case VibOp::Qsq: return O_Qsq(bra, ket);
    }
}

}  // namespace gci
