#include "gciMixedOperator.h"
#include "gciRun.h"
#include "gciConstants.h"


#include <stdexcept>
#include <valarray>
#include <fstream>

namespace gci {

MixedOperator::MixedOperator
        (const FCIdump &fcidump) : nMode(fcidump.parameter("NMODE", std::vector<int>{0})[0]),
                                   freq(fcidump.parameter("FREQ", std::vector<double>(nMode, 0))),
                                   Hel(constructOperator(fcidump)),
                                   m_inc_d1((bool) fcidump.parameter("INC_D1", std::vector<int>{0})[0]),
                                   m_inc_d2((bool) fcidump.parameter("INC_D2", std::vector<int>{0})[0]),
                                   m_inc_T1((bool) fcidump.parameter("INC_T1", std::vector<int>{0})[0]),
                                   m_inc_T2((bool) fcidump.parameter("INC_T1", std::vector<int>{0})[0]) {
    freq.resize(nMode);
//    std::transform(freq.begin(), freq.end(), freq.begin(), [](auto el) {return el * constants::CM_TO_AU;});
    zpe = 0.5 * std::accumulate(freq.cbegin(), freq.cend(), 0.0);
    std::cout << "MixedOperator: zpe = " << zpe << std::endl;
    auto file_exists = [](const std::string &fname) {
        if (std::ifstream{fname}.fail()) {
            std::cout << "Warning (MixedOperator): fcidump not found --" << fname << std::endl;
            return false;
        }
        return true;
    };
    auto store_fcidump = [&](const std::string &fname, const VibOp &vibOp) {
        if (file_exists(fname)) Hmix[vibOp.type].push_back(MixedOpTerm(vibOp, FCIdump(fname)));
    };
    if (m_inc_d1) {
        for (int iMode = 0; iMode < nMode; ++iMode) {
            std::string f = fcidump.fileName() + "_d1_" + std::to_string(iMode);
            store_fcidump(f, VibOp{VibOpType::Q, {iMode}});
        }
    }
    if (m_inc_d2) {
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int jMode = 0; jMode <= iMode; ++jMode) {
                std::string f = fcidump.fileName() + "_d2_" + std::to_string(iMode) + "_" + std::to_string(jMode);
                store_fcidump(f, VibOp{VibOpType::Qsq, {iMode, jMode}});
            }
        }
        for (auto &op : Hmix[VibOpType::Qsq]) {
            if (op.vibOp.mode[0] != op.vibOp.mode[1]) continue;
            auto i = op.vibOp.mode[0];
            op.Hel.m_O0 -= std::pow(freq[i], 2);
        }
    }
    if (m_inc_T1) {
        for (int iMode = 0; iMode < nMode; ++iMode) {
            std::string f = fcidump.fileName() + "_t1_" + std::to_string(iMode);
            store_fcidump(f, VibOp{VibOpType::dQ, {iMode}});
        }
    }
    if (m_inc_T2) {
        for (int iMode = 0; iMode < nMode; ++iMode) {
            std::string f = fcidump.fileName() + "_t2_" + std::to_string(iMode);
            if (file_exists(f)) {
                auto H_t2 = constructOperator(FCIdump(f));
                Hel += H_t2;
            }
        }
    }
}

double MixedOperator::O_Hvib(const HProduct &bra, const HProduct &ket) const {
//    Suppress error when operator is exactly zero
//    if (bra != ket) throw std::logic_error("HO operator is diagonal. Always 0.");
    if (bra != ket) return 0.0;
    double O_Hvib = zpe;
    for (const auto &el : bra) {
        auto iMode = el[0];
        auto iModal = el[1];
        O_Hvib += freq[iMode] * iModal;
    }
    return O_Hvib;
}

double MixedOperator::expectVal(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const {
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
    std::valarray<int> braOcc(0, nMode);
    std::valarray<int> ketOcc(0, nMode);
    for (const auto &modal : bra) braOcc[modal[0]] = modal[1];
    for (const auto &modal : ket) ketOcc[modal[0]] = modal[1];
    std::valarray<int> exc = std::abs(braOcc - ketOcc);
    int diff = exc.sum();
    double eQsq = 0.0;
    // Q_A * Q_B
    if (vibOp.mode[0] != vibOp.mode[1]) {
        eQsq = 1.0;
        for (auto mode : vibOp.mode) {
//    Suppress error when operator is exactly zero
//            if (exc[mode] != 1) throw std::logic_error("Always 0.");
            auto n = braOcc[mode] > ketOcc[mode] ? braOcc[mode] : ketOcc[mode];
            eQsq *= std::sqrt(n / (2.0 * freq[mode]));
        }
    } else {
        auto mode = vibOp.mode[0];
        // Q_A * Q_A
        if (diff == 0) {
            eQsq = -1.0 / freq[mode] * (braOcc[mode] + 0.5);
        } else if (diff == 2) {
            auto n = braOcc[mode] > ketOcc[mode] ? braOcc[mode] : ketOcc[mode];
            eQsq = 1.0 / (2 * freq[mode]) * std::sqrt(n * (n - 1));
        }
//    Suppress error when operator is exactly zero
//        else throw std::logic_error("Always 0.");
        else eQsq = 0.0;
    }
    eQsq=1;
    return eQsq;
}

double MixedOperator::QtypeOperator(const HProduct &bra, const HProduct &ket,
                                    const std::function<double(double, int, int)> &func,
                                    const int targetMode) const {
    std::valarray<int> braOcc(0, nMode);
    std::valarray<int> ketOcc(0, nMode);
    for (const auto &modal : bra) braOcc[modal[0]] = modal[1];
    for (const auto &modal : ket) ketOcc[modal[0]] = modal[1];
    std::valarray<int> exc = std::abs(braOcc - ketOcc);
    int diff = exc.sum();
//    Suppress error when operator is exactly zero
//    if (diff != 1) throw std::logic_error("Bra and ket are separated by more than 1 excitation. Always 0.");
    if (diff != 1) return 0.0;
    diff = braOcc[targetMode] - ketOcc[targetMode];
    int n = diff > 0 ? braOcc[targetMode] : ketOcc[targetMode];
    return func(freq[targetMode], n, diff);
}

}  // namespace gci
