#include "gciMixedOperator.h"
#include "gciRun.h"
#include "gciConstants.h"


#include <stdexcept>
#include <valarray>
#include <fstream>

namespace gci {
SymmetryMatrix::Operator constructOperatorT1(const FCIdump &dump) {
    std::vector<char> portableByteStream;
    int lPortableByteStream;
    int rank = 0;
#ifdef HAVE_MPI_H
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
#endif
    if (rank == 0) {
        int verbosity = 0;
        std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
        SymmetryMatrix::dim_t dim(8);
        for (const auto &s : orbital_symmetries) {
            dim.at(s - 1)++;
        }
        SymmetryMatrix::Operator result(SymmetryMatrix::dims_t{dim, dim, dim, dim}, 1, dump.parameter("IUHF")[0] > 0,
                                        {-1, -1}, {-1, -1}, 0, true, "Hamiltonian T1");

        dump.rewind();
        double value;
        FCIdump::integralType type;
        int i, j, k, l;
        auto &integrals_a = result.O1(true);
        integrals_a.assign(0);
        auto &integrals_b = result.O1(false);
        integrals_b.assign(0);
        if (verbosity > 0) {
            xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
            xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
        }
        unsigned int si,sj,sk,sl;
        size_t oi,oj,ok,ol;
        while ((type = dump.nextIntegral(si, oi, sj, oj, sk, ok, sl, ol, value)) != FCIdump::endOfFile) {
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
            if (si < sj || (si == sj && oi < oj)) {
                std::swap(oi, oj);
                std::swap(si, sj);
            }
            if (sk < sl || (sk == sl && ok < ol)) {
                std::swap(ok, ol);
                std::swap(sk, sl);
            }
            unsigned int sij = si ^sj;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

            if (type == FCIdump::I1a) {
                if (verbosity > 1) xout << "ha(" << i << "," << j << ") = " << value << std::endl;
                integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I1b) {
                if (verbosity > 1) xout << "hb(" << i << "," << j << ") = " << value << std::endl;
                integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I0)
                result.m_O0 = value;
        }
        if (verbosity > 0) xout << result << std::endl;
        portableByteStream = result.bytestream().data();
        lPortableByteStream = portableByteStream.size();
    }
#ifdef HAVE_MPI_H
    MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, MPI_COMM_COMPUTE);
#endif
    char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
    MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, MPI_COMM_COMPUTE);
#endif
    class bytestream bs(buf);
    auto result = SymmetryMatrix::Operator::construct(bs);
    if (rank != 0) free(buf);
    return result;
}

MixedOperator::MixedOperator
        (const FCIdump &fcidump) : nMode(fcidump.parameter("NMODE", std::vector<int>{0})[0]),
                                   freq(fcidump.parameter("FREQ", std::vector<double>(nMode, 0))),
                                   Hel(constructOperator(fcidump)),
                                   m_inc_d1((bool) fcidump.parameter("INC_D1", std::vector<int>{0})[0]),
                                   m_inc_d2((bool) fcidump.parameter("INC_D2", std::vector<int>{0})[0]),
                                   m_inc_T1((bool) fcidump.parameter("INC_T1", std::vector<int>{0})[0]),
                                   m_inc_T2((bool) fcidump.parameter("INC_T2", std::vector<int>{0})[0]) {
    freq.resize(nMode);
//    std::transform(freq.begin(), freq.end(), freq.begin(), [](auto el) {return el * constants::CM_TO_AU;});
    zpe = 0.5 * std::accumulate(freq.cbegin(), freq.cend(), 0.0);
    std::cout << "MixedOperator: zpe = " << zpe << std::endl;
    auto file_exists = [](const std::string &fname) {
        if (std::ifstream{fname}.fail()) {
            std::cout << "Warning (MixedOperator): fcidump not found   " << fname << std::endl;
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
            if (file_exists(f))
                Hmix[VibOpType::dQ].push_back(
                        MixedOpTerm(VibOp{VibOpType::dQ, {iMode}}, constructOperatorT1(FCIdump(f))));
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
            eQsq = 0.5 * 1.0 / freq[mode] * (braOcc[mode] + 0.5);
        } else if (diff == 2) {
            auto n = braOcc[mode] > ketOcc[mode] ? braOcc[mode] : ketOcc[mode];
            eQsq = 0.5 * 1.0 / (2. * freq[mode]) * std::sqrt(n * (n - 1));
        }
//    Suppress error when operator is exactly zero
//        else throw std::logic_error("Always 0.");
        else eQsq = 0.0;
    }
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
