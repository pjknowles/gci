#include "gciMixedOperatorSecondQuant.h"
#include "gciRun.h"

namespace gci {

auto file_exists(const std::string &fname) {
    if (std::ifstream{fname}.fail()) {
        std::cout << "Warning (MixedOperatorSecondQuant): fcidump not found   " << fname << std::endl;
        return false;
    }
    return true;
}

MixedOperatorSecondQuant::MixedOperatorSecondQuant(const FCIdump &fcidump) :
        nMode(fcidump.parameter("NMODE", std::vector<int>{0})[0]),
        nModal(fcidump.parameter("NMODAL", std::vector<int>{0})[0]),
        Hel(constructOperator(fcidump)),
        Hvib(constructHvib(fcidump.fileName(), nMode, nModal)),
        includeHel(fcidump.parameter("INCLUDE_HEL", std::vector<int>{0})[0]),
        includeO(fcidump.parameter("INCLUDE_O", std::vector<int>{0})[0]),
        includeK(fcidump.parameter("INCLUDE_K", std::vector<int>{0})[0]) {
    std::string f;
    std::string name;
    auto vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::even,
                                            ns_VibOperator::parity_t::even, name);
    if (includeHel) {
        name = "Hel[1]";
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal <= iModal; ++jModal) {
                    f = fcidump.fileName() + "_Hel_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructOperator(FCIdump(f));
                        vibOp.append(op, vibExc);
                    }
                }
            }
        }
        mixedHam.insert({name, vibOp});
    }
    if (includeO) {
        name = "O[1]";
        vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::odd,
                                           ns_VibOperator::parity_t::even, name);
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal <= iModal; ++jModal) {
                    f = fcidump.fileName() + "_O_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructOperator(FCIdump(f));
                        vibOp.append(op, vibExc);
                    }
                }
            }
        }
        mixedHam.insert({name, vibOp});
    }
    if (includeK) {
        f = fcidump.fileName() + "_K";
        if (file_exists(f)) {
            auto op = constructOperator(FCIdump(f));
            Hel += op;
        }
        vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::even,
                                           ns_VibOperator::parity_t::even, "K[1]");
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal <= iModal; ++jModal) {
                    f = fcidump.fileName() + "_K_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructOperator(FCIdump(f));
                        vibOp.append(op, vibExc);
                    }
                }
            }
        }
        if (mixedHam.find("Hel[1]") != mixedHam.end()) mixedHam.at("Hel[1]") += vibOp;
        else mixedHam.insert({"Hel[1]", vibOp});
    }
}

VibOperator<double> MixedOperatorSecondQuant::constructHvib(const std::string &fcidump_name, int nmode, int nmodal) {
    VibOperator<double> vibOp(nmode, nmodal, ns_VibOperator::parity_t::even, ns_VibOperator::parity_t::even, "Hvib");
    FCIdump dump(fcidump_name + "_Hvib");
    dump.rewind();
    double value;
    int i, j, k, l;
    while (dump.nextIntegral(i, j, k, l, value) != FCIdump::endOfFile) {
        // Assume only 1MC, no constants
        if (i == 0 || j == 0 || k == 0 || l != 0 || i > nmode) throw std::logic_error("Hvib file in the wrong format");
        if (j > nmodal || k > nmodal) continue;
        VibExcitation vibExc({{i - 1, j - 1, k - 1}});
        vibOp.append(value, vibExc);
    }
    return vibOp;
}

SymmetryMatrix::Operator MixedOperatorSecondQuant::constructOperatorAntisymm1el(const FCIdump &dump) {
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
        unsigned int si, sj, sk, sl;
        size_t oi, oj, ok, ol;
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

} // namespace gci
