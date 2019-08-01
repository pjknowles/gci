#include "gciMixedOperatorSecondQuant.h"
#include "gciRun.h"

#include <utility>

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
        includeLambda(fcidump.parameter("INCLUDE_LAMBDA", std::vector<int>{0})[0]),
        includeK(fcidump.parameter("INCLUDE_K", std::vector<int>{0})[0]) {
    std::string f;
    std::string name;
    auto vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                            ns_VibOperator::parity_t::even, name);
    if (includeHel) {
        name = "Hel[1]";
        vibOp.name = name;
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal <= iModal; ++jModal) {
                    f = fcidump.fileName() + "_Hel_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructOperator(FCIdump(f));
                        vibOp.append(op, vibExc);
                        if (iModal != jModal) {
                            vibOp.append(constructOperator(FCIdump(f)), VibExcitation({{iMode, jModal, iModal}}));
                        }
                    }
                }
            }
        }
        mixedHam.insert({name, vibOp});
    }
    if (includeLambda) {
        name = "Lambda[1]";
        vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                           ns_VibOperator::parity_t::even, name);
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal < iModal; ++jModal) {
                    f = fcidump.fileName() + "_Lambda_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructOperatorAntisymm1el(FCIdump(f));
                        vibOp.append(op, vibExc);
                        if (iModal != jModal) {
                            op = constructOperatorAntisymm1el(FCIdump(f)) * (-1.);
                            vibOp.append(op, VibExcitation({{iMode, jModal, iModal}}));
                        }
                    }
                }
            }
        }
        mixedHam.insert({name, vibOp});
    }
    if (includeK) {
        f = fcidump.fileName() + "_K";
        if (file_exists(f)) {
            K0 = std::make_shared<SymmetryMatrix::Operator>(constructK(FCIdump(f)));
        }
        name = "K[1]";
        vibOp = VibOperator<mixed_op_el_t>(nMode, nModal, ns_VibOperator::parity_t::even,
                                           ns_VibOperator::parity_t::even, name);
        for (int iMode = 0; iMode < nMode; ++iMode) {
            for (int iModal = 0; iModal < nModal; ++iModal) {
                for (int jModal = 0; jModal <= iModal; ++jModal) {
                    f = fcidump.fileName() + "_K_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                        "_" + std::to_string(jModal + 1);
                    if (file_exists(f)) {
                        VibExcitation vibExc({{iMode, iModal, jModal}});
                        auto &&op = constructK(FCIdump(f));
                        vibOp.append(op, vibExc);
                        if (iModal != jModal) {
                            op = constructK(FCIdump(f));
                            vibOp.append(op, VibExcitation({{iMode, jModal, iModal}}));
                        }
                    }
                }
            }
        }
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
                                        {-1, -1}, {-1, -1}, 0, true, "Hamiltonian Lambda[1]");

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
            if (si < sj || (si == sj && oi < oj)) {
                std::swap(oi, oj);
                std::swap(si, sj);
            }
            if (sk < sl || (sk == sl && ok < ol)) {
                std::swap(ok, ol);
                std::swap(sk, sl);
            }
            unsigned int sij = si ^sj;

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
    class memory::bytestream bs(buf);
    auto result = SymmetryMatrix::Operator::construct(bs);
    if (rank != 0) free(buf);
    return result;
}

//TODO: Use proper symmetry!
SymmetryMatrix::Operator MixedOperatorSecondQuant::constructK(const FCIdump &dump) {
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
        int n_tot = orbital_symmetries.size();
        for (const auto &s : orbital_symmetries) {
            dim.at(s - 1)++;
        }
//        SymmetryMatrix::Operator result(dim, 2, dump.parameter("IUHF")[0] > 0, 0, true, "Hamiltonian");
        SymmetryMatrix::Operator result(SymmetryMatrix::dims_t{dim, dim, dim, dim}, 1, dump.parameter("IUHF")[0] > 0,
                                        {0, 0}, {0, 0}, 0, true, "Hamiltonian K");

        dump.rewind();
        double value;
        FCIdump::integralType type;
        auto &integrals_a = result.O1(true);
        integrals_a.assign(0);
        auto &integrals_b = result.O1(false);
        integrals_b.assign(0);
        auto &integrals_aa = result.O2(true, true);
        auto &integrals_ab = result.O2(true, false);
        auto &integrals_bb = result.O2(false, false);
        if (verbosity > 0) {
            xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
            xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
            xout << "integral addresses " << &integrals_aa << " " << &integrals_ab << " " << &integrals_bb << std::endl;
            xout << "integral sizes " << integrals_aa.size() << " " << integrals_ab.size() << " " << integrals_bb.size()
                 << std::endl;
        }
        unsigned int si, sj, sk, sl;
        size_t oi, oj, ok, ol;
        while ((type = dump.nextIntegral(si, oi, sj, oj, sk, ok, sl, ol, value)) != FCIdump::endOfFile) {
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
//            if (si < sj || (si == sj && oi < oj)) {
//                std::swap(oi, oj);
//                std::swap(si, sj);
//            }
//            if (sk < sl || (sk == sl && ok < ol)) {
//                std::swap(ok, ol);
//                std::swap(sk, sl);
//            }
            unsigned int sij = si ^sj;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

            if (type == FCIdump::I2aa && false) {
//                (sij ? integrals_aa.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
//                 integrals_aa.smat(sij, si, oi, oj)->block(sk)[ok * (ok + 1) / 2 + ol]) = value;
//                (sij ? integrals_aa.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) :
//                 integrals_aa.smat(sij, sk, ok, ol)->block(si)[oi * (oi + 1) / 2 + oj]) = value;
                integrals_aa.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) = value;
            } else if (type == FCIdump::I2ab && false) {
//                (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
//                 integrals_ab.smat(sij, si, oi, oj)->block(sk)[ok * (ok + 1) / 2 + ol]) = value;
                integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) = value;
            } else if (type == FCIdump::I2bb && false) {
//                (sij ? integrals_bb.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
//                 integrals_bb.smat(sij, si, oi, oj)->block(sk)[ok * (ok + 1) / 2 + ol]) = value;
//                (sij ? integrals_bb.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) :
//                 integrals_bb.smat(sij, sk, ok, ol)->block(si)[oi * (oi + 1) / 2 + oj]) = value;
                integrals_bb.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) = value;
            } else if (type == FCIdump::I1a) {
//                integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
                integrals_a.block(si).at(oi * dim.at(si) + oj) = value;
            } else if (type == FCIdump::I1b) {
//                integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
                integrals_b.block(si).at(oi * dim.at(si) + oj) = value;
            } else if (type == FCIdump::I0)
                result.m_O0 = value;
        }
        if (verbosity > 0) xout << result << std::endl;
        if (verbosity > 1) xout << "int1:\n" << int1(result, 1) << std::endl;
        if (verbosity > 1) xout << "intJ:\n" << intJ(result, 1, 1) << std::endl;
        if (verbosity > 1) xout << "intK:\n" << intK(result, 1) << std::endl;
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
    class memory::bytestream bs(buf);
    auto result = SymmetryMatrix::Operator::construct(bs);
    if (rank != 0) free(buf);
    return result;
}
} // namespace gci
