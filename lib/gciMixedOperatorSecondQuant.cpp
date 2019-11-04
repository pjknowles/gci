#include "gciMixedOperatorSecondQuant.h"
#include "gciRun.h"
#include "gciUtils.h"
#include "gciPersistentOperator.h"

#include <utility>

namespace gci {

inline auto _fcidump_f(const Options &options) {return options.parameter("FCIDUMP", "");}

inline auto _nMode(const Options &options) {return options.parameter("NMODE", 0);}

inline auto _nModal(const Options &options) {return options.parameter("NMODAL", 0);}

inline auto _hdf5_fname(const Options &options) {
    auto fname_save = options.parameter("HAM_HDF5", "");
    auto fname_restart = options.parameter("HAM_HDF5_RESTART", "");
    if (fname_save.empty() && fname_restart.empty())
        throw std::runtime_error("MixedOperatorSecondQuant::MixedOperatorSecondQuant() hdf5_fname is empty");
    return fname_restart.empty() ? fname_save : fname_restart;
}

MixedOperatorSecondQuant::MixedOperatorSecondQuant(const Options &options) :
        m_fcidump_f(_fcidump_f(options)),
        hdf5_fname(_hdf5_fname(options)),
        restart(!options.parameter("HAM_HDF5_RESTART", "").empty()),
        hid_file(utils::open_hdf5_file(hdf5_fname, gci::mpi_comm_compute, !restart)),
        includeHel(options.parameter("INCLUDE_HEL", 0)),
        includeLambda(options.parameter("INCLUDE_LAMBDA", 0)),
        includeK(options.parameter("INCLUDE_K", 0)),
        includeD(options.parameter("INCLUDE_D", 0)),
        nMode(_nMode(options)),
        nModal(_nModal(options)),
        Hvib(constructHvib(_fcidump_f(options), _nMode(options), _nModal(options))) {
    FCIdump fcidump(m_fcidump_f);
    if (includeHel) initializeHel(fcidump);
    if (includeLambda) initializeLambda(fcidump);
    if (includeK) initializeK(fcidump);
    if (includeD) initializeD(fcidump);
}

MixedOperatorSecondQuant::~MixedOperatorSecondQuant() {H5Fclose(hid_file);}

template<class Op>
inline PersistentOperator create_persistentoperator(const std::string &fcidump, bool restart, std::string description,
                                                    hid_t hid_file, int root,
                                                    Op(*construct_op)(const FCIdump &, bool)) {
    if (restart) return PersistentOperator(hid_file, description);
    std::shared_ptr<Op> op;
    if (gci::parallel_rank == root) {
        op = std::make_shared<Op>(construct_op(FCIdump(fcidump), false));
        op->m_description = description;
        op->ensure_dirac();
    }
    auto p_op = PersistentOperator(op, description, root, hid_file);
    return p_op;
}

void MixedOperatorSecondQuant::initializeHel(const FCIdump &fcidump) {
    int i_operator = 0;
    std::string f = fcidump.fileName();
    if (utils::file_exists(f, "MixedOperatorSecondQuant::initializeHel(): fcidump not found, " + f)) {
        auto description = std::string("Hel[0]");
        int root = i_operator % gci::parallel_size;
        ++i_operator;
        auto p_op = create_persistentoperator(f, restart, description, hid_file, root, constructOperator);
        elHam.insert({description, std::move(p_op)});
    }
    std::string name = "Hel[1]";
    auto vibOp = VibOperator<hel_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                    ns_VibOperator::parity_t::even, name);
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 0; iModal < nModal; ++iModal) {
            for (int jModal = 0; jModal <= iModal; ++jModal) {
                f = fcidump.fileName() + "_Hel_" + std::to_string(iMode + 1) + "_" +
                    std::to_string(iModal + 1) + "_" + std::to_string(jModal + 1);
                if (!utils::file_exists(f, "MixedOperatorSecondQuant::initializeHel(): fcidump not found, " + f))
                    continue;
                VibExcitation vibExc({{iMode, iModal, jModal}});
                std::string description = name + "(" + std::to_string(iMode + 1) + ","
                                          + std::to_string(iModal + 1) + ","
                                          + std::to_string(jModal + 1) + ")";
                int root = i_operator % gci::parallel_size;
                ++i_operator;
                auto p_op = create_persistentoperator(f, restart, description, hid_file, root, constructOperator);
                vibOp.append(p_op, vibExc);
                if (iModal != jModal)
                    vibOp.append(p_op, VibExcitation({{iMode, jModal, iModal}}));
            }
        }
    }
    mixedHam.insert({name, vibOp});
}

void MixedOperatorSecondQuant::initializeLambda(const FCIdump &fcidump) {
    std::string f, name = "Lambda[1]";
    auto vibOp = VibOperator<hel_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                    ns_VibOperator::parity_t::even, name);
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 0; iModal < nModal; ++iModal) {
            for (int jModal = 0; jModal < iModal; ++jModal) {
                f = fcidump.fileName() + "_Lambda_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                    "_" + std::to_string(jModal + 1);
                if (!utils::file_exists(f, "MixedOperatorSecondQuant::initializeLambda(): fcidump not found, " + f))
                    continue;
                VibExcitation vibExc({{iMode, iModal, jModal}});
                std::string description = name + "(" + std::to_string(iMode + 1) + ","
                                          + std::to_string(iModal + 1) + ","
                                          + std::to_string(jModal + 1) + ")";
                auto op = std::make_shared<SymmetryMatrix::Operator>(
                        constructOperatorAntisymm1el(FCIdump(f), true));
                op->m_description = description;
                op->ensure_dirac();
                auto p_op = PersistentOperator(op, description, 0, hid_file, true);
                vibOp.append(p_op, vibExc);
                if (iModal != jModal) {
                    op = std::make_shared<SymmetryMatrix::Operator>(
                            constructOperatorAntisymm1el(FCIdump(f), true) * (-1.));
                    op->ensure_dirac();
                    p_op = PersistentOperator(op, description, 0, hid_file, true);
                    vibOp.append(p_op, VibExcitation({{iMode, jModal, iModal}}));
                }

            }
        }
    }
    mixedHam.insert({name, vibOp});
}

void MixedOperatorSecondQuant::initializeK(const FCIdump &fcidump) {
    std::string f = fcidump.fileName() + "_K";
    if (utils::file_exists(f, "MixedOperatorSecondQuant::initializeK(): fcidump not found, " + f)) {
        auto description = std::string("K[0]");
        auto op = std::make_shared<SymmetryMatrix::Operator>(constructK(FCIdump(f), true));
        op->m_description = description;
        op->ensure_dirac();
        auto p_op = PersistentOperator(op, description, 0, hid_file, true);
        elHam.insert({"K[0]", p_op});
    }
    std::string name = "K[1]";
    auto vibOp = VibOperator<hel_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                    ns_VibOperator::parity_t::even, name);
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 0; iModal < nModal; ++iModal) {
            for (int jModal = 0; jModal <= iModal; ++jModal) {
                f = fcidump.fileName() + "_K_" + std::to_string(iMode + 1) + "_" + std::to_string(iModal + 1) +
                    "_" + std::to_string(jModal + 1);
                if (!utils::file_exists(f, "MixedOperatorSecondQuant::initializeK(): fcidump not found, " + f))
                    continue;
                VibExcitation vibExc({{iMode, iModal, jModal}});
                std::string description = name + "(" + std::to_string(iMode + 1) + ","
                                          + std::to_string(iModal + 1) + ","
                                          + std::to_string(jModal + 1) + ")";
                auto op = std::make_shared<SymmetryMatrix::Operator>(constructK(FCIdump(f), true));
                op->m_description = description;
                op->ensure_dirac();
                auto p_op = PersistentOperator(op, description, 0, hid_file, true);
                vibOp.append(p_op, vibExc);
                if (iModal != jModal)
                    vibOp.append(p_op, VibExcitation({{iMode, jModal, iModal}}));
            }
        }
    }
    mixedHam.insert({name, vibOp});
}

void MixedOperatorSecondQuant::initializeD(const FCIdump &fcidump) {
    int i_operator = 0;
    std::string f = fcidump.fileName() + "_D";
    if (utils::file_exists(f, "MixedOperatorSecondQuant::initializeD(): fcidump not found, " + f)) {
        auto description = std::string("D[0]");
        int root = i_operator % gci::parallel_size;
        ++i_operator;
        auto p_op = create_persistentoperator(f, restart, description, hid_file, root, constructD);
        elHam.insert({description, std::move(p_op)});
    }
    std::string name = "D[1]";
    auto vibOp = VibOperator<hel_t>(nMode, nModal, ns_VibOperator::parity_t::none,
                                    ns_VibOperator::parity_t::even, name);
    for (int iMode = 0; iMode < nMode; ++iMode) {
        for (int iModal = 0; iModal < nModal; ++iModal) {
            for (int jModal = 0; jModal <= iModal; ++jModal) {
                f = fcidump.fileName() + "_D_" + std::to_string(iMode + 1) + "_" +
                    std::to_string(iModal + 1) + "_" + std::to_string(jModal + 1);
                if (!utils::file_exists(f, "MixedOperatorSecondQuant::initializeD(): fcidump not found, " + f))
                    continue;
                VibExcitation vibExc({{iMode, iModal, jModal}});
                std::string description = name + "(" + std::to_string(iMode + 1) + ","
                                          + std::to_string(iModal + 1) + ","
                                          + std::to_string(jModal + 1) + ")";
                int root = i_operator % gci::parallel_size;
                ++i_operator;
                auto p_op = create_persistentoperator(f, restart, description, hid_file, root, constructD);
                vibOp.append(p_op, vibExc);
                if (iModal != jModal)
                    vibOp.append(p_op, VibExcitation({{iMode, jModal, iModal}}));

            }
        }
    }
    mixedHam.insert({name, vibOp});
}

VibOperator<double> MixedOperatorSecondQuant::constructHvib(const std::string &fcidump_name, int nmode, int nmodal) {
    VibOperator<double> vibOp(nmode, nmodal, ns_VibOperator::parity_t::none, ns_VibOperator::parity_t::even, "Hvib");
    FCIdump dump(fcidump_name + "_Hvib");
    dump.rewind();
    double value;
    int iMode, jModal, kModal, l;
    while (dump.nextIntegral(iMode, jModal, kModal, l, value) != FCIdump::endOfFile) {
        // Assume only 1MC, no constants
        if (iMode == 0 || jModal == 0 || kModal == 0 || l != 0 || iMode > nmode)
            throw std::logic_error("Hvib file in the wrong format");
        if (jModal > nmodal || kModal > nmodal) continue;
        VibExcitation vibExc({{iMode - 1, jModal - 1, kModal - 1}});
        vibOp.append(value, vibExc);
        if (jModal != kModal) {
            vibExc.conjugate();
            vibOp.append(value, vibExc);
        }
    }
    return vibOp;
}

SymmetryMatrix::Operator MixedOperatorSecondQuant::constructOperatorAntisymm1el(const FCIdump &dump, bool collective) {
    std::vector<char> portableByteStream;
    int lPortableByteStream;
    auto rank = gci::parallel_rank;
    if (rank == 0 || !collective) {
        int verbosity = 0;
        std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
        SymmetryMatrix::dim_t dim(8);
        for (unsigned int s : orbital_symmetries) {
            dim.at(s - 1)++;
        }
        SymmetryMatrix::Operator result(SymmetryMatrix::dims_t{dim, dim, dim, dim}, 1, dump.parameter("IUHF")[0] > 0,
                                        {-1, -1}, {-1, -1}, 0, true, false, "Hamiltonian Lambda[1]");
        result.zero();

        dump.rewind();
        double value;
        FCIdump::integralType type;
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
            int phase = 1;
            if (si < sj || (si == sj && oi < oj)) {
                std::swap(oi, oj);
                std::swap(si, sj);
                phase *= -1;
            }
            if (sk < sl || (sk == sl && ok < ol)) {
                std::swap(ok, ol);
                std::swap(sk, sl);
                phase *= -1;
            }
            value *= phase;

            if (type == FCIdump::I1a) {
                if (verbosity > 1) xout << "ha(" << oi << "," << oj << ") = " << value << std::endl;
                if (si != sj) continue;
                integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I1b) {
                if (verbosity > 1) xout << "hb(" << oi << "," << oj << ") = " << value << std::endl;
                if (si != sj) continue;
                integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I0)
                result.m_O0 = value;
        }
        if (verbosity > 0) xout << result << std::endl;
        if (collective) {
            portableByteStream = result.bytestream().data();
            lPortableByteStream = portableByteStream.size();
        } else
            return result;
    }
    if (collective) {
#ifdef HAVE_MPI_H
        MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, mpi_comm_compute);
#endif
        char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
        MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, mpi_comm_compute);
#endif
        class memory::bytestream bs(buf);
        auto result = SymmetryMatrix::Operator::construct(bs);
        if (rank != 0) free(buf);
        return result;
    }
}

SymmetryMatrix::Operator MixedOperatorSecondQuant::constructK(const FCIdump &dump, bool collective) {
    std::vector<char> portableByteStream;
    int lPortableByteStream;
    int rank = gci::parallel_rank;
    if (rank == 0 || !collective) {
        int verbosity = 0;
        std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
        SymmetryMatrix::dim_t dim(8);
        for (unsigned int s : orbital_symmetries) {
            dim.at(s - 1)++;
        }
//        SymmetryMatrix::Operator result(dim, 2, dump.parameter("IUHF")[0] > 0, 0, true, "Hamiltonian");
        SymmetryMatrix::Operator result(SymmetryMatrix::dims_t{dim, dim, dim, dim}, 1, dump.parameter("IUHF")[0] > 0,
                                        {1, 1}, {-1, -1}, 0, true, false, "Hamiltonian K");
        result.zero();

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
            if (si < sj || (si == sj && oi < oj)) {
                std::swap(oi, oj);
                std::swap(si, sj);
            }
            if (sk < sl || (sk == sl && ok < ol)) {
                std::swap(ok, ol);
                std::swap(sk, sl);
            }
            if (type == FCIdump::I1a) {
                if (si != sj) continue;
                integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I1b) {
                if (si != sj) continue;
                integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
            } else if (type == FCIdump::I0)
                result.m_O0 = value;
        }
        if (verbosity > 0) xout << result << std::endl;
        if (verbosity > 1) xout << "int1:\n" << int1(result, 1) << std::endl;
        if (verbosity > 1) xout << "intJ:\n" << intJ(result, 1, 1) << std::endl;
        if (verbosity > 1) xout << "intK:\n" << intK(result, 1) << std::endl;
        if (collective) {
            portableByteStream = result.bytestream().data();
            lPortableByteStream = portableByteStream.size();
        } else
            return result;
    }
    if (collective) {
#ifdef HAVE_MPI_H
        MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, mpi_comm_compute);
#endif
        char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
        MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, mpi_comm_compute);
#endif
        class memory::bytestream bs(buf);
        auto result = SymmetryMatrix::Operator::construct(bs);
        if (rank != 0) free(buf);
        return result;
    }
}

SymmetryMatrix::Operator MixedOperatorSecondQuant::constructD(const FCIdump &dump, bool collective) {
    std::vector<char> portableByteStream;
    int lPortableByteStream;
    int rank = gci::parallel_rank;
    if (rank == 0 || !collective) {
        int verbosity = 0;
        std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
        SymmetryMatrix::dim_t dim(8);
        for (unsigned int s : orbital_symmetries) {
            dim.at(s - 1)++;
        }
//        SymmetryMatrix::Operator result(dim, 2, dump.parameter("IUHF")[0] > 0, 0, true, "Hamiltonian");
        SymmetryMatrix::Operator result(SymmetryMatrix::dims_t{dim, dim, dim, dim}, 2, dump.parameter("IUHF")[0] > 0,
                                        {-1, -1}, {-1, -1}, 0, true, false, "Hamiltonian D");
        result.zero();

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
            int phase = 1;
            if (si < sj || (si == sj && oi < oj)) {
                std::swap(oi, oj);
                std::swap(si, sj);
                phase *= -1;
            }
            if (sk < sl || (sk == sl && ok < ol)) {
                std::swap(ok, ol);
                std::swap(sk, sl);
                phase *= -1;
            }
            unsigned int sij = si ^sj;
            unsigned int skl = sk ^sl;
            value *= phase;

            if (type == FCIdump::I2aa) {
                if (sij != skl) continue;
                (sij ? integrals_aa.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
                 integrals_aa.smat(sij, si, oi, oj)->block(sk)[ok * (ok - 1) / 2 + ol]) = value;
                (sij ? integrals_aa.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) :
                 integrals_aa.smat(sij, sk, ok, ol)->block(si)[oi * (oi - 1) / 2 + oj]) = value;
            } else if (type == FCIdump::I2ab) {
                if (sij != skl) continue;
                (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
                 integrals_ab.smat(sij, si, oi, oj)->block(sk)[ok * (ok - 1) / 2 + ol]) = value;
            } else if (type == FCIdump::I2bb) {
                if (sij != skl) continue;
                (sij ? integrals_bb.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) :
                 integrals_bb.smat(sij, si, oi, oj)->block(sk)[ok * (ok - 1) / 2 + ol]) = value;
                (sij ? integrals_bb.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) :
                 integrals_bb.smat(sij, sk, ok, ol)->block(si)[oi * (oi - 1) / 2 + oj]) = value;
            } else if (type == FCIdump::I1a) {
                if (si != sj) continue;
                integrals_a.block(si).at(oi * (oi - 1) / 2 + oj) = value;
            } else if (type == FCIdump::I1b) {
                if (si != sj) continue;
                integrals_b.block(si).at(oi * (oi - 1) / 2 + oj) = value;
            } else if (type == FCIdump::I0)
                result.m_O0 = value;
        }
        if (verbosity > 0) xout << result << std::endl;
        if (verbosity > 1) xout << "int1:\n" << int1(result, 1) << std::endl;
        if (verbosity > 1) xout << "intJ:\n" << intJ(result, 1, 1) << std::endl;
        if (verbosity > 1) xout << "intK:\n" << intK(result, 1) << std::endl;
        if (collective) {
            portableByteStream = result.bytestream().data();
            lPortableByteStream = portableByteStream.size();
        } else
            return result;
    }
    if (collective) {
#ifdef HAVE_MPI_H
        MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, mpi_comm_compute);
#endif
        char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
        MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, mpi_comm_compute);
#endif
        class memory::bytestream bs(buf);
        auto result = SymmetryMatrix::Operator::construct(bs);
        if (rank != 0) free(buf);
        return result;
    }
}

bool MixedOperatorSecondQuant::connected(const HProduct &bra, const HProduct &ket) const {
    for (const auto &mixedTerm : mixedHam) {
        const auto &vibTensor = mixedTerm.second;
        for (const auto &vibEl : vibTensor.tensor) {
            auto &op = vibEl.second.oper;
            auto vibExc = VibExcitation{vibEl.second.exc};
            auto new_bra = ket.excite(vibExc);
            if (new_bra == bra) return true;
        }
    }
    return false;
}


} // namespace gci
