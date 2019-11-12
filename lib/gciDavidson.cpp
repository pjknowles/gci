#include "gciDavidson.h"
#include "gciWavefunction.h"
#include "gciMixedWavefunction.h"
#include "gciUtils.h"

#include <iomanip>
#include <fstream>
#include <ga.h>

namespace gci {
namespace run {


template<class t_Wavefunction, class t_Operator>
Davidson<t_Wavefunction, t_Operator>::Davidson(const t_Wavefunction &prototype, const t_Operator &_ham, Options opt)
        : prototype(prototype),
          ham(_ham),
          options(std::move(opt)),
          energyThreshold(options.parameter("TOL", 1e-13)),
          nState(options.parameter("NSTATE", 1)),
          maxIterations(options.parameter("MAXIT", 1000)),
          solverVerbosity(options.parameter("SOLVER_VERBOSITY", 1)),
          parallel_stringset(options.parameter("PARALLEL_STRINGSET")),
          restart_file(options.parameter("RESTART_FILE", "")),
          backup_file(options.parameter("BACKUP_FILE", "")) {
    solver.m_thresh = energyThreshold;
    solver.m_verbosity = solverVerbosity;
    solver.m_maxIterations = (unsigned int) maxIterations;
    solver.m_roots = (size_t) nState;
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::message() {
    if (GA_Nodeid() == 0) {
        std::cout << "Davidson eigensolver, maximum iterations=" << maxIterations;
        std::cout << "; number of states=" << nState;
        std::cout << "; energy threshold=" << std::scientific << std::setprecision(1) << energyThreshold << std::endl;
        std::cout << "size of Fock space=" << prototype.size() << std::endl;
        xout << std::fixed << std::setprecision(12);
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::printMatrix() {
    if (GA_Nodeid() != 0) return;
    t_Wavefunction w(ww[0]);
    t_Wavefunction action(gg[0]);
    w.allocate_buffer();
    action.allocate_buffer();
    auto n = w.size();
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0));
    if (GA_Nodeid() == 0) {
        std::cout << "Hamiltonian matrix (nxn):" << std::endl;
        std::cout << "n = " << n << std::endl;
        std::cout << " H = {";
    }
    for (size_t i = 0; i < n; ++i) {
        w.zero();
        w.set(i, 1.0);
        action.zero();
        action.operatorOnWavefunction(ham, w);
        if (GA_Nodeid() == 0)
            std::cout << "{";
        for (size_t j = 0; j < n; ++j) {
            H[i][j] = action.at(j);
            if (GA_Nodeid() == 0) {
                if (j == n - 1) std::cout << H[i][j] << "";
                else std::cout << H[i][j] << ",";
            }
        }
        if (i == n - 1) std::cout << "}};" << std::endl;
        else std::cout << "}," << std::endl;
    }
}

template<class t_Wavefunction>
void write_vec(const t_Wavefunction &w, const std::string &message) {
    size_t n = w.size();

    if (GA_Nodeid() == 0)
        std::cout << message << "  {";
    for (size_t i = 0; i < n; ++i) {
        auto v = w.at(i);
        if (GA_Nodeid() == 0)
            std::cout << v << ", ";
    }
    if (GA_Nodeid() == 0)
        std::cout << "}" << std::endl;
}

template<class t_Wavefunction>
typename std::enable_if<!std::is_base_of<Array, t_Wavefunction>::value>::type
davidson_read_write_array(t_Wavefunction &w, const std::string &fname, unsigned int i, hid_t id, bool save) { }

void davidson_read_write_array(Array &w, const std::string &fname, unsigned int i, hid_t id, bool save) {
    auto dataset = utils::open_or_create_hdf5_dataset(id, "result_" + std::to_string(i), H5T_NATIVE_DOUBLE, w.size());
    auto buffer = Array::LocalBuffer(w);
    hsize_t count[1] = {(hsize_t) buffer.size()};
    hsize_t offset[1] = {(hsize_t) buffer.lo};
    auto memspace = H5Screate_simple(1, count, nullptr);
    auto fspace = H5Dget_space(dataset);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    auto plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    if (save) {
        auto error = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, &buffer[0]);
        if (error < 0) throw std::runtime_error("davidson_read_write_array(): write failed");
    } else {
        auto error = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, fspace, plist_id, &buffer[0]);
        if (error < 0) throw std::runtime_error("davidson_read_write_array(): read failed");
    }
    H5Dclose(dataset);
    H5Sclose(fspace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
}

template<typename t_Wavefunction>
void davidson_read_write_wfn(std::vector<t_Wavefunction> &ww, const std::string &fname, bool save) {
    if (fname.empty())
        return;
    auto id = utils::open_hdf5_file(fname, mpi_comm_compute, save);
    for (auto i = 0ul; i < ww.size(); ++i) {
        davidson_read_write_array(ww[i], fname, i, id, save);
    }
    H5Fclose(id);
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::prepareGuess() { }

template<>
void Davidson<MixedWavefunction, MixedOperatorSecondQuant>::prepareGuess() {
    if (!restart_file.empty()) {
        davidson_read_write_wfn<MixedWavefunction>(ww, restart_file, false);
        return;
    }
    // Currently assumes only 1 mode
    auto prof = profiler->push("prepareGuess");
    auto nMode = options.parameter("NMODE", int(0));
    if (nMode != 1) return;
    for (auto &root : ww) root.zero();
    auto id = GA_Nodeid();
    if (id == 0) {
        std::cout << "Entered Davidson::prepareGuess()" << std::endl;
        Wavefunction w = prototype.wavefunctionAt(0, ww[0].m_child_communicator);
        // get number of electronic states from options
        auto nM = options.parameter("NMODAL", int(0));
        int nElSt = nState / nM + nState % nM ? 1 : 0;
        nElSt = options.parameter("NELSTATE_INIT", nElSt); // number of electronic states to search for
        auto modOptions = Options(options);
        modOptions.addParameter("NSTATE", nElSt);
        modOptions.addParameter("MAXIT", (int) 100);
        modOptions.addParameter("TOL", options.parameter("TOLGUESS", energyThreshold));
        modOptions.addParameter("BACKUP_FILE", "");
        modOptions.addParameter("RESTART_FILE", "");
        auto h = ham.elHam.at("Hel[0]").get();
        Davidson<Wavefunction, SymmetryMatrix::Operator> elecSolver(w, *h, modOptions);
        elecSolver.run();
        // get initial guess vectors from options
        auto nElStG = options.parameter("NELSTATE_GUESS",
                                        nElSt); // number of lowest energy electronic states to construct the guess from
        int guess_size = nState * nM * nElStG;
        std::vector<double> init_guess;
        auto guess_file = options.parameter("INIT_GUESS", std::string(""));
        if (guess_file.empty()) {
            init_guess.resize(guess_size, 0.);
            for (auto i = 0ul; i < nState; ++i) {
                auto iVibSt = i % nM;
                auto iElSt = i / nM;
                init_guess[i * nM * nElStG + iVibSt * nElStG + iElSt] = 1.0;
            }
        } else {
            xout << "Read guess from file, " << guess_file << std::endl;
            auto input = std::ifstream(guess_file);
            std::string val;
            while (std::getline(input, val))
                if (!val.empty()) {
                    init_guess.push_back(std::stod(val));
                }
            if (init_guess.size() != guess_size) GA_Error((char *) "initial guess vector is of the wrong size", 0);
        }
        for (auto i = 0ul; i < nState; ++i) {
            for (int iVibSt = 0; iVibSt < nM; ++iVibSt) {
                auto ind = i * nM * nElStG + iVibSt * nElStG;
                auto weight = init_guess.at(ind);
                auto wfn = weight * elecSolver.ww[0];
                for (auto j = 1ul; j < nElStG; ++j) {
                    weight = init_guess.at(ind + j);
                    wfn += weight * elecSolver.ww[0];
                }
                ww[i].put(iVibSt, wfn);
            }
        }
        std::cout << "Exit Davidson::prepareGuess()" << std::endl;
    }
    for (auto &root : ww) root.sync();
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::run() {
    auto prof = profiler->push("Davidson");
    message();
    initialize();
    prepareGuess();
    backup(ww);
//    printMatrix();
    for (unsigned int iteration = 1; iteration <= maxIterations; iteration++) {
        action();
        solver.addVector(ww, gg);
        backup(ww);
        update();
        if (solver.endIteration(ww, gg)) break;
    }
    if (maxIterations > 0) {
        xout << "energies: ";
        for (unsigned int i = 0; i < nState; ++i) xout << solver.eigenvalues()[i] << ", ";
        xout << std::endl;
    }
    std::cout << "Exit Davidson::run()" << std::endl;
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::initialize() {
    if (!diagonalH) diagonalH = std::make_shared<t_Wavefunction>(prototype, 0);
    diagonalH->allocate_buffer();
    diagonalH->diagonalOperator(ham);
    diag_minlocN = diagonalH->minlocN(nState);
    for (auto i = 0ul; i < nState; ++i)
        diag_val_at_minlocN.push_back(diagonalH->at(diag_minlocN[i]));
    std::vector<int> roots(nState, 0);
    for (unsigned int root = 0; root < nState; root++) {
        ww.emplace_back(prototype, 0);
        ww.back().allocate_buffer();
        ww.back().zero();
        auto n = diag_minlocN[root];
        if (GA_Nodeid() == 0)
            if (std::count(roots.begin(), roots.begin() + root, n) != 0)
                throw std::logic_error("Davidson::initialize duplicate guess vector, n =" + std::to_string(n));
        roots[root] = n;
        ww.back().set(n, 1.0);
        gg.emplace_back(prototype, 0);
        gg.back().allocate_buffer();
        gg.back().settilesize(
                options.parameter("TILESIZE", std::vector<int>(1, -1)).at(0),
                options.parameter("ALPHATILESIZE", std::vector<int>(1, -1)).at(0),
                options.parameter("BETATILESIZE", std::vector<int>(1, -1)).at(0));
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::action() {
    for (size_t k = 0; k < ww.size(); k++) {
        const t_Wavefunction &x = ww[k];
        t_Wavefunction &g = gg[k];
        g.zero();
        if (solver.active().at(k)) {
            auto prof = profiler->push("Hc");
            g.operatorOnWavefunction(ham, x, false);
        }
    }
}

template<>
void Davidson<MixedWavefunction, MixedOperatorSecondQuant>::action() {
    for (auto &g: gg) g.zero();
    DivideTasks(10000000000, 1, 1, gg[0].m_communicator);
    for (size_t k = 0; k < ww.size(); k++) {
        auto &x = ww[k];
        auto &g = gg[k];
        if (solver.active().at(k)) {
            auto prof = profiler->push("Hc");
            g.operatorOnWavefunction(ham, x, false, false);
        }
    }
    for (const auto &g: gg) g.sync();
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::update() {
    auto eigval = solver.eigenvalues();
    for (size_t state = 0; state < nState; state++) {
        t_Wavefunction &cw = ww[state];
        const t_Wavefunction &gw = gg[state];
        auto shift = -eigval[state] + 1e-10;
        shift += 2 * std::numeric_limits<value_type>::epsilon()
                 * std::max<value_type>(1, std::abs(diag_val_at_minlocN[state])); // to guard against zero
        cw.divide(&gw, diagonalH.get(), shift, true, true);
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::backup(std::vector<t_Wavefunction> &ww) {
    davidson_read_write_wfn<t_Wavefunction>(ww, backup_file, true);
}

template
class Davidson<gci::MixedWavefunction, gci::MixedOperatorSecondQuant>;

template
class Davidson<gci::Wavefunction, SymmetryMatrix::Operator>;

template
void davidson_read_write_wfn<gci::Array>(std::vector<Array> &ww, const std::string &fname, bool save);
}  // namespace run
}  // namespace gci
