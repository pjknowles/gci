#include "gciDavidson.h"

#include <iomanip>
#include <ga.h>

namespace gci {
namespace run {


template<class t_Wavefunction, class t_Operator>
Davidson<t_Wavefunction, t_Operator>::Davidson(t_Wavefunction &&prototype,
                                               t_Operator &&ham, Options &opt)
        : prototype(std::make_shared<t_Wavefunction>(std::move(prototype))),
          ham(std::make_shared<t_Operator>(std::move(ham))),
          options(opt),
          energyThreshold(options.parameter("TOL", 1e-13)),
          nState(options.parameter("NSTATE", 1)),
          maxIterations(options.parameter("MAXIT", 1000)),
          solverVerbosity(options.parameter("SOLVER_VERBOSITY", 1)),
          parallel_stringset(options.parameter("PARALLEL_STRINGSET")) {
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
        std::cout << "size of Fock space=" << prototype->size() << std::endl;
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
        action.operatorOnWavefunction(*ham, w);
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

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::prepareGuess() { }

template<>
void Davidson<MixedWavefunction, MixedOperatorSecondQuant>::prepareGuess() {
    // Currently assumes only 1 mode
    auto prof = profiler->push("prepareGuess");
    if (GA_Nodeid() == 0)
        std::cout << "Entered Davidson::prepareGuess()" << std::endl;
    auto nS = nState;
    auto nMode = options.parameter("NMODE", int(0));
    Wavefunction w = prototype->wavefunctionAt(0, mpi_comm_compute);
    if (nMode != 1) return;
    auto nM = options.parameter("NMODAL", int(0));
    auto n = nS / nM;
    n += (n * nM == nS) ? 0 : 1;
    n += (n == 1 && nS > 3 && w.size() > 3) ? 2 : 0;
    auto modOptions = Options(options);
    modOptions.addParameter("NSTATE", (int) n);
    // Modify options to choose the correct number of electronic states
    SymmetryMatrix::Operator *h = ham->elHam["Hel[0]"].get();
    Davidson<Wavefunction, SymmetryMatrix::Operator> elecSolver(std::move(w),
                                                                SymmetryMatrix::Operator{*h},
                                                                modOptions);
    elecSolver.run();
    // Loop over electronic states, loop over modals, set each element of the electronic wavefunction
    for (unsigned int root = 0; root < nState; root++) {
        auto ind_elec_state = root / nM;
        auto iVib = root - ind_elec_state * nM;
        ww[root].zero();
        ww[root].put(iVib, elecSolver.ww[ind_elec_state]);
    }
    MPI_Barrier(mpi_comm_compute);
    if (GA_Nodeid() == 0)
        std::cout << "Exit Davidson::prepareGuess()" << std::endl;
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::run() {
    auto prof = profiler->push("Davidson");
    message();
    initialize();
    prepareGuess();
//    printMatrix();
    for (unsigned int iteration = 1; iteration <= maxIterations; iteration++) {
        action();
        solver.addVector(ww, gg);
        update();
//        if (solver.endIteration(ww, gg) && iteration > 4) break;
        if (solver.endIteration(ww, gg)) break;
    }
    if (GA_Nodeid() == 0) {
        xout << "energies: ";
        for (unsigned int i = 0; i < nState; ++i) xout << solver.eigenvalues()[i] << ", ";
        xout << std::endl;
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::initialize() {
    if (!diagonalH) diagonalH = std::make_shared<t_Wavefunction>(*prototype, 0);
    diagonalH->allocate_buffer();
    diagonalH->diagonalOperator(*ham);
    //diagonalH->set(1.0);
    //std::vector<double> v(100);
    //for (auto i=0; i < 100; ++i){
    //    v[i] = diagonalH->at(i);
    //}
    //for (auto el: v){
    //    xout << el << ", " << std::endl;
    //}
    auto minLocs = diagonalH->minlocN(nState);
    std::vector<int> roots(nState, 0);
    for (unsigned int root = 0; root < nState; root++) {
        ww.emplace_back(*prototype, 0);
        ww.back().allocate_buffer();
        ww.back().zero();
        auto n = minLocs[root];
        if (GA_Nodeid() == 0)
            if (std::count(roots.begin(), roots.begin() + root, n) != 0)
                throw std::logic_error("Davidson::initialize duplicate guess vector, n =" + std::to_string(n));
        roots[root] = n;
        ww.back().set(n, 1.0);
        for (auto i =0; i <=root; ++i){
            std::cout << i << " " << root << " " <<  ww[i].dot(ww.back()) << std::endl;
        }
        gg.emplace_back(*prototype, 0);
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
            g.operatorOnWavefunction(*ham, x, false);
        }
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::update() {
    auto eigval = solver.eigenvalues();
    auto minLocs = diagonalH->minlocN(nState);
    for (size_t state = 0; state < nState; state++) {
        t_Wavefunction &cw = ww[state];
        const t_Wavefunction &gw = gg[state];
        auto shift = -eigval[state] + 1e-10;
        shift += 2 * std::numeric_limits<value_type>::epsilon() * std::max<value_type>(1, std::abs(
                diagonalH->at(minLocs[state]))); // to guard against zero
        cw.divide(&gw, diagonalH.get(), shift, true, true);
    }
}

}  // namespace run
}  // namespace gci


template
class gci::run::Davidson<gci::MixedWavefunction, gci::MixedOperatorSecondQuant>;

template
class gci::run::Davidson<gci::Wavefunction, SymmetryMatrix::Operator>;

