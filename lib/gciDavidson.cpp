#include "gciDavidson.h"

#include <iomanip>

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
    std::cout << "Davidson eigensolver, maximum iterations=" << maxIterations;
    std::cout << "; number of states=" << nState;
    std::cout << "; energy threshold=" << std::scientific << std::setprecision(1) << energyThreshold << std::endl;
    std::cout << "size of Fock space=" << prototype->size() << std::endl;
    xout << std::fixed << std::setprecision(12);
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::printMatrix() { }

template<>
void Davidson<MixedWavefunction, MixedOperator>::printMatrix() {
    MixedWavefunction w(ww[0]);
    MixedWavefunction action(gg[0]);
    w.allocate_buffer();
    action.allocate_buffer();
    auto n = w.size();
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0));
    std::cout << "Hamiltonian matrix (nxn):" << std::endl;
    std::cout << "n = " << n << std::endl;
    std::cout << " H = {";
    for (int i = 0; i < n; ++i) {
        w.zero();
        w.set(i, 1.0);
        action.zero();
        action.operatorOnWavefunction(*ham, w);
        std::cout << "{";
        for (int j = 0; j < n; ++j) {
            H[i][j] = action.at(j);
            if (j == n - 1) std::cout << H[i][j] << "";
            else std::cout << H[i][j] << ",";
        }
        if (i == n - 1) std::cout << "}};" << std::endl;
        else std::cout << "}," << std::endl;
    }
}

template<class t_Wavefunction>
void write_vec(const t_Wavefunction &wfn, const std::string &message) { }

template<>
void write_vec<MixedWavefunction>(const MixedWavefunction &w, const std::string &message) {
    size_t n = w.size();
    std::cout << message << "  ";
    for (int i = 0; i < n; ++i) {
        std::cout << w.at(i) << " ";
    }
    std::cout << std::endl;
}

//template<class t_Wavefunction, class t_Operator>
//void printHFmatrixElements(const t_Wavefunction &wfn, const t_Operator &ham) { }

template<class t_Wavefunction, class t_Operator>
void printHFmatrixElements(const t_Wavefunction &wfn, const t_Operator &ham, int n) { }

/*!
 * @brief Prints expectation values of all terms in the Hamiltonian for specified HF determinant
 * @param wfn
 * @param ham
 * @param n Index of electronic determinant
 */
template<>
void printHFmatrixElements<MixedWavefunction, MixedOperator>(const MixedWavefunction &wfn, const MixedOperator &ham, int n) {
    MixedOperator hamMod(ham);
    for (auto &op : hamMod.Hmix[VibOpType::Qsq]) {
        if (op.vibOp.mode[0] != op.vibOp.mode[1]) continue;
        auto i = op.vibOp.mode[0];
        op.Hel.m_O0 += std::pow(hamMod.freq[i], 2);
    }
    auto matEls = wfn.hfMatElems(hamMod, n);
    std::cout << "Matrix elements over Determinant n = " << n << std::endl;
    for (const auto &el : matEls) {
        std::cout << el.first << " = " << el.second << std::endl;
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::run() {
    auto prof = profiler->push("Davidson");
    message();
    initialize();
//    printMatrix();
//    printHFmatrixElements(*prototype, *ham, 0);
//    printHFmatrixElements(*prototype, *ham, 1);
//    printHFmatrixElements(*prototype, *ham, 2);
//    printHFmatrixElements(*prototype, *ham, 3);
//    printHFmatrixElements(*prototype, *ham, 4);
    for (auto iteration = 1; iteration <= (size_t) maxIterations; iteration++) {
        action();
        solver.addVector(ww, gg, active);
        update();
        if (solver.endIteration(ww, gg, active)) break;
    }
    if (solver.m_verbosity > 0)
        xout << "Number of actions of matrix on vector = " << solver.m_actions << std::endl;
//    for (auto root = 0; root < nState; root++) {
//        write_vec(ww[root],
//                  "eigenvector " + std::to_string(root) + " active=" + std::to_string(active[root]) + " converged=" +
//                  std::to_string(solver.errors()[root]) + ":");
//    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::initialize() {
    if (!diagonalH) diagonalH = std::make_shared<t_Wavefunction>(*prototype, 0);
    diagonalH->allocate_buffer();
    diagonalH->diagonalOperator(*ham);
    std::vector<int> roots(nState, 0);
    for (int root = 0; root < nState; root++) {
        active.emplace_back(true);
        ww.push_back(t_Wavefunction(*prototype, 0));
        ww.back().allocate_buffer();
        ww.back().zero();
        auto n = diagonalH->minloc(root + 1);
        if (std::count(roots.begin(), roots.begin() + root, n) != 0)
            throw std::logic_error("Davidson::initialize duplicate guess vector, n =" + std::to_string(n));
        roots[root] = n;
        ww.back().set(n, 1.0);
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
//        write_vec<t_Wavefunction>(x, "eigenvector applied on Hc " + std::to_string(k) + " before update");
        t_Wavefunction &g = gg[k];
        g.zero();
        if (active[k]) {
            auto prof = profiler->push("Hc");
            g.operatorOnWavefunction(*ham, x, parallel_stringset);
        }
//        write_vec<t_Wavefunction>(g, "residual berfore addVector " + std::to_string(k) + " ");
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::update() {
    auto eigval = solver.eigenvalues();
    for (size_t state = 0; state < nState; state++) {
        t_Wavefunction &cw = ww[state];
        const t_Wavefunction &gw = gg[state];
//        write_vec<t_Wavefunction>(cw, "eigenvector " + std::to_string(state) + " before update");
//        write_vec<t_Wavefunction>(gw, "residual " + std::to_string(state) + " ");
        auto shift = -eigval[state] + 1e-10;
        shift += 2 * std::numeric_limits<value_type>::epsilon() * std::max<value_type>(1, std::abs(
                diagonalH->at(diagonalH->minloc(state + 1)))); // to guard against zero
        cw.divide(&gw, diagonalH.get(), shift, true, true);
//        write_vec<t_Wavefunction>(cw, "eigenvector " + std::to_string(state) + " after update");
    }
}


}  // namespace run
}  // namespace gci


