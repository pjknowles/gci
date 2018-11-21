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
          energyThreshold(options.parameter("TOL", 1e-8)),
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
    xout << std::fixed << std::setprecision(8);
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::run() {
    auto prof = profiler->push("Davidson");
    message();
    initialize();
    for (auto iteration = 1; iteration <= (size_t) maxIterations; iteration++) {
        action();
        solver.addVector(ww, gg, active);
        update();
        if (solver.endIteration(ww, gg, active)) break;
    }
    if (solver.m_verbosity > 0)
        xout << "Number of actions of matrix on vector = " << solver.m_actions << std::endl;
//    for (auto root = 0; root < options.nState; root++) {
//        m_wavefunctions.push_back(std::make_shared<Wavefunction>(ww[root]));
//        m_wavefunctions.back()->m_properties["ENERGY"] = solver.eigenvalues()[root];
//      if (options.parameter("DENSITY",0)>0)
//        m_wavefunctions.back()->density = m_wavefunctions.back()->density(options.parameter("DENSITY",0));
//    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::initialize() {
    if (!diagonalH) diagonalH = std::make_shared<t_Wavefunction>(*prototype, 0);
    diagonalH->allocate_buffer();
    diagonalH->diagonalOperator(*ham);
    for (int root = 0; root < nState; root++) {
        active.emplace_back(true);
        ww.emplace_back(*prototype);
        ww.back().allocate_buffer();
        ww.back().zero();
        ww.back().set(diagonalH->minloc(root + 1), 1.0);
        gg.emplace_back(*prototype);
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
        if (active[k]) {
            auto prof = profiler->push("Hc");
            g.operatorOnWavefunction(*ham, x, parallel_stringset);
        }
    }
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::update() {
    for (size_t state = 0; state < nState; state++) {
        t_Wavefunction &cw = ww[state];
        const t_Wavefunction &gw = gg[state];
        auto shift = -solver.eigenvalues()[state] + 1e-10;
        shift += 2 * std::numeric_limits<value_type>::epsilon() * std::max<value_type>(1, std::abs(
                diagonalH->at(diagonalH->minloc(state + 1)))); // to guard against zero
        cw.divide(&gw, diagonalH.get(), shift, true, true);
    }
}


}  // namespace run
}  // namespace gci


