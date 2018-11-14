#include "gciDavidson.h"

#include <iomanip>

#include <Operator.h>
#include <IterativeSolver.h>
#include "gciWavefunction.h"
#include "gciMixedWavefunction.h"
#include "gciMixedOperator.h"

namespace gci {
namespace run {


std::ostream &operator<<(std::ostream &os, DavidsonOptions const &obj) {
    os << "Davidson eigensolver, maximum iterations=" << obj.maxIterations;
    os << "; number of states=" << obj.nState;
    os << "; energy threshold=" << std::scientific << std::setprecision(1) << obj.energyThreshold << std::endl;
    return os;
}

template<class t_Wavefunction, class t_Operator>
Davidson<t_Wavefunction, t_Operator>::Davidson(const std::shared_ptr<t_Wavefunction> &prototype,
                                               const std::shared_ptr<t_Operator> &ham, Options &opt)
        : prototype(prototype), ham(ham), options(opt), solver(),
          Pcoeff(), P(), Presid(ham, P) {
    solver.m_thresh = options.energyThreshold;
    solver.m_verbosity = options.solverVerbosity;
    solver.m_maxIterations = (unsigned int) options.maxIterations;
    solver.m_roots = (size_t) options.nState;
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::run() {
    auto prof = profiler->push("Davidson");
    xout << options << std::endl;
    Presidual<t_Wavefunction, t_Operator> Presid(ham, P);
    initialize();

    for (size_t p = 0; p < initialNP; p++) {
        auto det1 = diagonalH.minloc(p + 1);
        P.emplace_back(Pvector{{det1, 1}});
    }
    for (auto iteration = 1; iteration <= (size_t) options.maxIterations; iteration++) {
        //TODO I don't understand why one would want to rebuild the P space at every iteration.
        if (options.pSpaceRebuild && iteration != 1) { // clear out the P space and rebuild it
            solver.clearP();
            P.clear();
            Pcoeff.clear();
        }
        buildPspace();
        Presid(Pcoeff, gg); // augment residual with contributions from P space
        std::vector<double> shift;
        for (auto root = 0; root < options.nState; root++) {
            shift.push_back(-solver.eigenvalues()[root] + 1e-10);
        }
        update(ww, gg, shift, true);
        if (solver.endIteration(ww, gg, active)) break;
        resid(ww, gg, active);
        solver.addVector(ww, gg, active, Pcoeff);
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
    if (!diagonalH) diagonalH = std::make_shared<t_Wavefunction>(*prototype);
    diagonalH->diagonalOperator(ham);
    update = updater<t_Wavefunction>(*diagonalH, false);
    resid = residual<t_Wavefunction, t_Operator>(*ham, false);
    for (int root = 0; root < options.nState; root++) {
        active.push_back(true);
        ww.emplace_back(prototype);
        ww.back().allocate_buffer();
        gg.emplace_back(prototype);
        gg.back().allocate_buffer();
        //FIXME Tile size should have been set in prototype already!
        gg.back().settilesize(
                options.global().parameter("TILESIZE", std::vector<int>(1, -1)).at(0),
                options.global().parameter("ALPHATILESIZE", std::vector<int>(1, -1)).at(0),
                options.global().parameter("BETATILESIZE", std::vector<int>(1, -1)).at(0));
    }
    initialNP = (int) std::min((size_t) options.pSpaceInitial, ww.front().size());
    maxNP = (int) std::min((size_t) options.pSpace, ww.front().size());
}

template<class t_Wavefunction, class t_Operator>
void Davidson<t_Wavefunction, t_Operator>::buildPspace() {
    auto nP = P.size();
    if (maxNP <= nP) return;
    Presid(Pcoeff, gg); // augment residual with contributions from P space
    auto newP = solver.suggestP(ww, gg, active, (maxNP - nP));
    for (const auto &pp : P) {
        newP.erase(std::remove(newP.begin(), newP.end(), pp.begin()->first),
                   newP.end());
    } // remove anything already in P
    for (const auto &det1 : newP) {
        P.emplace_back(Pvector{{det1, 1}});
    }
    const auto newNP = P.size();
    if (solver.m_verbosity > 1 && newNP > nP)
        xout << "Adding " << newNP - nP << " P-space configurations (total " << newNP << ")" << std::endl;
    std::vector<double> addHPP(newNP * (newNP - nP), (double) 0);
    for (size_t p0 = nP; p0 < newNP; p0++) {
        Wavefunction wsparse(prototype);
        wsparse.m_sparse = true;
        Wavefunction gsparse(prototype);
        gsparse.m_sparse = true;
        wsparse.set(P[p0].begin()->first, (double) 1);
        gsparse.operatorOnWavefunction(ham, wsparse);
        for (size_t p1 = 0; p1 < newNP; p1++) {
            auto jdet1 = P[p1].begin()->first;
            if (gsparse.buffer_sparse.count(jdet1))
                addHPP[p1 + (p0 - nP) * newNP] = gsparse.buffer_sparse.at(jdet1);
        }
    }
    Pcoeff.resize(newNP);
    solver.addP(std::vector<Pvector>(P.begin() + nP, P.end()), addHPP.data(), ww, gg, Pcoeff);
}

}  // namespace run
}  // namespace gci


