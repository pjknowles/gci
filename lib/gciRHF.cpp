#include "gciRun.h"
#include "gciRHF.h"

double Run::RHF(const SymmetryMatrix::Operator &hamiltonian, const State &prototype,
                double thresh, int maxIterations) {
    profiler->start("RHF:preamble");
    std::vector<int> symmetries;
    for (const auto &s : prototype.orbitalSpace->orbital_symmetries) {
        symmetries.push_back(s + 1);
    }
    SymmetryMatrix::dim_t dim;
    for (const auto s: *prototype.orbitalSpace) dim.push_back(s);
    SymmetryMatrix::Operator C(dim, 1, false, prototype.symmetry, false, "MO");
    SymmetryMatrix::dim_t occ(8, 0);
    std::vector<int> occInt(8, 0);
    occInt = options.parameter("OCC", std::vector<int>{8});
    for (unsigned j = 0; j < 8; ++j) occ[j] = (size_t) occInt[j];
    xout << "OCC ";
    for (int j = 0; j < 8; ++j) {xout << occ[j] << ",";}
    xout << std::endl;
    nm_RHF::Density density(dim, occ, symmetries, prototype.symmetry);
    profiler->stop("RHF:preamble");
    profiler->start("RHF:algorithm");
    double energy, energyPrev = 0.0, err;
    std::cout << "Iter " << " E " << " dE " << " res " << std::endl;
    for (int iIter = 0; iIter < maxIterations; ++iIter) {
        density.update();
//        energy = nm_RHF::electronicEnergy(density.P, hamiltonian);
        auto F = nm_RHF::electronicEnergy(density.P, hamiltonian, energy);
        // identity FP == PF
        auto res = (F.O1() * density.P.O1()) - (density.P.O1() * F.O1());
        err = 0.0;
        for (unsigned jSym = 0; jSym < 8; ++jSym) {
            if (res.dimension(jSym) == 0) continue;
            auto r = std::max_element(res.block(jSym).begin(), res.block(jSym).end());
            if (abs(*r) > err) err = *r;
        }
        std::cout << iIter << "    " << energy << "    " << energy - energyPrev << "    " << err << std::endl;
        if (err < thresh) break;
        energyPrev = energy;
        SymmetryMatrix::SMat eigVal({dim}, SymmetryMatrix::parityNone, -1, "Eigenvalues");
        F.O1().ev(eigVal, &density.Cmat);
    }
    profiler->stop("RHF:algorithm");
}

nm_RHF::Density::Density(SymmetryMatrix::dim_t &dim, SymmetryMatrix::dim_t &occ, std::vector<int> &symmetries,
                         int symmetry) :
        dim(dim), occ(occ), symmetries(symmetries), symmetry(symmetry),
        Cmat(SymmetryMatrix::SMat({dim, dim}, SymmetryMatrix::parityNone)),
        Csplice(SymmetryMatrix::SMat({dim, occ}, SymmetryMatrix::parityNone)),
        P(SymmetryMatrix::Operator(dim, 1, false, (unsigned) symmetry, true, "density")) {
    Cmat.setIdentity();
    update();
}

void nm_RHF::Density::update() {
    Csplice.splice(Cmat);
    P.O1(true) = 2 * (Csplice * SymmetryMatrix::transpose(Csplice));
}

SymmetryMatrix::Operator nm_RHF::electronicEnergy(const SymmetryMatrix::Operator &P, const SymmetryMatrix::Operator &Hel, double &energy) {
    auto F = Hel.fock(P, true, "Fock operator");
    energy = Hel.m_O0 + 0.5 * (P.O1() & (Hel.O1() + F.O1()));
    return F;
};
