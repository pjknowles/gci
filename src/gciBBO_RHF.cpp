#include "gciRun.h"
#include "gciBBO_RHF.h"

#include <valarray>

#include "gciOperatorBBO.h"

namespace gci {

void Run::BBO_RHF(const State &prototype) {
    // Preamble
    profiler->start("BBO_RHF:preamble");
    int maxIterations = options.parameter("HF_MAX_ITER", std::vector<int>{50}).at(0);
    double thresh = options.parameter("HF_E_THRESH", std::vector<double>{1.0e-8}).at(0);
    OperatorBBO molHam(options); // molecular Hamiltonian
    std::vector<int> symmetries;
    for (auto s : prototype.orbitalSpace->orbital_symmetries) {symmetries.push_back(s + 1);}
    dim_t dim;
    for (auto s: *prototype.orbitalSpace) dim.push_back(s);
    dim_t occ(8, 0);
    std::vector<int> occInt = options.parameter("OCC", std::vector<int>{8});
    for (int i = 0; i < 8; ++i) {occ[i] = (size_t) occInt[i];}
    nm_BBO_RHF::Density density(dim, occ, symmetries, prototype.symmetry);
    std::vector<SMat> U; // vibrational modals
    for (int iMod = 0; iMod < molHam.m_nMode; ++iMod) {
        dim.assign(8, 0);
        dim[molHam.m_symMode[iMod] - 1] = (size_t) molHam.m_nModal;
        SMat Umat({dim, dim}, parityNone, molHam.m_symMode[iMod], "Modals for mode " + std::to_string(iMod + 1));
        Umat.setIdentity();
        U.push_back(Umat);
    }
    profiler->stop("BBO_RHF:preamble");

    profiler->start("BBO_RHF:algorithm");
    std::valarray<double> energy(4, 0), energyPrev(4, 0);
    nm_BBO_RHF::writeFormat();
    // Energy and convergence check first. Fock matrix evaluation second.
    for (int iIter = 0; iIter < maxIterations; ++iIter) {
        density.update();
        molHam.energy(density.P, U, energy);
        nm_BBO_RHF::writeIter(iIter, energy, energyPrev);
        auto res = std::abs(energy - energyPrev);
        if (res.max() < thresh) {
            xout << "Convergence achieved " << std::endl;
            break;
        };
        energyPrev = energy;
        nm_BBO_RHF::solveFock(molHam, density, U);
    }
    profiler->stop("RHF");
}

void nm_BBO_RHF::writeFormat() {
    std::cout << "Iter " << " Etot " << " dE " << std::endl;
    std::cout << "     " << " Eel " << " dE " << std::endl;
    std::cout << "     " << " Evib " << " dE " << std::endl;
    std::cout << "     " << " Eint " << " dE " << std::endl;
}

void nm_BBO_RHF::writeIter(int iIter, std::valarray<double> &energy, std::valarray<double> &energyPrev) {
    std::cout << iIter << "    " << energy[0] << "    " << energy[0] - energyPrev[0] << std::endl;
    std::cout << "             " << energy[1] << "    " << energy[1] - energyPrev[1] << std::endl;
    std::cout << "             " << energy[2] << "    " << energy[2] - energyPrev[2] << std::endl;
    std::cout << "             " << energy[3] << "    " << energy[3] - energyPrev[3] << std::endl;
}

void nm_BBO_RHF::solveFock(OperatorBBO &molHam, Density &density, std::vector<SMat> &U) {
    {
        // Electronic degrees of freedom
        Operator F = molHam.electronicFock(density.P, U);
        SMat eigVal({density.dim}, parityNone, -1, "Eigenvalues");
        F.O1().ev(eigVal, &density.Cmat);
        density.update();
    }
    // Vibrational degrees of freedom
    for (int iMode = 0; iMode < molHam.m_nMode; ++iMode) {
        Operator F = molHam.vibrationalFock(density.P, U[iMode], iMode);
        SMat eigVal({U[iMode].dimensions()[1]}, parityNone, -1, "Eigenvalues");
        F.O1().ev(eigVal, &U[iMode]);
    }
}

} // namespace gci
