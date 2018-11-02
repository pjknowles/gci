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
        SMat Umat({dim, dim}, parityNone, 0, "Modals for mode " + std::to_string(iMod + 1));
//        SMat Umat({dim, dim}, parityNone, molHam.m_symMode[iMod], "Modals for mode " + std::to_string(iMod + 1));
        Umat.setIdentity();
//        nm_BBO_RHF::rotate(Umat, 0.5);
        U.push_back(Umat);
    }
//    nm_BBO_RHF::rotate(density.Cmat, 0.5);
    density.update();
    profiler->stop("BBO_RHF:preamble");

    profiler->start("BBO_RHF:algorithm");
    std::valarray<double> energy(0.0, (size_t) 2 + molHam.m_nMode * 2), energyPrev(0.0,
                                                                                   (size_t) 2 + molHam.m_nMode * 2);
    nm_BBO_RHF::writeFormat();
    // Energy and convergence check first. Fock matrix evaluation second.
    for (int iIter = 0; iIter < maxIterations; ++iIter) {
        density.update();
        molHam.energy(density.P, U, energy);
        nm_BBO_RHF::writeIter(iIter, energy, energyPrev, molHam.m_nMode);
        std::valarray<double> res(0.0, energy.size());
        res = std::abs(energy - energyPrev);
        if (res.max() < thresh && iIter >= 5) {
            xout << "Convergence achieved " << std::endl;
            break;
        };
        energyPrev = energy;
        nm_BBO_RHF::solveFock(molHam, density, U);
    }
    profiler->stop("BBO_RHF:algorithm");
// Let's change the occupancy of second mode and recalculate the energy
    if (false) {
        molHam.m_vibOcc[1] = 1;
        molHam.energy(density.P, U, energy);
        nm_BBO_RHF::writeIter(-1, energy, energyPrev, molHam.m_nMode);
        for (int iIter = 0; iIter < maxIterations; ++iIter) {
            density.update();
            molHam.energy(density.P, U, energy);
            nm_BBO_RHF::writeIter(iIter, energy, energyPrev, molHam.m_nMode);
            std::valarray<double> res(0.0, energy.size());
            res = std::abs(energy - energyPrev);
            if (res.max() < thresh) {
                xout << "Convergence achieved " << std::endl;
                break;
            };
            energyPrev = energy;
            nm_BBO_RHF::solveFock(molHam, density, U);
        }
    }
}

void nm_BBO_RHF::writeFormat() {
    xout << "Iter " << " Etot " << " dE " << std::endl;
    xout << "     " << " Eel " << " dE " << std::endl;
    xout << "     " << " Evib " << std::endl;
    xout << "     " << " diff " << std::endl;
    xout << "     " << " Eint " << std::endl;
    xout << "     " << " diff " << std::endl;
    xout << std::endl;
}

void nm_BBO_RHF::writeIter(int iIter, std::valarray<double> &energy, std::valarray<double> &energyPrev, int nMode) {
    xout << iIter << "     Etot " << energy[0] << ", " << energy[0] - energyPrev[0] << ", " << std::endl;
    xout << "      Eel " << energy[1] << ", " << energy[1] - energyPrev[1] << ", " << std::endl;
//    std::cout << "      " << energy[2] << "    " << energy[2] - energyPrev[2] << std::endl;
//    std::cout << "      " << energy[3] << "    " << energy[3] - energyPrev[3] << std::endl;
    std::valarray<double> diff = energy - energyPrev;
    xout << "      Evib ";
    for (auto iter = std::begin(energy) + 2; iter < std::begin(energy) + 2 + nMode; ++iter) {
        xout << *iter << ", ";
    }
    xout << std::endl;
    xout << "      diff ";
    for (auto iter = std::begin(diff) + 2; iter < std::begin(diff) + 2 + nMode; ++iter) {
        xout << *iter << ", ";
    }
    xout << std::endl;
    xout << "      Eint ";
    for (auto iter = std::begin(energy) + nMode + 2; iter < std::begin(energy) + 2 + 2 * nMode; ++iter) {
        xout << *iter << ", ";
    }
    xout << std::endl;
    xout << "      diff ";
    for (auto iter = std::begin(diff) + nMode + 2; iter < std::begin(diff) + 2 + 2 * nMode; ++iter) {
        xout << *iter << ", ";
    }
    xout << "\n" << std::endl;
}

void nm_BBO_RHF::solveFock(OperatorBBO &molHam, Density &density, std::vector<SMat> &U) {
    {
        // Electronic degrees of freedom
        Operator F = molHam.electronicFock(density.P, U, density.Cmat);
//        xout << SymmetryMatrix::transpose(density.Cmat) * F.O1(true) * density.Cmat << std::endl;
        SMat eigVal({density.dim}, parityNone, -1, "Eigenvalues");
        SMat cMat(density.Cmat);
        F.O1().ev(eigVal, &cMat /*&density.Cmat*/, nullptr, nullptr, "lapack", "ascending");
        density.Cmat = cMat;
//        xout << "Cmat " << cMat << std::endl;
        density.update();
//        xout << SymmetryMatrix::transpose(cMat) * F.O1(true) * cMat << std::endl;
//        xout << density.Cmat << std::endl;
    }
    // Vibrational degrees of freedom
    for (int iMode = 0; iMode < molHam.m_nMode; ++iMode) {
        Operator F = molHam.vibrationalFock(density.P, U[iMode], iMode);
        SMat eigVal({U[iMode].dimensions()[1]}, parityNone, -1, "Eigenvalues");
        F.O1().ev(eigVal, &U[iMode]);
//        xout << eigVal.m_description;
//        xout << SymmetryMatrix::transpose(U[iMode]) * F.O1(true) * U[iMode] << std::endl;
    }
//    if(true){
//        xout << molHam.transformedVibHam(molHam.m_Hvib[1], U[1]);
//    }
}

void nm_BBO_RHF::rotate(SMat &mat, double ang) {
// For now lets do it only for the modals. The
    SMat rot(mat);
    rot.setIdentity();
    int sym = rot.symmetry(), n = (int) rot.dimension(sym);
    std::vector<std::valarray<double> > vecs;
    for (int i = 0; i < n; ++i) {
        vecs.emplace_back(0.0, n);
        vecs[i][i] = 1.0;
    }
    for (int i = 0; i < n - 1; ++i) {
        vecs[i] = std::cos(ang) * vecs[i] + std::sin(ang) * vecs[i + 1];
        vecs[i + 1] = -std::sin(ang) * vecs[i] + std::cos(ang) * vecs[i + 1];
    }
// Gram-Schmidt to orthonormalise
    auto norm = [](std::valarray<double> &v1, std::valarray<double> &v2) {return (v1 * v2).sum();};
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            double d = norm(vecs[i], vecs[j]);
            vecs[i] -= d * vecs[j];
        }
        double d = norm(vecs[i], vecs[i]);
        vecs[i] /= std::sqrt(d);
    }
    memory::vector<double> v = rot.block(sym);
    for (int i = 0; i < n; ++i) {
        for (int jj = i * n, j = 0; jj < (i + 1) * n; ++jj, ++j) {
            v[jj] = vecs[i][j];
        }
    }
//    for(auto i: v) xout << ", " << i;
//    xout << std::endl;
    /*
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double d = norm(vecs[i], vecs[j]);
            xout << i << " " << j << " " << d << std::endl;
        }
    }
    for (int i = 0; i < n; ++i) {
        xout << vecs[0].size();
        for (auto &v:vecs[i]) {xout << " , " << v;}
        xout << std::endl;
    }
    */
}


} // namespace gci
