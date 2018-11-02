#include "gciOperatorBBO.h"
#include "gciConstants.h"
#include "gciRun.h"

#include <utility>
#include <numeric>
#include <FCIdump.h>

namespace gci {

OperatorBBO::OperatorBBO(const Options &options, std::string description) :
        m_description(std::move(description)),
        m_nMode(options.parameter("NMODE", (int) 0)),
        m_nModal(options.parameter("NMODAL", (int) 0)),
        m_vibOcc(options.parameter("VIBOCC", std::vector<int>(m_nMode, 0))),
        m_symMode(options.parameter("SYMMODE", std::vector<int>(m_nMode, 1))),
        m_freq(options.parameter("FREQ", std::vector<double>(m_nMode, 0))),
        m_nmDisp(options.parameter("NM_DISP", std::vector<double>(m_nMode, 0))),
        m_fcidump(options.parameter("FCIDUMP", std::vector<std::string>(1, "")).at(0)),
        m_Hel(constructOperator(FCIdump(m_fcidump))) {
    for (auto iter = m_freq.begin(); iter < m_freq.end(); ++iter) *iter *= Constants::CM_TO_AU;
    std::cout << "FREQ = ";
    for (auto iter = m_freq.begin(); iter < m_freq.end(); ++iter) std::cout << *iter << ", ";
    std::cout << std::endl;
    m_Hel.m_description = "Electronic Hamiltonian at the reference geometry";
    // initiate interaction and vibrational Hamiltonian
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        dim_t dimension(8);
        dimension[m_symMode[iMode] - 1] = (unsigned int) m_nModal;
        std::vector<int> orbital_symmetries(m_nModal, m_symMode[iMode]);
        std::string descr;
        descr = "Vibrational Hamiltonian: mode = " + std::to_string(iMode);
        Operator hVib(dimension, 1, false, (unsigned int) m_symMode[iMode] - 1, true, descr);
        descr = "Vibrational component of Interaction Hamiltonian: mode = " + std::to_string(iMode);
        Operator hIntVib(dimension, 1, false, (unsigned int) m_symMode[iMode] - 1, true, descr);
        hVib.O1(true).assign(0.0);
        hIntVib.O1(true).assign(0.0);
        for (int jModal = 0; jModal < m_nModal; ++jModal) {
            double *elVib_jj = &(hVib.element(jModal, m_symMode[iMode] - 1, jModal, m_symMode[iMode] - 1));
            *elVib_jj = m_freq[iMode] * (jModal + 0.5);
            double *elIntVib_jj = &(hIntVib.element(jModal, m_symMode[iMode] - 1, jModal, m_symMode[iMode] - 1));
            *elIntVib_jj = -m_nmDisp[iMode];
            if (jModal > 0) {
                double *elIntVib_jjm1 = &(hIntVib.element(jModal, m_symMode[iMode] - 1, jModal - 1,
                                                          m_symMode[iMode] - 1));
                *elIntVib_jjm1 = std::sqrt((double) jModal) / std::sqrt(2.0 * m_freq[iMode]);
            }
        }
        std::string fcidumpN = m_fcidump + std::to_string(iMode + 1);
        Operator hIntEl(constructOperator(FCIdump(fcidumpN)));
        hIntEl.m_O0 -= m_freq[iMode] * m_nmDisp[iMode];
        m_Hel.m_O0 -= 0.5 * m_freq[iMode] * std::pow(m_nmDisp[iMode], 2.0);
//        xout << hVib << std::endl;
//        xout << hIntVib << std::endl;
        m_Hvib.push_back(hVib);
        m_HintEl.push_back(hIntEl);
        m_HintVib.push_back(hIntVib);
    }
}


double OperatorBBO::transformedVibHamElement(const Operator &hamiltonian, const SMat &U, int r, int s, int symm) {
    dim_t dim(8, 0);
    dim[symm] = 1;
    SMat Usplice({U.dimensions()[1], dim}, parityNone, symm);
    SMat UspliceT({U.dimensions()[1], dim}, parityNone, symm);
    dim_t offset(8, 0);
    dim_t zeroOffset(8, 0);
    offset[symm] = (size_t) s;
    Usplice.splice(U, {zeroOffset, offset});
    offset.assign(8, 0);
    offset[symm] = (size_t) r;
    UspliceT.splice(U, {zeroOffset, offset});
    UspliceT.transpose();
    SMat el = (UspliceT * hamiltonian.O1(true) * Usplice);
    return el.block((unsigned) symm)[0];
}

Operator OperatorBBO::transformedVibHam(const Operator &hamiltonian, const SMat &U) {
    Operator H(hamiltonian);
    H.O1(true) = SymmetryMatrix::transpose(U) * hamiltonian.O1(true) * U;
    return H;
}

void OperatorBBO::energy(const Operator &density, const std::vector<SMat> &U, std::valarray<double> &energy) {
    energy = 0;
    // Pure electronic energy
    nm_RHF::electronicEnergy(density, m_Hel, energy[1]);
    // Pure vibrational energy
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        energy[2 + iMode] += transformedVibHamElement(m_Hvib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode],
                                                      m_symMode[iMode] - 1);
//        xout << m_Hvib[iMode] << std::endl;
//        xout << transformedVibHam(m_Hvib[iMode], U[iMode]);
    }
    // Interaction energy
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        double d, o;
        nm_RHF::electronicEnergy(density, m_HintEl[iMode], d);
        o = transformedVibHamElement(m_HintVib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode],
                                     m_symMode[iMode] - 1);
        energy[2 + m_nMode + iMode] += d * o;
//        xout << "d " << d << ", o " << o << std::endl;
    }
    energy[0] = std::accumulate(std::begin(energy) + 1, std::end(energy), 0.0);
}

Operator OperatorBBO::electronicFock(const Operator &P, std::vector<SMat> &U, const SMat Cmat) {
    Operator F = m_Hel.fock(P, true, "Fock operator");
    std::cout << SymmetryMatrix::transpose(Cmat) * F.O1(true) * Cmat << std::endl;
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        Operator Fint = m_HintEl[iMode].fock(P, true, "interaction component of Fock operator");
        double o = transformedVibHamElement(m_HintVib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode],
                                            m_symMode[iMode] - 1);
//        xout << Fint << std::endl;
        F.O1(true) += o * Fint.O1(true);
        std::cout << "mode " << iMode << std::endl;
        std::cout << SymmetryMatrix::transpose(Cmat) * Fint.O1(true) * Cmat << std::endl;
    }
    return F;
}

Operator OperatorBBO::vibrationalFock(const Operator &P, const SMat &U, int iMode) {
    dim_t dimension(8);
    dimension[m_symMode[iMode] - 1] = (unsigned int) m_nModal;
    std::vector<int> orbital_symmetries(m_nModal, m_symMode[iMode]);
    Operator F(dimension, 1, false, (unsigned) m_symMode[iMode] - 1, true,
               "Vibrational Fock matrix, mode = " + std::to_string(iMode));
//    F.O1(true) = transformedVibHam(m_Hvib[iMode], U).O1(true);
    F.O1(true) = m_Hvib[iMode].O1(true);
    double d;
    nm_RHF::electronicEnergy(P, m_HintEl[iMode], d);
    F.O1(true) += 0.5 * d * m_HintVib[iMode].O1(true);
//    xout << F << std::endl;
    return F;
}

std::ostream &operator<<(std::ostream &os, const OperatorBBO &obj) {
    os << "OperatorBBO: " << obj.m_description << std::endl;
    os << "nMode = " << obj.m_nMode << std::endl;
    os << "nModal = " << obj.m_nModal << std::endl;
    os << "vibOcc = ";
    for (auto iter = obj.m_vibOcc.begin(); iter < obj.m_vibOcc.end(); ++iter) os << *iter << " ";
    os << std::endl << "freq = ";
    for (auto iter = obj.m_freq.begin(); iter < obj.m_freq.end(); ++iter) os << *iter << " ";
    os << std::endl;
//    os << obj.m_Hel << std::endl;
    for (int iMode = 0; iMode < obj.m_nMode; ++iMode) {
        os << obj.m_Hvib[iMode] << std::endl;
//        os << obj.m_HintEl[iMode] << std::endl;
        os << obj.m_HintVib[iMode] << std::endl;
    }
    return os;
}

} // namespace gci
