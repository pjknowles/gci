#include "gciOperatorBBO.h"
#include "gciConstants.h"
#include <utility>

#include "FCIdump.h"

namespace gci {

OperatorBBO::OperatorBBO(Options &options, std::string description) :
        m_description(std::move(description)),
        m_nMode(options.parameter("NMODE", (int) 0)),
        m_nModal(options.parameter("NMODAL", (int) 0)),
        m_vibOcc(options.parameter("VIBOCC", std::vector<int>(m_nMode, 0))),
        m_symMode(options.parameter("SYMMODE", std::vector<int>(m_nMode, 1))),
        m_freq(options.parameter("FREQ", std::vector<double>(m_nMode, 0))),
        m_fcidump(options.parameter("FCIDUMP", std::vector<std::string>(1, "")).at(0)),
        m_Hel(Operator::construct(FCIdump(m_fcidump))) {
    for (auto iter = m_freq.begin(); iter < m_freq.end(); ++iter) *iter *= Constants::CM_TO_AU;
    m_Hel.m_description = "Electronic Hamiltonian at the reference geometry";
    // initiate interaction and vibrational Hamiltonian
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        dim_t dimension(8);
        dimension[m_symMode[iMode] - 1] = (unsigned int) m_nModal;
        std::vector<int> orbital_symmetries(m_nModal, m_symMode[iMode]);
        std::string descr;
        descr = "Vibrational Hamiltonian: mode = " + std::to_string(iMode);
        Operator hVib(dimension, orbital_symmetries, 1, false, (unsigned int) m_symMode[iMode] - 1, true, true, descr);
        descr = "Vibrational component of Interaction Hamiltonian: mode = " + std::to_string(iMode);
        Operator hIntVib(dimension, orbital_symmetries, 1, false, (unsigned int) m_symMode[iMode] - 1, true, true,
                         descr);
        hVib.O1(true).assign(0.0);
        hIntVib.O1(true).assign(0.0);
        for (int jModal = 0; jModal < m_nModal; ++jModal) {
            double *v = &(hVib.O1(true).block((unsigned) m_symMode[iMode] - 1)[jModal * (jModal + 1) / 2 + jModal]);
            v[0] = m_freq[iMode] * (jModal + 0.5);
            if (jModal > 0) {
                v = &(hIntVib.O1(true).block((unsigned) m_symMode[iMode] - 1)[jModal * (jModal + 1) / 2 + jModal - 1]);
                v[0] = std::sqrt((double) jModal) / std::sqrt(2.0 * m_freq[iMode]);
            }
        }
        std::string fcidumpN = m_fcidump + std::to_string(iMode + 1);
        Operator hIntEl(Operator::construct(FCIdump(fcidumpN)));
        m_Hvib.push_back(hVib);
        m_HintEl.push_back(hIntEl);
        m_HintVib.push_back(hIntVib);
    }
}

void OperatorBBO::energy(Operator &density, std::vector<SMat> &U, std::valarray<double> &energy) {
    energy = 0;
    // Pure electronic energy
    nm_RHF::electronicEnergy(density, m_Hel, energy[1]);
    // Pure vibrational energy
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        energy[2] = transformedVibHamElement(m_Hvib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode]);
    }
    // Interaction energy
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        double d, o;
        nm_RHF::electronicEnergy(density, m_HintEl[iMode], d);
        o = transformedVibHamElement(m_HintVib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode]);
        energy[3] += d * o;
    }
    energy[0] = energy[1] + energy[2] + energy[3];
}


double OperatorBBO::transformedVibHamElement(Operator &hamiltonian, SMat &U, int r, int s) {
    //! @todo assert that there is only one symmetry block
    auto symm = U.symmetry();
    dim_t dim{8, 0};
    dim[symm] = 1;
    SMat Usplice({U.dimensions()[1], dim}, parityNone, symm);
    SMat UspliceT({dim, U.dimensions()[1]}, parityNone, symm);
    dim_t rowOffset(8, 0), colOffset(8, 0);
    colOffset[symm] = (size_t) s;
    Usplice.splice(U, rowOffset, colOffset);
    colOffset.assign(8, 0);
    rowOffset[symm] = (size_t) r;
    UspliceT.splice(U, rowOffset, colOffset);
    SMat el = (UspliceT * hamiltonian.O1(true) * Usplice);
    //! @todo check that this is a 1 by 1 matrix
    return el.block((unsigned) symm)[0];
}

Operator OperatorBBO::transformedVibHam(Operator &hamiltonian, SMat &U) {
    SMat Ut = U;
    Ut.transpose();
    Operator H(hamiltonian);
    H.O1(true) = Ut * hamiltonian.O1(true) * U;
    return H;
}

Operator OperatorBBO::electronicFock(Operator &P, std::vector<SMat> &U) {
    Operator F = m_Hel.fock(P, true, "Fock operator");
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        Operator Fint = m_HintEl[iMode].fock(P, true, "interaction component of Fock operator");
        double o = transformedVibHamElement(m_HintVib[iMode], U[iMode], m_vibOcc[iMode], m_vibOcc[iMode]);
        F.O1(true) += o * Fint.O1(true);
    }
    return F;
}

Operator OperatorBBO::vibrationalFock(Operator &P, SMat &U, int iMode) {
    dim_t dimension(8);
    dimension[m_symMode[iMode] - 1] = (unsigned int) m_nModal;
    std::vector<int> orbital_symmetries(m_nModal, m_symMode[iMode]);
    Operator F(dimension, orbital_symmetries, 1, false, (unsigned) m_symMode[iMode], true, true,
               "Vibrational Fock matrix, mode = " + std::to_string(iMode));
    F.O1(true) = transformedVibHam(m_HintVib[iMode], U).O1(true);
    double d;
    nm_RHF::electronicEnergy(P, m_HintEl[iMode], d);
    F.O1(true) *= 0.5 * d;
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
