#include "gciMixedWavefunction.h"

#include <assert.h>

namespace gci {
MixedWavefunction::MixedWavefunction(const State &state, int nMode, int nModal, int modeCoupling)
        : m_nMode(nMode), m_nModal(nModal), m_modeCoupling(modeCoupling), m_elDim(0), m_vibDim(0),
          m_vibExcLvlDim(modeCoupling + 1, 0), m_dimension(0) {
    // Zeroth order
    m_wfn.emplace_back(state);
    m_elDim = m_wfn[0].size();
    m_dimension = m_elDim;
    // Vibrationally coupled states
    if (m_modeCoupling > 2) {
        throw std::logic_error("modeCoupling > 2 is not supported yet");
    }
    setVibDim();
    m_wfn.resize(m_vibDim, m_wfn[0]);
    m_dimension += m_vibDim * m_elDim;
}

void MixedWavefunction::setVibDim() {
    // This can be generalized
    m_vibExcLvlDim[0] = (size_t) 0;
    if (m_modeCoupling >= 1) {
        m_vibExcLvlDim[1] = (size_t) m_nMode * (m_nModal - 1);
    }
    if (m_modeCoupling >= 2) {
        m_vibExcLvlDim[2] = (size_t) m_vibExcLvlDim[1] * (m_nMode - 1) / 2 * (m_nModal - 1);
    }
    m_vibDim = std::accumulate(m_vibExcLvlDim.cbegin(), m_vibExcLvlDim.cend(), (size_t) 0);
}

// TODO Implement unit testing
size_t MixedWavefunction::indexVibWfn(HartreeProduct phi) {
    if (phi.empty() == 0) return 0;
    else if (phi.size() == 1) {
        if (phi[0][0] == 0) return indexVibWfn({}); // Excitation is a ground state
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        return (size_t) iMode * (m_nModal - 1) + iModal;
    } else if (phi.size() == 2) {
        if (phi[0][0] == 0) return indexVibWfn({phi[1]}); // First excitation is a ground state
        if (phi[1][0] == 0) return indexVibWfn({phi[0]}); // Second excitation is a ground state
        if (phi[0][0] == phi[1][0]) throw std::logic_error("Double excitation of the same mode is not allowed");
        if (phi[0][0] < phi[1][0]) phi[0].swap(phi[1]);
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        auto jMode = phi[1][0];
        auto jModal = phi[1][1];
        auto index = (size_t) m_vibExcLvlDim[1];
        index += iMode * (iMode + 1) / 2 * (m_nModal - 1) * (m_nModal - 1);
        index += (iModal - 1) * m_nModal + (jModal - 1);
        return index;
    } else throw std::logic_error("Vibrational basis is currently restricted to double excitations only");
}

bool MixedWavefunction::compatible(const MixedWavefunction &w2) const {
    bool sameSize = (m_wfn.size() == w2.m_wfn.size());
    if (!sameSize) return sameSize;
    bool sameElectronicWfn = true;
    for (int i = 0; i < m_wfn.size(); ++i) {
        sameElectronicWfn = sameElectronicWfn && m_wfn[i].compatible(w2.m_wfn[i]);
    }
    bool sameVibBasis = (m_nMode == w2.m_nMode) && (m_nModal == w2.m_nModal) && (m_modeCoupling == w2.m_modeCoupling);
    return sameSize && sameElectronicWfn && sameVibBasis;
}

void MixedWavefunction::allocate_buffer() {
    for (auto &el: m_wfn) el.allocate_buffer();
}

void MixedWavefunction::axpy(double a, const LinearAlgebra::vector<double> &other) {
    const auto &x = dynamic_cast <const MixedWavefunction &> (other);
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].axpy(a, x.m_wfn[i]);
    }
}

void MixedWavefunction::scal(double a) {
    for (auto &el: m_wfn) el.scal(a);
}

double MixedWavefunction::dot(const LinearAlgebra::vector<double> &other) const {
    return (*this) * ((dynamic_cast<const MixedWavefunction &>(other)));
}

void MixedWavefunction::zero() {
    for (auto &el : m_wfn) el.zero();
}

MixedWavefunction &MixedWavefunction::operator*=(const double &value) {
    for (auto &el : m_wfn) el *= value;
    return *this;
}

MixedWavefunction &MixedWavefunction::operator+=(const MixedWavefunction &other) {
    if (!compatible(other))
        throw std::domain_error("Attempting to add incompatible MixedWavefunction objects");
    for (auto i = 0; i < m_wfn.size(); ++i) m_wfn[i] += other.m_wfn[i];
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-=(const MixedWavefunction &other) {
    if (!compatible(other))
        throw std::domain_error("Attempting to subtract incompatible MixedWavefunction objects");
    for (auto i = 0; i < m_wfn.size(); ++i) m_wfn[i] -= other.m_wfn[i];
    return *this;
}

MixedWavefunction &MixedWavefunction::operator+=(double value) {
    for (auto &el: m_wfn) el += value;
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-=(double value) {
    for (auto &el: m_wfn) el -= value;
    return *this;
}

MixedWavefunction &MixedWavefunction::operator-() {
    for (auto &el: m_wfn) el = -el;
    return *this;
}

MixedWavefunction &MixedWavefunction::operator/=(const MixedWavefunction &other) {
    if (!compatible(other))
        throw std::domain_error("Attempting to divide incompatible MixedWavefunction objects");
    for (auto i = 0; i < m_wfn.size(); ++i) m_wfn[i] /= other.m_wfn[i];
    return *this;
}

MixedWavefunction &MixedWavefunction::addAbsPower(const MixedWavefunction &c, double k, double factor) {
    if (!compatible(c)) throw std::domain_error("MixedWavefunction::addAbsPower vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].addAbsPower(c.m_wfn[i], k, factor);
    }
    return *this;
}

void MixedWavefunction::times(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b) {
    assert(a != nullptr && b != nullptr && "Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    assert(compatible(aa) && compatible(bb) && "All vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].times(&aa.m_wfn[i], &bb.m_wfn[i]);
    }
}

void MixedWavefunction::divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b,
                               double shift, bool append, bool negative) {
    assert(a != nullptr && b != nullptr && "Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    assert(compatible(aa) && compatible(bb) && "All vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].divide(&aa.m_wfn[i], &bb.m_wfn[i], shift, append, negative);
    }
}

double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    if (!w1.compatible(w2))
        throw std::domain_error("Attempt to form scalar product between incompatible MixedWavefunction objects");
    double result = 0.0;
    for (int i = 0; i < w1.size(); ++i) {
        result += w1.m_wfn[i] * w2.m_wfn[i];
    }
    return result;
}

void MixedWavefunction::operatorOnWavefunction(const MixedOperator &hamiltonian, const MixedWavefunction &w,
                                               bool parallel_stringset) {
    // Vibrational ground state
    m_wfn[0].operatorOnWavefunction(hamiltonian.Hel, w.m_wfn[0], parallel_stringset);
    // Vibrational CIS
    Wavefunction wfn_scaled(m_wfn[0]);
    for (auto iMode = 0; iMode < m_nMode; ++iMode) {
        for (auto iModal = 0; iModal < m_nModal; ++iModal) {
            auto iWfn = indexVibWfn({{iMode, iModal}});
            // Hvib
            m_wfn[iWfn].axpy(O_Hvib(hamiltonian, iMode, iModal), w.m_wfn[iWfn]);
            // Hel
            m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel, w.m_wfn[iWfn], parallel_stringset);
            1
            for (int diff = -1; diff < 2; diff += 2) {
                auto jModal = iModal + diff;
                if (jModal < 0 || jModal >= m_nModal) continue;
                auto jWfn = indexVibWfn({{iMode, jModal}});
                // H'el(A)
                wfn_scaled = O_Q(hamiltonian, iMode, iModal, jModal) * m_wfn[jWfn];
                m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel_A[iMode], wfn_scaled, parallel_stringset);
                // Hs(A)
                wfn_scaled = O_dQ(hamiltonian, iMode, iModal, jModal) * m_wfn[jWfn];
                m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hs_A[iMode], wfn_scaled, parallel_stringset);
            }
        }
    }
    // TODO implement Hel(A,B) matrix elements.
    for (auto iMode = 0; iMode < m_nMode; ++iMode) {
        for (auto iModal = 0; iModal < m_nModal; ++iModal) {
            for (auto jMode = 0; jMode < m_nMode; ++jMode) {
                for (auto jModal = 0; jModal < m_nModal; ++jModal) {
                    for (const auto jWfn : connectedVibBasis({{iMode, iModal}, {jMode, jModal}})) {

                    }
                }
            }
        }
    }
}

void MixedWavefunction::diagonalOperator(const MixedOperator &hamiltonian) {
    m_wfn[0].diagonalOperator(hamiltonian.Hel);
    for (auto iMode = 0, iWfn = 0; iMode < m_nMode; ++iMode) {
        for (auto iModal = 0; iModal < m_nModal; ++iModal, ++iWfn) {
            if (iWfn == 0) continue;
            // Hel
            m_wfn[iWfn] = m_wfn[0];
            // Hvib
            m_wfn[iWfn] += O_Hvib(hamiltonian, iMode, iModal);
        }
    }

}

double MixedWavefunction::O_Hvib(const MixedOperator &hamiltonian, int mode, int modal) {
    return 0.5 * hamiltonian.freq[mode] * modal;
}

double MixedWavefunction::O_Q(const MixedOperator &hamiltonian, int mode, int iModal, int jModal) {
    assert(std::abs(iModal - jModal) == 1 && "Operator out of range");
    double n = iModal > jModal ? iModal : jModal;
    return std::sqrt(n / (2.0 * hamiltonian.freq[mode]));
}

double MixedWavefunction::O_dQ(const MixedOperator &hamiltonian, int mode, int iModal, int jModal) {
    assert(std::abs(iModal - jModal) == 1 && "Operator out of range");
    double n = iModal > jModal ? iModal : jModal;
    int sign = iModal > jModal ? -1 : 1;
    return sign * std::sqrt(0.5 * hamiltonian.freq[mode] * n);
}

void MixedWavefunction::set(const double val) {
    for (auto &el:m_wfn) el.set(val);
}

size_t MixedWavefunction::minloc(size_t n) const {
    return 0;
}

double MixedWavefunction::at(size_t offset) const {
    return 0;
}

MixedWavefunction operator+(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result += w2;
}

MixedWavefunction operator-(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result -= w2;
}

MixedWavefunction operator/(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    MixedWavefunction result = w1;
    return result /= w2;
}

MixedWavefunction operator*(const MixedWavefunction &w1, const double &value) {
    MixedWavefunction result = w1;
    return result *= value;
}

MixedWavefunction operator*(const double &value, const MixedWavefunction &w1) {
    MixedWavefunction result = w1;
    return result *= value;
}

}  // namespace gci
