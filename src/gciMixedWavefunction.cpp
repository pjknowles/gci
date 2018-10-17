#include "gciMixedWavefunction.h"
#include "gciHProductSet.h"

namespace gci {

MixedWavefunction::MixedWavefunction(const State &state, const VibSpace &vibSpace)
        : m_vibSpace(vibSpace), m_vibBasis(vibSpace), m_elDim(0), m_dimension(0) {
    // Zeroth order
    m_wfn.emplace_back(state);
    m_elDim = m_wfn[0].size();
    m_dimension = m_elDim;
    m_vibBasis.generateFullSpace();
    m_wfn.resize(m_vibBasis.vibDim(), m_wfn[0]);
    m_dimension += m_vibBasis.vibDim() * m_elDim;
}


std::vector<HProduct> MixedWavefunction::connectedVibBasis(const HProduct &phi) {

}

bool MixedWavefunction::compatible(const MixedWavefunction &w2) const {
    bool sameSize = (m_wfn.size() == w2.m_wfn.size());
    if (!sameSize) return sameSize;
    bool sameElectronicWfn = true;
    for (int i = 0; i < m_wfn.size(); ++i) {
        sameElectronicWfn = sameElectronicWfn && m_wfn[i].compatible(w2.m_wfn[i]);
    }
    bool sameVibBasis = (m_vibSpace == w2.m_vibSpace);
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
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    if (!compatible(aa) || !compatible(bb)) throw std::logic_error("Vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].times(&aa.m_wfn[i], &bb.m_wfn[i]);
    }
}

void MixedWavefunction::divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b,
                               double shift, bool append, bool negative) {
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    if (!compatible(aa) || !compatible(bb)) throw std::logic_error("Vectors must be compatible");
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
    for (const auto &bra : m_vibBasis) {
        auto iWfn = indexVibWfn(bra);
        // Hvib
        m_wfn[iWfn].axpy(O_Hvib(hamiltonian, bra), w.m_wfn[iWfn]);
        // Hel
        m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel, w.m_wfn[iWfn], parallel_stringset);
        // Q and dQ operators connect CIS and CISD basis
        // Hint * Q[iMode]
        for (const auto &ket : connectedVibBasis(bra, vibOp::Q, iMode)) {
            auto jWfn = indexVibWfn(ket);
            // H'el(A)
            wfn_scaled = O_Q(hamiltonian, bra, ket) * m_wfn[jWfn];
            m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel_A[iMode], wfn_scaled, parallel_stringset);
            // Hs(A)
            wfn_scaled = O_dQ(hamiltonian, bra, ket) * m_wfn[jWfn];
            m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hs_A[iMode], wfn_scaled, parallel_stringset);
        }
    }

// Vibrational CIS
    Wavefunction wfn_scaled(m_wfn[0]);
    for (auto iMode = 0; iMode < m_nMode; ++iMode) {
        for (auto iModal = 0; iModal < m_nModal; ++iModal) {
            auto bra = HartreeProduct({{iMode, iModal}});
            auto iWfn = indexVibWfn(bra);
            // Hvib
            m_wfn[iWfn].axpy(O_Hvib(hamiltonian, bra), w.m_wfn[iWfn]);
            // Hel
            m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel, w.m_wfn[iWfn], parallel_stringset);
            // Q and dQ operators connect CIS and CISD basis
            // Hint * Q[iMode]
            for (const auto &ket : connectedVibBasis(bra, vibOp::Q, iMode)) {
                auto jWfn = indexVibWfn(ket);
                // H'el(A)
                wfn_scaled = O_Q(hamiltonian, bra, ket) * m_wfn[jWfn];
                m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel_A[iMode], wfn_scaled, parallel_stringset);
                // Hs(A)
                wfn_scaled = O_dQ(hamiltonian, bra, ket) * m_wfn[jWfn];
                m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hs_A[iMode], wfn_scaled, parallel_stringset);
            }
        }
    }
// Vibrational CID component
    for (auto iMode = 0; iMode < m_nMode; ++iMode) {
        for (auto iModal = 1; iModal < m_nModal; ++iModal) {
            for (auto jMode = iMode + 1; jMode < m_nMode; ++jMode) {
                for (auto jModal = 1; jModal < m_nModal; ++jModal) {
                    auto bra = HartreeProduct({{iMode, iModal}, {jMode, jModal}});
                    auto iWfn = indexVibWfn(bra);
                    // Hvib
                    m_wfn[iWfn].axpy(O_Hvib(hamiltonian, bra), w.m_wfn[iWfn]);
                    // Hel
                    m_wfn[iWfn].operatorOnWavefunction(hamiltonian.Hel, w.m_wfn[iWfn], parallel_stringset);
                    // Hint * Q[iMode]
                    for (const auto ket : connectedVibBasis(bra, "Q", 0)) {
                    }
                    // Hint * Q[1]
                    for (const auto ket : connectedVibBasis(bra, "Q", 1)) {
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
            m_wfn[iWfn] += O_Hvib(hamiltonian, iMode);
        }
    }

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
