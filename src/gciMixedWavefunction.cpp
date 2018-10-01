#include "gciMixedWavefunction.h"

namespace gci {
MixedWavefunction::MixedWavefunction(const State &state, int nMode, int nModal, int modeCoupling)
        : m_nMode(nMode), m_nModal(nModal), m_modeCoupling(modeCoupling), m_elecWfnSize(0), m_dimension(0) {
    assert((modeCoupling == 1 || modeCoupling == 2) && "Current Hamiltonian assumes 1 or 2 mode coupling");
    // Zeroth order
    m_wfn.emplace_back(state);
    m_elecWfnSize = m_wfn[0].size();
    m_dimension = m_elecWfnSize;
    // First order coupling
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        for (int iModal = 0; iModal < m_nModal; ++iModal) {
            m_wfn.emplace_back(state);
        }
    }
    m_dimension += m_nMode * m_nModal * m_elecWfnSize;
    // Second order coupling
    if (m_modeCoupling > 1) {
        throw std::logic_error("modeCoupling > 1 is not yet implemented");
    }
}

void MixedWavefunction::operatorOnWavefunction(const gci::Operator &h, const Wavefunction &w, bool parallel_stringset) {

}

bool MixedWavefunction::compatible(const MixedWavefunction &w2) const {
    bool b_size = (m_wfn.size() == w2.m_wfn.size());
    if (!b_size) return b_size;
    bool b_electronic_wfn = true;
    for (int i = 0; i < m_wfn.size(); ++i) {
        b_electronic_wfn = b_electronic_wfn && m_wfn[i].compatible(w2.m_wfn[i]);
    }
    bool b_vib_basis = (m_nMode == w2.m_nMode) && (m_nModal == w2.m_nModal) && (m_modeCoupling == w2.m_modeCoupling);
    return b_size && b_electronic_wfn && b_vib_basis;
}

void MixedWavefunction::allocate_buffer() {
    for(auto & el: m_wfn) el.allocate_buffer();
}

void MixedWavefunction::axpy(double a, const LinearAlgebra::vector<double> &other) {
    const auto &x = dynamic_cast <const MixedWavefunction &> (other);
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].axpy(a, x.m_wfn[i]);
    }
}

//std::tuple<std::vector<size_t>, std::vector<double> >
//MixedWavefunction::select(const LinearAlgebra::vector<double> &measure,
//                          const size_t maximumNumber,
//                          const double threshold) const {
//
//}

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
        m_wfn[i].times(aa.m_wfn[i], bb.m_wfn[i]);
    }
}

void MixedWavefunction::divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b,
                               double shift, bool append, bool negative) {
    assert(a != nullptr && b != nullptr && "Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    assert(compatible(aa) && compatible(bb) && "All vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].divide(aa.m_wfn[i], bb.m_wfn[i], shift, append, negative);
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
