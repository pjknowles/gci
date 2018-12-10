#include "gciMixedWavefunction.h"
#include "gciHProductSet.h"

#include <limits>
#include <utility>

namespace gci {

MixedWavefunction::MixedWavefunction(const Options &options)
        : m_vibSpace(options.parameter("NMODE", 0), options.parameter("NMODAL", 1),
                     options.parameter("VIB_EXC_LVL", 1)), m_vibBasis(m_vibSpace), m_elDim(0), m_dimension(0) {
    // Zeroth order
    auto state = State(options);
    m_wfn.emplace_back(state);
    m_elDim = m_wfn[0].size();
    m_vibBasis.generateFullSpace();
    m_wfn.resize(m_vibBasis.vibDim(), m_wfn[0]);
    m_dimension = m_vibBasis.vibDim() * m_elDim;
}

bool MixedWavefunction::empty() const {
    return m_wfn.empty();
}

void MixedWavefunction::allocate_buffer() {
    for (auto &el: m_wfn) el.allocate_buffer();
}

std::vector<size_t> MixedWavefunction::minlocN(size_t n) const {
    if (empty()) return {};
    // Find lowest n elements per vibrational basis. Store the absolute index and value
    std::vector<size_t> result;
    std::vector<std::pair<size_t, value_type>> minVals;
    for (int iVib = 0, offset = 0; iVib < m_vibBasis.vibDim(); ++iVib, offset += m_elDim) {
        int n_el = n > m_elDim ? m_elDim : n;
        result = m_wfn[iVib].minlocN(n_el);
        for (const auto &ind : result) {
            minVals.emplace_back(offset + ind, m_wfn[iVib].at(ind));
        }
    }
    // Search for the lowest n among all
    std::stable_sort(minVals.begin(), minVals.end(),
                     [](const auto &el1, const auto &el2) {return el1.second < el2.second;});
    result.resize(n);
    std::transform(minVals.begin(), std::next(minVals.begin(), n), result.begin(),
                   [](const auto &el) {return el.first;});
    return result;
}

size_t MixedWavefunction::minloc(size_t n) const {
    return minlocN(n).back();
}

double MixedWavefunction::at(size_t ind) const {
    if (ind >= m_dimension) throw std::logic_error("Out of bounds");
    size_t vibInd = ind / m_elDim;
    size_t elInd = ind - vibInd * m_elDim;
    return m_wfn[vibInd].at(elInd);
}

std::string MixedWavefunction::str() const {
    std::string out = "MixedWavefunction elements: ";
    for (size_t i = 0; i < m_dimension; ++i) {
        out += std::to_string(at(i)) + ", ";
    }
    return out;
}

void MixedWavefunction::operatorOnWavefunction(const MixedOperator &ham, const MixedWavefunction &w,
                                               bool parallel_stringset) {
    Wavefunction wfn_scaled(m_wfn[0]);
    for (const auto &bra : m_vibBasis) {
        auto iWfn = m_vibBasis.index(bra);
        auto expVal = ham.expectVal(bra, bra, {VibOpType::HO, {}});
        // Pure vibrational and electronic operators
        m_wfn[iWfn].axpy(expVal, w.m_wfn[iWfn]);
        m_wfn[iWfn].operatorOnWavefunction(ham.Hel, w.m_wfn[iWfn], parallel_stringset);
        // all mixed vibrational - electronic operators
        for (const auto &hamTerm : ham) {
            for (const auto &op : hamTerm.second) {
                auto connectedSet = HProductSet(bra, m_vibBasis.vibSpace(), op.vibOp);
                for (const auto &ket : connectedSet) {
                    auto jWfn = m_vibBasis.index(ket);
                    double expVal = ham.expectVal(bra, ket, op.vibOp);
                    if (std::abs(expVal) < 1.e-12) continue;
                    wfn_scaled = expVal * w.m_wfn[jWfn];
                    m_wfn[iWfn].operatorOnWavefunction(op.Hel, wfn_scaled, parallel_stringset);
                }
            }
        }
    }
}

void MixedWavefunction::diagonalOperator(const MixedOperator &ham, bool parallel_stringset) {
    Wavefunction wfn_scaled(m_wfn[0]);
    for (const auto &bra : m_vibBasis) {
        auto iWfn = m_vibBasis.index(bra);
        // Pure vibrational and electronic operators
        Wavefunction &wfn_diagonal = m_wfn[iWfn];
        wfn_diagonal.diagonalOperator(ham.Hel);
        auto w = ham.expectVal(bra, bra, {VibOpType::HO, {}});
        wfn_diagonal += w;
        // all mixed vibrational - electronic operators
        for (const auto &hamTerm : ham) {
            for (const auto &op : hamTerm.second) {
                auto expVal = ham.expectVal(bra, bra, op.vibOp);
                if (std::abs(expVal) < 1.e-12) continue;
                wfn_scaled = m_wfn[iWfn];
                wfn_scaled.diagonalOperator(op.Hel);
                wfn_scaled *= expVal;
                wfn_diagonal += wfn_scaled;
            }
        }
    }
}

bool MixedWavefunction::compatible(const MixedWavefunction &w2) const {
    bool sameSize = (m_wfn.size() == w2.m_wfn.size());
    if (!sameSize) return sameSize;
    bool sameVibBasis = (m_vibSpace == w2.m_vibSpace);
    bool sameElectronicWfn = true;
    for (int i = 0; i < m_wfn.size(); ++i) {
        sameElectronicWfn = sameElectronicWfn && m_wfn[i].compatible(w2.m_wfn[i]);
    }
    return sameSize && sameElectronicWfn && sameVibBasis;
}

void MixedWavefunction::axpy(double a, const MixedWavefunction &other) {
    for (int i = 0; i < m_wfn.size(); ++i) {
        m_wfn[i].axpy(a, other.m_wfn[i]);
    }
}

void MixedWavefunction::scal(double a) {
    for (auto &el: m_wfn) el.scal(a);
}

double MixedWavefunction::dot(const MixedWavefunction &other) const {
    return (*this) * ((dynamic_cast<const MixedWavefunction &>(other)));
}

void MixedWavefunction::zero() {
    for (auto &el : m_wfn) el.zero();
}

void MixedWavefunction::set(const double val) {
    for (auto &el:m_wfn) el.set(val);
}

void MixedWavefunction::set(const size_t ind, const double val) {
    if (ind >= m_dimension) throw std::logic_error("Out of bounds");
    size_t vibInd = ind / m_elDim;
    size_t elInd = ind - vibInd * m_elDim;
    m_wfn[vibInd].set(elInd, val);
}

MixedWavefunction &MixedWavefunction::operator*=(const double &value) {
    for (auto &el : m_wfn) el *= value;
    return *this;
}

MixedWavefunction &MixedWavefunction::operator+=(const MixedWavefunction &other) {
    if (!compatible(other)) throw std::domain_error("Attempting to add incompatible MixedWavefunction objects");
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
    if (!compatible(other)) throw std::domain_error("Attempting to divide incompatible MixedWavefunction objects");
    for (auto i = 0; i < m_wfn.size(); ++i) m_wfn[i] /= other.m_wfn[i];
    return *this;
}

MixedWavefunction &MixedWavefunction::addAbsPower(const MixedWavefunction &c, double k, double factor) {
    if (!compatible(c)) throw std::domain_error("MixedWavefunction::addAbsPower vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) m_wfn[i].addAbsPower(c.m_wfn[i], k, factor);
    return *this;
}

void MixedWavefunction::times(const MixedWavefunction *a, const MixedWavefunction *b) {
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    if (!compatible(aa) || !compatible(bb)) throw std::logic_error("Vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) m_wfn[i].times(&aa.m_wfn[i], &bb.m_wfn[i]);
}

void MixedWavefunction::divide(const MixedWavefunction *a, const MixedWavefunction *b,
                               double shift, bool append, bool negative) {
    if (a == nullptr || b == nullptr) throw std::logic_error("Vectors cannot be null");
    const auto aa = dynamic_cast<const MixedWavefunction &> (*a);
    const auto bb = dynamic_cast<const MixedWavefunction &> (*b);
    if (!compatible(aa) || !compatible(bb)) throw std::logic_error("Vectors must be compatible");
    for (int i = 0; i < m_wfn.size(); ++i) m_wfn[i].divide(&aa.m_wfn[i], &bb.m_wfn[i], shift, append, negative);
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
