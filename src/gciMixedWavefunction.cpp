#include "gciMixedWavefunction.h"

#include <algorithm>

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
    generateVibrationalSpace();
    m_wfn.resize(m_vibDim, m_wfn[0]);
    m_dimension += m_vibDim * m_elDim;
}

void MixedWavefunction::setVibDim() {
    // This can be generalized
    auto vibExcLvlDim = m_vibExcLvlDim;
    vibExcLvlDim[0] = (size_t) 1;
    if (m_modeCoupling >= 1) {
        vibExcLvlDim[1] = (size_t) m_nMode * (m_nModal - 1);
    }
    if (m_modeCoupling >= 2) {
        vibExcLvlDim[2] = (size_t) vibExcLvlDim[1] * (m_nMode - 1) / 2 * (m_nModal - 1);
    }
    m_vibDim = vibExcLvlDim[0];
    m_vibExcLvlDim[0] = vibExcLvlDim[0];
    for (int iExc = 1; iExc < vibExcLvlDim.size(); ++iExc) {
        m_vibDim += vibExcLvlDim[iExc];
        m_vibExcLvlDim[iExc] = m_vibDim;
    }
}

void MixedWavefunction::generateVibrationalSpace() {
    setVibDim();
    m_vibBasis.reserve(m_vibDim);
    m_vibBasis.emplace_back(std::initializer_list<HartreeProduct::t_Modal>{{}});
    size_t iWfn = 1;
    for (int exc = 1; exc < m_vibExcLvlDim.size() - 1; ++exc) {
        for (auto iBase = m_vibExcLvlDim[exc - 1]; iBase < m_vibExcLvlDim[exc]; ++iBase) {
            // Taking each product at previous level, apply single excitations on unique modes
            auto base = m_vibBasis[iBase];
            auto excitedModes = base.excitedModes();
            int startMode = base.empty() ? 0 : excitedModes.front() + 1;
            for (int iMode = startMode; iMode < m_nMode ; ++iMode){
                auto* prevProduct = &base;
                for(int iModal = 1; iModal < m_nModal; ++iModal){
                    m_vibBasis[iWfn] = prevProduct->excite(iMode);
                    prevProduct = &(m_vibBasis[iWfn]);
                    ++iWfn;
                }
            }
        }
        if (iWfn != m_vibExcLvlDim[exc]) throw std::runtime_error("Generated inconsistent number of Hartree products");
    }
}

MixedWavefunction::HartreeProduct::HartreeProduct(t_Product phi) : m_prod(std::move(phi)) {
    reorder(m_prod);
}

void MixedWavefunction::HartreeProduct::reorder(t_Product &prod) {
    // remove any modalIndex == 0
    auto newEnd = std::remove_if(prod.begin(), prod.end(), [](t_Modal &a) {return a[1] == 0;});
    prod.erase(newEnd, prod.end());
    // sort by modeIndex
    std::sort(prod.begin(), prod.end(), [](t_Modal &a, t_Modal &b) {return a[0] < b[0];});
    // if there are any duplicates of modeIndex, throw an error
    for (auto modal = m_prod.begin(); modal < m_prod.end() - 1; ++modal) {
        if (modal[0][0] == modal[1][0])
            throw std::logic_error("HartreeProduct with more than a single modal occupied per mode");
    }
}

MixedWavefunction::HartreeProduct MixedWavefunction::HartreeProduct::excite(int iMode) const {
//    if (empty()) return t_Product({{iMode, 1}});
    auto excitedProduct = t_Product(m_prod);
    auto modal = std::find_if(excitedProduct.begin(), excitedProduct.end(),
                              [&iMode](t_Modal &el) {return (el[0] == iMode);});
    if (modal == m_prod.end()) {
        excitedProduct.push_back({{iMode, 1}});
    } else {
        (*modal)[1] += 1;
    }
    return excitedProduct;
}

std::vector<int> MixedWavefunction::HartreeProduct::excitedModes() const {
    std::vector<int> modes(size());
    std::transform(m_prod.begin(), m_prod.end(), modes.begin(), [](const t_Modal &in) {return in[0];});
    return modes;
}

// TODO Implement unit testing
size_t MixedWavefunction::indexVibWfn(const HartreeProduct &phi) {
    if (phi.empty() == 0) return 0;
    else if (phi.size() == 1) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        return (size_t) iMode * (m_nModal - 1) + iModal;
    } else if (phi.size() == 2) {
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

std::vector<MixedWavefunction::HartreeProduct> MixedWavefunction::connectedVibBasis(const HartreeProduct &phi) {

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

}
// Vibrational CIS
Wavefunction wfn_scaled(m_wfn[0]);
for (
auto iMode = 0;
iMode<m_nMode;
++iMode) {
for (
auto iModal = 0;
iModal<m_nModal;
++iModal) {
auto bra = HartreeProduct({{iMode, iModal}});
auto iWfn = indexVibWfn(bra);
// Hvib
m_wfn[iWfn].
axpy(O_Hvib(hamiltonian, bra), w
.m_wfn[iWfn]);
// Hel
m_wfn[iWfn].
operatorOnWavefunction(hamiltonian
.Hel, w.m_wfn[iWfn], parallel_stringset);
// Q and dQ operators connect CIS and CISD basis
// Hint * Q[iMode]
for (
const auto &ket
:
connectedVibBasis(bra, vibOp::Q, iMode
)) {
auto jWfn = indexVibWfn(ket);
// H'el(A)
wfn_scaled = O_Q(hamiltonian, bra, ket) * m_wfn[jWfn];
m_wfn[iWfn].
operatorOnWavefunction(hamiltonian
.Hel_A[iMode], wfn_scaled, parallel_stringset);
// Hs(A)
wfn_scaled = O_dQ(hamiltonian, bra, ket) * m_wfn[jWfn];
m_wfn[iWfn].
operatorOnWavefunction(hamiltonian
.Hs_A[iMode], wfn_scaled, parallel_stringset);
}
}
}
// Vibrational CID component
for (
auto iMode = 0;
iMode<m_nMode;
++iMode) {
for (
auto iModal = 1;
iModal<m_nModal;
++iModal) {
for (
auto jMode = iMode + 1;
jMode<m_nMode;
++jMode) {
for (
auto jModal = 1;
jModal<m_nModal;
++jModal) {
auto bra = HartreeProduct({{iMode, iModal}, {jMode, jModal}});
auto iWfn = indexVibWfn(bra);
// Hvib
m_wfn[iWfn].
axpy(O_Hvib(hamiltonian, bra), w
.m_wfn[iWfn]);
// Hel
m_wfn[iWfn].
operatorOnWavefunction(hamiltonian
.Hel, w.m_wfn[iWfn], parallel_stringset);
// Hint * Q[iMode]
for (
const auto ket
:
connectedVibBasis(bra,
"Q", 0)) {
}
// Hint * Q[1]
for (
const auto ket
:
connectedVibBasis(bra,
"Q", 1)) {
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

double MixedWavefunction::O_Hvib(const MixedOperator &hamiltonian, const HartreeProduct &bra) {
    double O_Hvib = 0.0;
    for (const auto &el :bra) {
        auto iMode = el[0];
        auto iModal = el[1];
        O_Hvib += 0.5 * hamiltonian.freq[iMode] * iModal;
    }
    return O_Hvib;
}

double MixedWavefunction::O_Q(const MixedOperator &hamiltonian, const HartreeProduct &bra, const HartreeProduct &ket) {
    if (bra.size() != ket.size()) throw std::logic_error("Bra and Ket are of different mode-coupling level. Always 0.");
    double O_Q = 0.0;
    auto nExc = bra.size();
    for (int i = 0; i < nExc; ++i) {
        if (bra[i][0] != ket[i][0])
            throw std::logic_error("Bra and Ket are of different mode-coupling level. Always 0.");
        auto mode = bra[i][0], iModal = bra[i][1], jModal = ket[i][1];
        auto diff = iModal - jModal;
        if (diff == 0) continue;
        if (std::abs(diff) != 1)
            throw std::logic_error("Bra and Ket are separated by more than 1 excitation. Always 0.");
        double n = diff > 0 ? iModal : jModal;
        O_Q += std::sqrt(n / (2.0 * hamiltonian.freq[mode]));
    }
    return O_Q;
}

double MixedWavefunction::O_dQ(const MixedOperator &hamiltonian, const HartreeProduct &bra, const HartreeProduct &ket) {
    if (std::abs(iModal - jModal) != 1) throw std::logic_error("Operator out of range");
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
