#include "gciHProductSet.h"

#include <stdexcept>

namespace gci {

HProductSet::HProductSet(const VibSpace &vibSpace, const HProduct &bra, const VibOp &vibOp)
        : m_vibSpace(vibSpace), m_basis(), m_vibDim(0), m_vibExcLvlDim(), m_connectedSet(true) {
    switch (vibOp.type) {
        case VibOpType::HO: {
            m_basis.emplace_back(bra);
            break;
        }
        case VibOpType::Q:
        case VibOpType::dQ: generateQcoupledSpace(bra, vibOp);
        case VibOpType::Qsq: generateQsqCoupledSpace(bra, vibOp);
    }
    m_vibDim = m_basis.size();
}

void HProductSet::setVibDim() {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    // This can be generalized
    auto vibExcLvlDim = m_vibExcLvlDim;
    vibExcLvlDim[0] = (size_t) 1;
    if (m_vibSpace.modeCoupling >= 1) {
        vibExcLvlDim[1] = (size_t) m_vibSpace.nMode * (m_vibSpace.nModal - 1);
    }
    if (m_vibSpace.modeCoupling >= 2) {
        vibExcLvlDim[2] = (size_t) vibExcLvlDim[1] * (m_vibSpace.nMode - 1) / 2 * (m_vibSpace.nModal - 1);
    }
    m_vibDim = vibExcLvlDim[0];
    m_vibExcLvlDim[0] = vibExcLvlDim[0];
    for (int iExc = 1; iExc < vibExcLvlDim.size(); ++iExc) {
        m_vibDim += vibExcLvlDim[iExc];
        m_vibExcLvlDim[iExc] = m_vibDim;
    }
}

void HProductSet::generateFullSpace() {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    setVibDim();
    m_basis.reserve(m_vibDim);
    m_basis.emplace_back(HProduct::t_Modal{{}});
    unsigned long iWfn = 1;
    for (int exc = 1; exc < m_vibExcLvlDim.size() - 1; ++exc) {
        for (auto iBase = m_vibExcLvlDim[exc - 1]; iBase < m_vibExcLvlDim[exc]; ++iBase) {
// Taking each product at previous level, apply single excitations to generate unique products (see HartreeProducts ordering)
            auto base = m_basis[iBase];
            auto excitedModes = base.excitedModes();
            auto endMode = base.empty() ? m_vibSpace.nMode : excitedModes.front() + 1;
            for (int iMode = 0; iMode < endMode; ++iMode) {
                auto *prevProduct = &base;
                for (int iModal = 1; iModal < m_vibSpace.nModal; ++iModal) {
                    auto newProd = HProduct(*prevProduct);
                    newProd.raise(iMode);
                    m_basis.emplace_back(newProd);
                    prevProduct = &(m_basis.back());
                    ++iWfn;
                }
            }
        }
        if (iWfn != m_vibExcLvlDim[exc])
            throw std::runtime_error("Generated inconsistent number of Hartree products");
    }
}

size_t HProductSet::index(const HProduct &phi) {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    size_t index = 0;
    if (phi.empty()) index = 0;
    else if (phi.excLvl() == 1) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        index = m_vibExcLvlDim[0];
        index += (size_t) iMode * m_vibSpace.nModal + (iModal - 1);
    } else if (phi.excLvl() == 2) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        auto jMode = phi[1][0];
        auto jModal = phi[1][1];
        index = (size_t) m_vibExcLvlDim[1];
        index += iMode * (iMode + 1) / 2 * (m_vibSpace.nModal - 1) * (m_vibSpace.nModal - 1);
        index += (iModal - 1) * (m_vibSpace.nModal - 1) + (jModal - 1);
    } else throw std::logic_error("Index to more than doubly excited basis is not implemented yet.");
    return index;
}

void HProductSet::generateQcoupledSpace(const HProduct &bra, const VibOp &vibOp) {
    if (vibOp.type != VibOpType::Q && vibOp.type != VibOpType::dQ) throw std::logic_error("Mismatch of operator type");
    // Rase and lower by one
    auto excLvl = bra.excLvl(vibOp.mode[0]);
    if (excLvl >= 1) {
        auto newProd = bra;
        newProd.lower(vibOp.mode[0]);
        m_basis.emplace_back(newProd);
    }
    if (excLvl < m_vibSpace.nModal - 1) {
        auto newProd = bra;
        newProd.raise(vibOp.mode[0]);
        m_basis.emplace_back(newProd);
    }
}

void HProductSet::generateQsqCoupledSpace(const HProduct &bra, const VibOp &vibOp) {
    if (vibOp.type != VibOpType::Qsq) throw std::logic_error("Mismatch of operator type");
    if (vibOp.mode.size() == 1) {
        // Q_A * Q_A
        auto excLvl = bra.excLvl(vibOp.mode[0]);
        if (excLvl >= 2){
            auto newProd = bra;
            newProd.changeModal(vibOp.mode[0],-2);
            m_basis.emplace_back(newProd);
        }
        m_basis.push_back(bra);
        if (excLvl < m_vibSpace.nModal - 2) {
            auto newProd = bra;
            newProd.changeModal(vibOp.mode[0],+2);
            m_basis.emplace_back(newProd);
        }
    } else {
        // Q_A * Q_B
        auto excLvlA = bra.excLvl(vibOp.mode[0]);
        auto excLvlB = bra.excLvl(vibOp.mode[1]);
        // lower, lower
        if (excLvlA >= 1 && excLvlB >= 1){
            auto newProd = bra;
            newProd.lower(vibOp.mode[0]);
            newProd.lower(vibOp.mode[1]);
            m_basis.emplace_back(newProd);
        }
        // raise, lower
        if (excLvlA < m_vibSpace.nModal - 1 && excLvlB >= 1){
            auto newProd = bra;
            newProd.raise(vibOp.mode[0]);
            newProd.lower(vibOp.mode[1]);
            m_basis.emplace_back(newProd);
        }
        // lower, raise
        if (excLvlA >= 1 && excLvlB < m_vibSpace.nModal - 1){
            auto newProd = bra;
            newProd.lower(vibOp.mode[0]);
            newProd.raise(vibOp.mode[1]);
            m_basis.emplace_back(newProd);
        }
        // raise, raise
        if (excLvlA < m_vibSpace.nModal - 1 && excLvlB < m_vibSpace.nModal - 1){
            auto newProd = bra;
            newProd.raise(vibOp.mode[0]);
            newProd.raise(vibOp.mode[1]);
            m_basis.emplace_back(newProd);
        }
    }
}

}  // namespace gci

































