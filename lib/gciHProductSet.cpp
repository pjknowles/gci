#include "gciHProductSet.h"

#include <stdexcept>
#include <numeric>

namespace gci {

HProductSet::HProductSet(const HProduct &bra, const VibSpace &vibSpace, const VibOp &vibOp)
        : m_vibSpace(vibSpace), m_basis(), m_vibDim(0), m_excLvlDim(vibSpace.excLvl + 1, 0),
          m_connectedSet(true) {
    switch (vibOp.type) {
        case VibOpType::HO:m_basis.emplace_back(bra);
            break;
        case VibOpType::Q:
        case VibOpType::dQ: generateQcoupledSpace(bra, vibOp);
            break;
        case VibOpType::Qsq: generateQsqCoupledSpace(bra, vibOp);
            break;
    }
    m_vibDim = m_basis.size();
}

void HProductSet::setVibDim() {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    // This can be generalized
    m_excLvlDim[0] = (size_t) 1;
    if (m_vibSpace.excLvl >= 1) {
        m_excLvlDim[1] = (size_t) m_vibSpace.nMode * (m_vibSpace.nModal - 1);
    }
    if (m_vibSpace.excLvl >= 2) {
        m_excLvlDim[2] = (size_t) m_excLvlDim[1] * (m_vibSpace.nMode - 1) / 2 * (m_vibSpace.nModal - 1);
    }
    m_vibDim = std::accumulate(m_excLvlDim.begin(), m_excLvlDim.end(), 0ul);
}

void HProductSet::generateFullSpace() {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    setVibDim();
    m_basis.reserve(m_vibDim);
    m_basis.emplace_back(HProduct{});
    unsigned long iWfn = 1;
    for (int exc = 1; exc < m_excLvlDim.size(); ++exc) {
        auto iStart = std::accumulate(m_excLvlDim.begin(), m_excLvlDim.begin() + exc - 1, 0ul);
        auto iEnd = std::accumulate(m_excLvlDim.begin(), m_excLvlDim.begin() + exc, 0ul);
        for (auto iBase = iStart; iBase < iEnd; ++iBase) {
            // Taking each product at previous level, apply single excitations to generate
            // all unique products (see HartreeProduct's ordering)
            auto base = m_basis[iBase];
            auto excitedModes = base.excitedModes();
            auto endMode = base.empty() ? m_vibSpace.nMode : excitedModes.front();
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
        if (iWfn != std::accumulate(m_excLvlDim.begin(), m_excLvlDim.begin() + exc + 1, 0ul))
            throw std::runtime_error("Generated inconsistent number of Hartree products");
    }
}

size_t HProductSet::index(const HProduct &phi) const {
    if (m_connectedSet) throw std::logic_error("Cannot be called for a connectedSet.");
    size_t index = 0;
    if (phi.empty()) index = 0;
    else if (phi.excLvl() == 1) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        index = m_excLvlDim[0];
        index += (size_t) iMode * (m_vibSpace.nModal - 1) + (iModal - 1);
    } else if (phi.excLvl() == 2) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        auto jMode = phi[1][0];
        auto jModal = phi[1][1];
        // previous exc levels
        index = (size_t) m_excLvlDim[0] + m_excLvlDim[1] - 1;
        // full subset of products with jMode -= 1
        index += jMode * (jMode - 1) / 2 * (m_vibSpace.nModal - 1) * (m_vibSpace.nModal - 1);
        // full subset of products with jModal -= 1
        index += (jModal - 1) * jMode * (m_vibSpace.nModal - 1);
        // the rest of subspace for the partial i (j is fixed)
        index += iMode * (m_vibSpace.nModal - 1) + iModal;
    } else throw std::logic_error("Index to more than doubly excited basis is not implemented yet.");
    return index;
}

namespace {
void check_and_emplace(std::vector<HProduct> &basis, HProduct &bra, const VibSpace &vibSpace) {
    if (bra.withinSpace(vibSpace)) {
        basis.emplace_back(bra);
    }
}
}

void HProductSet::generateQcoupledSpace(const HProduct &bra, const VibOp &vibOp) {
    if (vibOp.type != VibOpType::Q && vibOp.type != VibOpType::dQ) throw std::logic_error("Mismatch of operator type");
    // Raise and lower by one
    {
        auto newProd = bra;
        newProd.lower(vibOp.mode[0]);
        check_and_emplace(m_basis, newProd, m_vibSpace);
    }
    {
        auto newProd = bra;
        newProd.raise(vibOp.mode[0]);
        check_and_emplace(m_basis, newProd, m_vibSpace);
    }
}

void HProductSet::generateQsqCoupledSpace(const HProduct &bra, const VibOp &vibOp) {
    if (vibOp.type != VibOpType::Qsq) throw std::logic_error("Mismatch of operator type");
    if (vibOp.mode[0] == vibOp.mode[1]) {
        // Q_A * Q_A
        auto excLvl = bra.excLvl(vibOp.mode[0]);
        if (excLvl >= 2) {
            auto newProd = bra;
            newProd.changeModal(vibOp.mode[0], -2);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
        m_basis.push_back(bra);
        if (excLvl < m_vibSpace.nModal - 2) {
            auto newProd = bra;
            newProd.changeModal(vibOp.mode[0], +2);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
    } else if (m_vibSpace.excLvl >= 2) {
        // Q_A * Q_B
        // lower, lower
        {
            auto newProd = bra;
            newProd.lower(vibOp.mode[0]);
            newProd.lower(vibOp.mode[1]);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
        // lower, raise
        {
            auto newProd = bra;
            newProd.lower(vibOp.mode[0]);
            newProd.raise(vibOp.mode[1]);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
        // raise, lower
        {
            auto newProd = bra;
            newProd.raise(vibOp.mode[0]);
            newProd.lower(vibOp.mode[1]);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
        // raise, raise
        {
            auto newProd = bra;
            newProd.raise(vibOp.mode[0]);
            newProd.raise(vibOp.mode[1]);
            check_and_emplace(m_basis, newProd, m_vibSpace);
        }
    }
}

}  // namespace gci

































