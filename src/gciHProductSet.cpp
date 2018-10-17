#include "gciHProductSet.h"

#include <stdexcept>

namespace gci {

HProductSet::HProductSet(const VibSpace &vibSpace, const HProduct &bra, MixedOperator::VibOp vibOp)
        : m_vibSpace(vibSpace), m_basis(), m_vibDim(0), m_vibExcLvlDim() {
    switch (vibOp) {
        case MixedOperator::VibOp::HO: {
            m_basis.emplace_back(bra);
            break;
        }
        case MixedOperator::VibOp::Q:
        case MixedOperator::VibOp::dQ: {
        }
        case MixedOperator::VibOp::Qsq: {
        }
    }
    m_vibDim = m_basis.size();
}

void HProductSet::setVibDim() {
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
                    m_basis.emplace_back(prevProduct->excite(iMode));
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
    size_t index = 0;
    if (phi.empty()) index = 0;
    else if (phi.modeCouplingLvl() == 1) {
        auto iMode = phi[0][0];
        auto iModal = phi[0][1];
        index = m_vibExcLvlDim[0];
        index += (size_t) iMode * m_vibSpace.nModal + (iModal - 1);
    } else if (phi.modeCouplingLvl() == 2) {
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

}  // namespace gci
