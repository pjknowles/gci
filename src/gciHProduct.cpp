#include "gciHProduct.h"

#include <algorithm>
#include <stdexcept>

namespace gci {
HProduct::HProduct(t_Product phi) : m_prod(std::move(phi)) {
    reorder(m_prod);
    check(m_prod);
}

int HProduct::excLvl(int iMode) const {
    auto index = std::find_if(m_prod.cbegin(), m_prod.cend(), [iMode](auto modal) {return modal[0] == iMode;});
    if (index == m_prod.end()) return 0;
    else return (*index)[1];
}

void HProduct::reorder(t_Product &prod) const {
    // remove any modalIndex == 0
    auto newEnd = std::remove_if(prod.begin(), prod.end(), [](t_Modal &modal) {return modal[1] == 0;});
    prod.erase(newEnd, prod.end());
    // sort by modeIndex
    std::sort(prod.begin(), prod.end(), [](t_Modal &modalA, t_Modal &modalB) {return modalA[0] < modalB[0];});
}

void HProduct::check(t_Product &prod) const {
    if (std::any_of(m_prod.begin(), m_prod.end(), [](t_Modal &modal) {return modal[0] < 0 || modal[1] < 0;}))
        throw std::logic_error("Product indices must be positive integers.");
    for (auto modal = m_prod.begin(); modal != m_prod.end() - 1; ++modal) {
        if (modal[0][0] == modal[1][0])
            throw std::logic_error("HProduct with more than a single modal occupied per mode.");
    }
}

std::vector<int> HProduct::excitedModes() const {
    std::vector<int> modes(excLvl());
    std::transform(m_prod.begin(), m_prod.end(), modes.begin(), [](const t_Modal &in) {return in[0];});
    return modes;
}

bool HProduct::withinSpace(const VibSpace &vibSpace) {
    return !std::any_of(m_prod.begin(), m_prod.end(),
                        [&vibSpace](const t_Modal &modal) {
                            return modal[0] >= vibSpace.nMode || modal[1] >= vibSpace.nModal;
                        });
}

void HProduct::changeModal(const int iMode, const int diff) {
    auto modal = std::find_if(m_prod.begin(), m_prod.end(),
                              [&iMode](t_Modal &el) {return (el[0] == iMode);});
    if (modal == m_prod.end()) {
        m_prod.emplace_back(t_Modal({iMode, diff}));
        reorder(m_prod);
    } else {
        (*modal)[1] += diff;
        check(m_prod);
    }
}

}  // namespace gci
