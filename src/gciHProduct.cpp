#include "gciHProduct.h"

#include <algorithm>
#include <stdexcept>

HProduct::HProduct(t_Product phi) : m_prod(std::move(phi)) {
    reorder(m_prod);
}


void HProduct::reorder(t_Product &prod) {
    // remove any modalIndex == 0
    auto newEnd = std::remove_if(prod.begin(), prod.end(), [](t_Modal &a) {return a[1] == 0;});
    prod.erase(newEnd, prod.end());
    // sort by modeIndex
    std::sort(prod.begin(), prod.end(), [](t_Modal &a, t_Modal &b) {return a[0] < b[0];});
    // if there are any duplicates of modeIndex, throw an error
    for (auto modal = m_prod.begin(); modal < m_prod.end() - 1; ++modal) {
        if (modal[0][0] == modal[1][0])
            throw std::logic_error("HProduct with more than a single modal occupied per mode");
    }
}

HProduct HProduct::excite(const int iMode) const {
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


std::vector<int> HProduct::excitedModes() const {
    std::vector<int> modes(modeCouplingLvl());
    std::transform(m_prod.begin(), m_prod.end(), modes.begin(), [](const t_Modal &in) {return in[0];});
    return modes;
}

bool HProduct::withinSpace(const VibSpace &vibSpace) {
    bool inside = modeCouplingLvl() <= vibSpace.modeCoupling;
    for (const auto &el : m_prod) {
        inside = inside && (el[0] < vibSpace.nMode && el[1] < vibSpace.nModal);
    }
    return inside;
}
