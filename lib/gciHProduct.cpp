#include "gciHProduct.h"

#include <algorithm>
#include <stdexcept>

namespace gci {
HProduct::HProduct(const t_Product &phi) : m_prod(phi) {
    // Catch empty initialization with HProduct({ {} }) or HProduct({ })
    if (m_prod.size() == 1) if (m_prod[0].empty()) m_prod = t_Product{};
    // Catch erroneous empty initialization with HProduct({ { {} } }). The extra {} initializes int to 0.
    if (m_prod.size() == 1) if (m_prod[0].size() == 1 && m_prod[0][0] == 0) m_prod = t_Product{};
    order();
    check();
}

int HProduct::excLvl(int iMode) const {
    auto index = std::find_if(m_prod.cbegin(), m_prod.cend(), [iMode](const auto modal) {return modal[0] == iMode;});
    if (index == m_prod.end()) return 0;
    else return (*index)[1];
}

void HProduct::order() {
    if (empty()) return;
    // remove any modalIndex == 0
    auto newEnd = std::remove_if(m_prod.begin(), m_prod.end(), [](const t_Modal &modal) {return modal[1] == 0;});
    m_prod.erase(newEnd, m_prod.end());
    // sort by modeIndex
    std::sort(m_prod.begin(), m_prod.end(), [](t_Modal &modalA, t_Modal &modalB) {return modalA[0] < modalB[0];});
}

void HProduct::check() const {
    if (empty()) return;
    if (std::any_of(m_prod.cbegin(), m_prod.cend(), [](const t_Modal &modal) {return modal[0] < 0 || modal[1] < 0;}))
        throw std::logic_error("Product indices must be positive integers.");
    for (auto modal = m_prod.begin(); modal != m_prod.end() - 1; ++modal) {
        if (modal[0][0] == modal[1][0])
            throw std::logic_error("HProduct with more than a single modal occupied per mode.");
    }
}

std::vector<int> HProduct::excitedModes() const {
    std::vector<int> modes(excLvl());
    std::transform(m_prod.cbegin(), m_prod.cend(), modes.begin(), [](const t_Modal &in) {return in[0];});
    return modes;
}

bool HProduct::withinSpace(const VibSpace &vibSpace) {
    return !std::any_of(m_prod.cbegin(), m_prod.cend(),
                        [&vibSpace](const t_Modal &modal) {
                            return modal[0] >= vibSpace.nMode || modal[1] >= vibSpace.nModal;
                        });
}

void HProduct::changeModal(const int iMode, const int diff) {
    auto modal = std::find_if(m_prod.begin(), m_prod.end(),
                              [&iMode](const t_Modal &el) {return (el[0] == iMode);});
    if (modal == m_prod.end()) {
        m_prod.emplace_back(t_Modal({iMode, diff}));
    } else {
        (*modal)[1] += diff;
    }
    order();
    check();
}

std::ostream &operator<<(std::ostream &os, HProduct const &obj) {
    os << " {";
    for (const auto &el : obj) os << "{" << el[0] << "," << el[1] << "},";
    os << "} ";
    return os;
}

}  // namespace gci
