#include "gciHProduct.h"

#include <algorithm>
#include <stdexcept>

namespace gci {
HProduct::HProduct(product_t phi) : m_prod(std::move(phi)) {
    // Catch empty initialization with HProduct({ {} }) or HProduct({ })
    if (m_prod.size() == 1) if (m_prod[0].empty()) m_prod = product_t{};
    // Catch erroneous empty initialization with HProduct({ { {} } }). The extra {} initializes int to 0.
    if (m_prod.size() == 1) if (m_prod[0].size() == 1 && m_prod[0][0] == 0) m_prod = product_t{};
    order();
    check();
}

int HProduct::excLvl(int iMode) const {
    auto index = std::find_if(m_prod.cbegin(), m_prod.cend(), [iMode](const auto modal) {return modal[0] == iMode;});
    if (index == m_prod.end()) return 0;
    else return (*index)[1];
}

//TODO Testing!
HProduct HProduct::excite(const VibExcitation &exc) const {
    auto newProduct = HProduct(*this);
    // if any of the modes are different, excitation is void
    // for each mode, if annihilation operator does not match the modal,
    // set it to -1 and exit.
    for (int i = 0; i < exc.mc_lvl; ++i) {
        auto &exc_i = exc.excitations[i];
        bool found = false;
        for (auto &modal : newProduct) {
            if (modal[0] != exc_i[0]) continue;
            if (exc_i[2] != modal[1]) {
                found = true;
                modal[1] = -1;
                break;
            } else {
                found = true;
                modal[1] = exc_i[1];
            }
        }
        if (!found) {
            newProduct.m_prod.emplace_back(modal_t({-1, -1}));
            break;
        }
    }
    return newProduct;
}

void HProduct::order() {
    if (empty()) return;
    // remove any modalIndex == 0
    auto newEnd = std::remove_if(m_prod.begin(), m_prod.end(), [](const modal_t &modal) {return modal[1] == 0;});
    m_prod.erase(newEnd, m_prod.end());
    // sort by modeIndex
    std::sort(m_prod.begin(), m_prod.end(), [](modal_t &modalA, modal_t &modalB) {return modalA[0] < modalB[0];});
}

void HProduct::check() const {
    if (empty()) return;
    if (std::any_of(m_prod.cbegin(), m_prod.cend(), [](const modal_t &modal) {return modal[0] < 0 || modal[1] < 0;}))
        throw std::logic_error("Product indices must be positive integers.");
    for (auto modal = m_prod.begin(); modal != m_prod.end() - 1; ++modal) {
        if (modal[0][0] == modal[1][0])
            throw std::logic_error("HProduct with more than a single modal occupied per mode.");
    }
}

std::vector<int> HProduct::excitedModes() const {
    std::vector<int> modes(excLvl());
    std::transform(m_prod.cbegin(), m_prod.cend(), modes.begin(), [](const modal_t &in) {return in[0];});
    return modes;
}

bool HProduct::withinSpace(const VibSpace &vibSpace) const {
    if (excLvl() > vibSpace.excLvl) return false;
    return !std::any_of(m_prod.cbegin(), m_prod.cend(),
                        [&vibSpace](const modal_t &modal) {
                            return modal[0] >= vibSpace.nMode || modal[1] >= vibSpace.nModal || modal[1] < 0;
                        });
}

void HProduct::changeModal(const int iMode, const int diff) {
    auto modal = std::find_if(m_prod.begin(), m_prod.end(),
                              [&iMode](const modal_t &el) {return (el[0] == iMode);});
    if (modal == m_prod.end()) {
        m_prod.emplace_back(modal_t({iMode, diff}));
    } else {
        (*modal)[1] += diff;
    }
    order();
//    check();
}

std::ostream &operator<<(std::ostream &os, HProduct const &obj) {
    os << " {";
    for (auto const &el : obj) os << "{" << el[0] << "," << el[1] << "},";
    os << "} ";
    return os;
}

}  // namespace gci
