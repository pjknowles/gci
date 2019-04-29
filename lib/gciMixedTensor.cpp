#include "gciMixedTensor.h"

namespace gci {

void MixedTensor::append(SymmetryMatrix::Operator &&op, const VibExcitation &vibExc) {
    auto h = hash(vibExc);
    auto &&el = tensor_el_t({op, vibExc});
    tensor.insert({h, el});
}

SymmetryMatrix::Operator &MixedTensor::at(const VibExcitation &exc) {
    auto h = hash(exc);
    return tensor.at(h).first;
}

size_t MixedTensor::hash(const VibExcitation &exc, int nMode, int nModal, parity_t hermiticity, parity_t exchange) {
    if (exc.mc_lvl() == 1) {
        if (hermiticity == parity_t::even || hermiticity == parity_t::odd) return hash_mc1_sym(exc, nModal);
        else return hash_mc1_nosym(exc, nModal);
    } else {
        throw std::logic_error("Hash number for MC > 1 not implemented yet");
    }
}

size_t MixedTensor::hash_mc1_sym(const VibExcitation &exc, int nModal) {
    size_t A = exc.modes()[0];
    size_t n = std::max(exc.excitations()[0].first, exc.excitations()[0].second);
    size_t m = std::min(exc.excitations()[0].first, exc.excitations()[0].second);
    size_t h_per_mode = nModal * (nModal + 1) / 2;
    return (A - 1) * h_per_mode + n * (n - 1) / 2 + m;
}

size_t MixedTensor::hash_mc1_nosym(const VibExcitation &exc, int nModal) {
    size_t A = exc.modes()[0];
    size_t i = exc.excitations()[0].first;
    size_t j = exc.excitations()[0].second;
    size_t h_per_mode = nModal * nModal;
    return (A - 1) * h_per_mode + (i - 1) * nModal + j;
}

size_t MixedTensor::hash_mc1_nosym_old(const VibExcitation &exc, int nModal) {
    size_t A = exc.modes()[0];
    size_t i = exc.excitations()[0].first;
    size_t j = exc.excitations()[0].second;
    size_t h_per_mode = nModal * nModal;
    if (i >= j) return (A - 1) * h_per_mode + (i - 1) * (i - 1) + j;
    else return (A - 1) * h_per_mode + (j - 1) * (j - 1) + 2 * j - i;
}
} // namespace gci
