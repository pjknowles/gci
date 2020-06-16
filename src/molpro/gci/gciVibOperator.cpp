#include "gciVibOperator.h"

namespace molpro {
namespace gci {

namespace ns_VibOperator {
size_t hash(const VibExcitation &exc, int nMode, int nModal, parity_t hermiticity, parity_t exchange) {
  if (exc.mc_lvl == 1) {
    if (hermiticity == parity_t::even || hermiticity == parity_t::odd)
      return hash_mc1_sym(exc, nModal);
    else
      return hash_mc1_nosym(exc, nModal);
  } else {
    throw std::logic_error("Hash number for MC > 1 not implemented yet");
  }
}

size_t hash_mc1_sym(const VibExcitation &exc, int nModal) {
  size_t A = 1 + exc.excitations[0][0];
  size_t n = 1 + std::max(exc.excitations[0][1], exc.excitations[0][2]);
  size_t m = 1 + std::min(exc.excitations[0][1], exc.excitations[0][2]);
  size_t h_per_mode = nModal * (nModal + 1) / 2;
  return (A - 1) * h_per_mode + n * (n - 1) / 2 + m;
}

size_t hash_mc1_nosym(const VibExcitation &exc, int nModal) {
  size_t A = 1 + exc.excitations[0][0];
  size_t i = 1 + exc.excitations[0][1];
  size_t j = 1 + exc.excitations[0][2];
  size_t h_per_mode = nModal * nModal;
  return (A - 1) * h_per_mode + (i - 1) * nModal + j;
}

size_t hash_mc1_nosym_old(const VibExcitation &exc, int nModal) {
  size_t A = 1 + exc.excitations[0][0];
  size_t i = 1 + exc.excitations[0][1];
  size_t j = 1 + exc.excitations[0][2];
  size_t h_per_mode = nModal * nModal;
  if (i >= j)
    return (A - 1) * h_per_mode + (i - 1) * (i - 1) + j;
  else
    return (A - 1) * h_per_mode + (j - 1) * (j - 1) + 2 * j - i;
}
} //  namespace ns_VibOperator
} // namespace gci
} // namespace molpro
