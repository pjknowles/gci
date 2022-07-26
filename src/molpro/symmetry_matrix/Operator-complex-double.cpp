#include <complex>

template <class T> class Constants {
public:
  const T s_zero = 0;
  T return_value;
};

using T = std::complex<double>;

static Constants<T> constants; // FIXME have to find a better way than this, because it's not thread-safe
#include "Operator-implementation.h"

namespace molpro {
template <> molpro::bytestream Operator_<T>::bytestream() const { throw std::logic_error("unimplemented"); }
} // namespace molpro

template class molpro::Operator_<T>;
