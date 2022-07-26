template <class T> class Constants {
public:
  const T s_zero = 0;
  T return_value;
};

using T = float;

static Constants<T> constants; // FIXME have to find a better way than this, because it's not thread-safe
#include "Operator-implementation.h"

namespace SymmetryMatrix {
template <> memory::bytestream Operator_<T>::bytestream() { throw std::logic_error("unimplemented"); }
} // namespace SymmetryMatrix

template class SymmetryMatrix::Operator_<T>;
