

template <class T> class Constants {
public:
  const T s_zero = 0;
  T return_value;
};

using T = double;

static Constants<T> constants; // FIXME have to find a better way than this, because it's not thread-safe
#include "Operator-implementation.h"

template class molpro::Operator_<T>;
