#include "SMatMat-implementation.h"

using T = molpro::SMat_<std::complex<double>>;

namespace molpro {
template <> typename molpro::bytestream SMatMat_<T>::bytestream(bool data) { throw std::logic_error("unimplemented"); }

template <> SMatMat_<T>::SMatMat_(const char* dump, typename T::value_type* buffer) {
  throw std::logic_error("unimplemented");
}
} // namespace molpro

template class molpro::SMatMat_<T>;
