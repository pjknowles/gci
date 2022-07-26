#include "SMatMat-implementation.h"

using T = SymmetryMatrix::SMat_<float>;
namespace SymmetryMatrix {
template <> memory::bytestream SMatMat_<T>::bytestream(bool data) { throw std::logic_error("unimplemented"); }
} // namespace SymmetryMatrix
template class SymmetryMatrix::SMatMat_<T>;
