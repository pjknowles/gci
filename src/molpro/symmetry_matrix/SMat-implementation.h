#include "SMat.h"
#include <molpro/iostream.h>

static constexpr bool sMat_debug = false;

#define ErrorExit(Why) throw std::logic_error(#Why);
#include <iostream>
typedef int32_t FORTINT; // assumes LP64
typedef double FORTDBL;
typedef FORTINT const& FINTARG;
typedef FORTDBL const& FDBLARG;

extern "C" {
#ifdef SMAT_USE_BLAS
void dgemm(char const& TransA, char const& TransB, FINTARG M, FINTARG N, FINTARG K, FORTDBL const& Alpha,
           FORTDBL const* A, FINTARG lda, FORTDBL const* B, FINTARG ldb, FORTDBL const& Beta, FORTDBL* C, FINTARG ldc);

void dgemv(char const& Trans, FINTARG M, FINTARG N, FORTDBL const& Alpha, FORTDBL const* A, FINTARG lda,
           FORTDBL const* X, FINTARG incx, FORTDBL const& Beta, FORTDBL* Y, FINTARG incy);

#endif
}

using namespace molpro;

namespace molpro {

template <class T>
SMat_<T>::SMat_(dims_t dimensions, parity_t parity, int symmetry, bool diagonal, std::string description)
    : m_dimensions(std::move(dimensions)), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(true), m_transposed(false), m_diagonal(diagonal) {
  if (sMat_debug)
    molpro::cout << "constructor 1" << std::endl;
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  for (size_t axis = 0; axis < rank(); axis++)
    while (this->m_dimensions[axis].size() < 8)
      this->m_dimensions[axis].push_back(0);
  //   this->buffer = new molpro::array<double>(size());
  size_t n; // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(n = size());
  this->m_buffer = &this->m_bufferp.get()[0];
  m_managed_buffer = true;
}

template <class T>
SMat_<T>::SMat_(dims_t dimensions, molpro::array<value_type>& buffer, parity_t parity, int symmetry, bool diagonal,
                std::string description)
    : m_dimensions(std::move(dimensions)), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(false), m_transposed(false), m_diagonal(false),
      m_buffer(&buffer) {
  if (sMat_debug)
    molpro::cout << "constructor 2" << std::endl;
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  for (size_t axis = 0; axis < rank(); axis++)
    while (this->m_dimensions[axis].size() < 8)
      this->m_dimensions[axis].push_back(0);
  if (buffer.size() != size())
    ErrorExit("bad attached buffer size");
}

template <class T>
SMat_<T>::SMat_(dims_t dimensions, value_type* buffer, parity_t parity, int symmetry, bool diagonal,
                std::string description)
    : m_dimensions(std::move(dimensions)), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(false), m_transposed(false), m_diagonal(diagonal) {
  if (sMat_debug)
    molpro::cout << "constructor 3" << std::endl;
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  for (size_t axis = 0; axis < rank(); axis++)
    while (this->m_dimensions[axis].size() < 8)
      this->m_dimensions[axis].push_back(0);
  //     this->buffer = new molpro::array<T>(buffer,size()); // no checks on size of buffer!
  size_t n;                                                                 // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(buffer, n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
}

template <class T>
SMat_<T>::SMat_(std::string space, molpro::parity_t parity, int symmetry, bool diagonal, std::string description)
    : m_dimensions(SymmetryMatrix::spaces(std::move(space))), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(true), m_transposed(false), m_diagonal(diagonal) {
  size_t n;                                                         // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
}

template <class T>
SMat_<T>::SMat_(std::string space, molpro::array<T>& buffer, molpro::parity_t parity, int symmetry, bool diagonal,
                std::string description)
    : m_dimensions(SymmetryMatrix::spaces(std::move(space))), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(false), m_transposed(false), m_diagonal(diagonal),
      m_buffer(&buffer) {}

template <class T>
SMat_<T>::SMat_(std::string space, T* buffer, molpro::parity_t parity, int symmetry, bool diagonal,
                std::string description)
    : m_dimensions(SymmetryMatrix::spaces(std::move(space))), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)), m_managed_buffer(false), m_transposed(false), m_diagonal(diagonal) {
  size_t n;                                                                 // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(buffer, n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
}

template <class T>
SMat_<T>::SMat_(SMat_<T> const* source, molpro::parity_t parity, int symmetry, unsigned int rank, bool diagonal,
                std::string description)
    : m_dimensions(source->m_dimensions), m_parity(parity == parityUnspecified ? source->m_parity : parity),
      m_symmetry(symmetry == 9 ? source->m_symmetry : symmetry), m_description(std::move(description)),
      m_managed_buffer(true), m_transposed(false), m_diagonal(diagonal) {
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  if (rank > 2)
    m_transposed = not source->m_transposed;
  if (rank == 1 && m_dimensions.size() == 2)
    m_dimensions.pop_back();
  if (rank == 2 && m_dimensions.size() == 1)
    m_dimensions.push_back(m_dimensions[0]);
  size_t n;                                                         // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(n = size()); // no checks on size of buffer!
  m_buffer = &this->m_bufferp.get()[0];
  if (sMat_debug)
    molpro::cout << "constructor 4, buffer at " << &(*m_buffer)[0] << std::endl;
}

template <class T>
SMat_<T>::SMat_(SMat_<T> const* source, molpro::array<T>& buffer, molpro::parity_t parity, int symmetry,
                unsigned int rank, bool diagonal, std::string description)
    : m_dimensions(source->m_dimensions), m_parity(parity == parityUnspecified ? source->m_parity : parity),
      m_symmetry(symmetry == 9 ? source->m_symmetry : symmetry), m_description(std::move(description)),
      m_managed_buffer(false), m_transposed(source->m_transposed), m_diagonal(diagonal), m_buffer(&buffer) {
  if (sMat_debug)
    molpro::cout << "constructor 5" << std::endl;
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  if (rank > 2)
    m_transposed = not source->m_transposed;
  if (rank == 1 && m_dimensions.size() == 2)
    m_dimensions.pop_back();
  if (rank == 2 && m_dimensions.size() == 1)
    m_dimensions.push_back(m_dimensions[0]);
}

template <class T>
SMat_<T>::SMat_(SMat_<T> const* source, T* buffer, molpro::parity_t parity, int symmetry, unsigned int rank,
                bool diagonal, std::string description)
    : m_dimensions(source->m_dimensions), m_parity(parity == parityUnspecified ? source->m_parity : parity),
      m_symmetry(symmetry == 9 ? source->m_symmetry : symmetry), m_description(std::move(description)),
      m_managed_buffer(false), m_transposed(source->m_transposed), m_diagonal(diagonal) {
  if (sMat_debug)
    molpro::cout << "constructor 6" << std::endl;
  if (symmetry < 0)
    throw std::logic_error("obsolete symmetry=-1: use diagonal instead");
  if (rank > 2)
    m_transposed = not source->m_transposed;
  if (rank == 1 && m_dimensions.size() == 2)
    m_dimensions.pop_back();
  if (rank == 2 && m_dimensions.size() == 1)
    m_dimensions.push_back(m_dimensions[0]);
  size_t n;                                                                 // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(buffer, n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
}

template <class T>
SMat_<T>::SMat_(SMat_<T> const& source, int option)
    : m_dimensions(source.m_dimensions), m_parity(source.m_parity), m_symmetry(source.m_symmetry),
      m_description(source.m_description), m_managed_buffer(true), m_transposed(source.m_transposed),
      m_diagonal(source.m_diagonal) {
  if (sMat_debug)
    molpro::cout << source.str("source in constructor 7", 1);
  size_t n; // For dumb pgc++ - bug 5086
  this->m_bufferp = std::make_shared<molpro::array<T>>(n = size());
  this->m_buffer = &this->m_bufferp.get()[0];
  std::copy(source.m_buffer->begin(), source.m_buffer->end(), this->m_buffer->begin());
  if (sMat_debug)
    molpro::cout << str("this in constructor 7", 1);
}

template <class T> SMat_<T>::SMat_(const char* dump, T* buffer) {
  throw std::logic_error("Construction of non-double SMat_ from bytestream is not implemented");
}

template <class T> molpro::array<T> SMat_<T>::block(unsigned int block_symmetry) const {
  if (sMat_debug)
    molpro::cout << block_symmetry << " " << block_offset(block_symmetry) << " " << block_size(block_symmetry) << " "
                 << &((*m_buffer)[block_offset(block_symmetry)]) << " " << (block_dimensions(block_symmetry)[0]) << " "
                 << block_dimensions(block_symmetry)[1] << std::endl;
  //  auto result = molpro::array<T>(&((*m_buffer)[block_offset(block_symmetry)]), block_size(block_symmetry));
  //  molpro::cout << result.data() << std::endl;
  //  return result;
  return molpro::array<T>(&((*m_buffer)[block_offset(block_symmetry)]), block_size(block_symmetry));
}

#ifdef EIGEN_CORE_H
using namespace Eigen;
// the result will unconditionally be an Eigen Matrix where the row index is from block_symmetry
// if the matrix is represented internally as transpose, then number of rows is m_dimensions[1][block_symmetry],
// otherwise it is m_dimensions[0][block_symmetry]
template <class T> typename SMat_<T>::M SMat_<T>::blockMap(unsigned int block_symmetry) {
  Mconst r = static_cast<const SMat_&>(*this).blockMap(block_symmetry);
  return M(const_cast<T*>(r.data()), r.rows(), r.cols(), Stride<Dynamic, Dynamic>(r.outerStride(), r.innerStride()));
}

template <class T>
typename SMat_<T>::Mconst SMat_<T>::blockMap(unsigned int block_symmetry) const // FIXME - push to header
// FIXME implement transposition opaquely, consistently across the whole class
{
  if (rank() == 1) // throw std::logic_error("SMat::blockMap called on a vector"); // FIXME this should be nevertheless
                   // OK and mappable
    return Mconst(&((*m_buffer)[block_offset(block_symmetry)]), m_dimensions[0][block_symmetry], 1);

  unsigned int block = m_transposed ? block_symmetry ^ m_symmetry : block_symmetry;
  unsigned int row_axis = m_transposed ? 1 : 0;
  unsigned int col_axis = 1 - row_axis;
  unsigned int row_sym = block_symmetry;
  unsigned int col_sym = row_sym ^ m_symmetry;
  size_t nrows = m_dimensions[row_axis][row_sym];
  size_t ncols = m_dimensions[col_axis][col_sym];
  bool tr = m_parity == parityNone
                ? m_transposed
                : (row_sym < col_sym); // whether the Eigen matrix map will turn out to be transposed or regular
  size_t row_stride = tr ? ncols : 1;
  size_t col_stride = tr ? 1 : nrows;
  //  std::cout << "SMat::blockMap block="<<block<<", block_symmetry="<<block_symmetry<<", tr="<<tr<<",
  //  col_stride="<<col_stride<<", row_stride="<<row_stride<<std::endl;
  if (blockMapPossible(block_symmetry))
    return Mconst(&((*m_buffer)[block_offset(block)]), nrows, ncols, Stride<Dynamic, Dynamic>(col_stride, row_stride));
  std::cout << "block_symmetry=" << block_symmetry << ", block=" << block << std::endl;
  throw std::logic_error("SMat::blockMap called for a symmetry block that cannot be memory-mapped");
}

template <class T>
bool SMat_<T>::blockNotEmpty(unsigned int block_symmetry) const // FIXME - push to header
{
  if (rank() == 1)
    return m_dimensions[0][block_symmetry] > 0;

  unsigned int row_axis = m_transposed ? 1 : 0;
  unsigned int col_axis = 1 - row_axis;
  unsigned int row_sym = block_symmetry;
  unsigned int col_sym = row_sym ^ m_symmetry;
  size_t nrows = m_dimensions[row_axis][row_sym];
  size_t ncols = m_dimensions[col_axis][col_sym];
  return nrows * ncols > 0;
}

template <class T>
bool SMat_<T>::blockMapPossible(unsigned int block_symmetry) const // FIXME - push to header
{
  if (rank() == 1)
    throw std::logic_error("SMat::blockMap called on a vector");
  unsigned int block = m_transposed ? block_symmetry ^ m_symmetry : block_symmetry;
  return (rank() == 1 || m_parity == parityNone || block > (block ^ m_symmetry) ||
          (m_parity != parityNone && m_symmetry > 0));
}

template <class T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> SMat_<T>::blockCopy(unsigned int block_symmetry) const {
  unsigned int block = m_transposed ? block_symmetry ^ m_symmetry : block_symmetry;
  if (blockMapPossible(block_symmetry))
    return Eigen::Matrix<T, Eigen::Dynamic,
                         Eigen::Dynamic>(
        blockMap(block_symmetry)); // TODO check that this really returns reference not copy
                                   //  std::cout << "SMat::blockCopy "<<block_symmetry<<" has to construct"<<std::endl;

  //   need to make a copy
  unsigned int row_axis = m_transposed ? 1 : 0;
  unsigned int col_axis = 1 - row_axis;
  unsigned int row_sym = block_symmetry;
  unsigned int col_sym = row_sym ^ m_symmetry;
  size_t nrows = m_dimensions[row_axis][row_sym];
  size_t ncols = m_dimensions[col_axis][col_sym];
  bool tr = m_transposed;
  if (m_parity == parityEven && row_sym < col_sym)
    tr = !tr;
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> e(nrows, ncols);
  size_t o = block_offset(block);
  if (m_symmetry == 0 && m_parity == parityEven) {
    for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j <= i; j++)
        e(j, i) = e(i, j) = (*m_buffer)[o++];
  } else if (m_symmetry == 0 && m_parity == parityOdd && m_transposed) {
    for (size_t i = 0; i < nrows; i++) {
      for (size_t j = 0; j <= i; j++)
        e(i, j) = -(e(j, i) = (*m_buffer)[o++]);
      e(i, i) = 0;
    }
  } else if (m_symmetry == 0 && m_parity == parityOdd) {
    for (size_t i = 0; i < nrows; i++) {
      for (size_t j = 0; j <= i; j++)
        e(j, i) = -(e(i, j) = (*m_buffer)[o++]);
      e(i, i) = 0;
    }
  } else if (m_symmetry == 0 && m_parity == parityOddPacked && m_transposed) {
    for (size_t i = 0; i < nrows; i++) {
      for (size_t j = 0; j < i; j++)
        e(i, j) = -(e(j, i) = (*m_buffer)[o++]);
      e(i, i) = 0;
    }
  } else if (m_symmetry == 0 && m_parity == parityOddPacked) {
    for (size_t i = 0; i < nrows; i++) {
      for (size_t j = 0; j < i; j++)
        e(j, i) = -(e(i, j) = (*m_buffer)[o++]);
      e(i, i) = 0;
    }
  } else if (m_transposed &&
             (m_parity == parityOdd || m_parity == parityOddPacked)) { // odd parity, non-symmetric, transpose block
    for (size_t j = 0; j < ncols; j++)
      for (size_t i = 0; i < nrows; i++)
        e(i, j) = -(*m_buffer)[o++];
  } else if (m_parity == parityOdd || m_parity == parityOddPacked) { // odd parity, non-symmetric, transpose block
    for (size_t i = 0; i < nrows; i++)
      for (size_t j = 0; j < ncols; j++)
        e(i, j) = -(*m_buffer)[o++];
  } else // should be never
    throw std::logic_error("Unimplemented case in SMat_::blockCopy");
  return e;
}

template <class T>
void SMat_<T>::blockImport(unsigned int block_symmetry, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& e) {
  unsigned int block = m_transposed ? block_symmetry ^ m_symmetry : block_symmetry;
  unsigned int row_axis = m_transposed ? 1 : 0;
  unsigned int col_axis = 1 - row_axis;
  unsigned int row_sym = block_symmetry;
  unsigned int col_sym = row_sym ^ m_symmetry;
  int nrows = m_dimensions[row_axis][row_sym];
  int ncols = m_dimensions[col_axis][col_sym];
  bool tr = m_transposed;
  if (m_parity == parityEven && row_sym < col_sym)
    tr = !tr;
  assert(e.rows() == nrows && e.cols() == ncols);
  if (blockMapPossible(block_symmetry)) {
    blockMap(block_symmetry) = e;
    return;
  }
  size_t o = block_offset(block);
  if (m_symmetry == 0 && m_parity == parityEven) {
    for (auto i = 0; i < nrows; i++)
      for (auto j = 0; j <= i; j++)
        (*m_buffer)[o++] = e(j, i);
  } else if (m_symmetry == 0 && m_parity == parityOdd && m_transposed) {
    for (auto i = 0; i < nrows; i++) {
      for (auto j = 0; j <= i; j++)
        (*m_buffer)[o++] = e(j, i);
    }
  } else if (m_symmetry == 0 && m_parity == parityOdd) {
    for (auto i = 0; i < nrows; i++) {
      for (auto j = 0; j <= i; j++)
        (*m_buffer)[o++] = e(i, j);
    }
  } else if (m_symmetry == 0 && m_parity == parityOddPacked && m_transposed) {
    for (auto i = 0; i < nrows; i++) {
      for (auto j = 0; j < i; j++)
        (*m_buffer)[o++] = e(j, i);
    }
  } else if (m_symmetry == 0 && m_parity == parityOddPacked) {
    for (auto i = 0; i < nrows; i++) {
      for (auto j = 0; j < i; j++)
        (*m_buffer)[o++] = e(i, j);
    }
  } else if (m_parity == parityOdd || m_parity == parityOddPacked) { // odd parity, non-symmetric, transpose block
    for (auto i = 0; i < nrows; i++)
      for (auto j = 0; j < ncols; j++)
        (*m_buffer)[o++] = -e(i, j);
  } else // should be never
    throw std::logic_error("Unimplemented case in SMat_::blockImport");
}

template <class T>
typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> SMat_<T>::blockV(unsigned int block_symmetry) const {
  if (rank() == 2 and not m_diagonal)
    throw std::logic_error("SMat::blockV called on a matrix");
  Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> result(&((*m_buffer)[block_offset(block_symmetry)]),
                                                  m_dimensions[0][block_symmetry]);
  return result;
}
#endif

template <class T> SMat_<T>::~SMat_() {
  if (sMat_debug) {
    molpro::cout << str("object SMat_ destructor", 6);
    memory_print_status();
    molpro::cout << "use count " << m_bufferp.use_count() << std::endl;
  }
  //  if (managed_buffer && buffer != nullptr)
  if (m_managed_buffer)
    m_bufferp.reset();
  if (sMat_debug) {
    memory_print_status();
    molpro::cout << "use count " << m_bufferp.use_count() << std::endl;
  }
}

template <class T> size_t SMat_<T>::size() const {
  size_t result = 0;
  if (m_diagonal)
    result = std::accumulate(m_dimensions[0].begin(), m_dimensions[0].end(), 0);
  else if (rank() == 1)
    result = m_dimensions[0][m_symmetry];
  else
    result = block_offset(8);
  return result;
}

template <class T> size_t SMat_<T>::block_size(unsigned int block_symmetry) const {
  //  molpro::cout << "block_size parity="<<m_parity<<", symmetry="<<m_symmetry<<", block_symmetry="<<block_symmetry<<",
  //  dimensions="<<m_dimensions[0][block_symmetry]<<m_dimensions[1][block_symmetry^m_symmetry]<<std::endl;
  if (m_diagonal)
    return std::accumulate(m_dimensions[0].begin(), m_dimensions[0].end(), 0);
  else if (rank() == 2) {
    if (m_parity != 0 && m_symmetry == 0)
      return m_dimensions[0][block_symmetry] *
             (m_dimensions[0][block_symmetry] + (m_parity == parityOddPacked ? -1 : 1)) / 2;
    else
      return m_dimensions[0][block_symmetry] * m_dimensions[1][block_symmetry ^ m_symmetry];
  } else if (rank() == 1) {
    return m_dimensions[0][block_symmetry];
  }
  return 0;
}

template <class T> size_t SMat_<T>::block_offset(unsigned int block_symmetry) const {
  size_t result = 0;
  if (m_diagonal)
    for (unsigned int ks = 0; ks < block_symmetry; ks++)
      result += m_dimensions[0][ks];
  else if (rank() == 2) {
    unsigned int bs = block_symmetry, bs2 = m_symmetry ^ block_symmetry;
    if (bs < max_symmetry_ && m_parity != 0 && bs2 > bs) {
      bs2 = block_symmetry;
      bs = m_symmetry ^ bs2;
    }
    for (unsigned int ks = 0; ks < bs; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      if (m_parity != 0 && m_symmetry == 0)
        result += m_dimensions[0][ks] * (m_dimensions[0][ks] + (m_parity == parityOddPacked ? -1 : 1)) / 2;
      else if (ks > ls || m_parity == 0)
        result += m_dimensions[0][ks] * m_dimensions[1][ls];
    }
  } else if (rank() == 1) {
    if (block_symmetry > (unsigned int)m_symmetry)
      result = m_dimensions[0][m_symmetry]; // usual use: block_symmetry=9 for total size
  }
  return result;
}

template <class T> std::vector<size_t> SMat_<T>::block_dimensions(unsigned int block_symmetry) const {
  dim_t result;
  if (m_diagonal) {
    result.resize(1, 0);
    result[0] = std::accumulate(m_dimensions[0].begin(), m_dimensions[0].end(), 0);
  } else if (rank() == 2) {
    size_t transpose = ((block_symmetry ^ (unsigned int)m_symmetry) > block_symmetry && m_parity != 0) ? 1 : 0;
    if (transpose) {
      result.push_back(m_dimensions[0][block_symmetry ^ m_symmetry]);
      result.push_back(m_dimensions[0][block_symmetry]);
      result.push_back(1);
    } else {
      result.push_back(m_dimensions[0][block_symmetry]);
      result.push_back(m_dimensions[1][block_symmetry ^ m_symmetry]);
      result.push_back(0);
    }
  } else if (rank() == 1) {
    result.resize(1, 0);
    if (block_symmetry == (unsigned int)m_symmetry)
      result[0] = m_dimensions[0][m_symmetry];
    result.push_back(1);
    result.push_back(0);
  }
  return result;
}

template <class T> bool SMat_<T>::block_transposed(unsigned int block_symmetry) const {
  return rank() != 1 && (block_symmetry ^ (unsigned int)m_symmetry) >= block_symmetry && m_parity != 0;
}

template <class T> T SMat_<T>::trace() const {
  if ((rank() != 2) || (m_symmetry != 0) || (m_parity < 0))
    return (T)0;
  auto result = (T)0;
  for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
    T* buff = &(*m_buffer)[block_offset(ks)];
    for (size_t i = 0; i < m_dimensions[0][ks]; i++) {
      result += *buff;
      buff += (m_diagonal) ? 1 : (m_parity > 0) ? i + 2 : m_dimensions[0][ks] + 1;
    }
  }
  return result;
}

template <class T> typename SMat_<T>::scalar_type SMat_<T>::norm() const {
  scalar_type result = 0;
  for (size_t k = 0; k < size(); k++)
    result += (*m_buffer)[k] * (*m_buffer)[k];
  if (rank() == 2 && m_parity != 0 && !m_diagonal)
    result *= 2;
  if (rank() == 2 && m_parity > 0 && !m_diagonal)
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      T* buff = &(*m_buffer)[block_offset(ks)];
      for (size_t i = 0; i < m_dimensions[0][ks]; i++) {
        result -= *buff;
        buff += (m_parity > 0) ? i + 2 : m_dimensions[0][ks] + 1;
      }
    }
  return std::sqrt(result);
}

template <class T> void SMat_<T>::assign(T value) {
  if (m_buffer == nullptr)
    return;
  for (size_t i = 0; i < size(); i++)
    (*m_buffer)[i] = value;
}

template <class T> void SMat_<T>::setIdentity() {
  if (m_symmetry != 0)
    ErrorExit("SMat::setIdentity has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::setIdentity has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::setIdentity has received a non-square matrix");
  assign(0);
  for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
    T* buff = &(*m_buffer)[block_offset(ks)];
    for (size_t i = 0; i < m_dimensions[0][ks]; i++) {
      *buff = 1;
      buff += m_diagonal ? 1 : (m_parity > 0) ? i + 2 : m_dimensions[0][ks] + 1;
    }
  }
}

template <class T> std::string SMat_<T>::str(std::string title, int level, int precision) const {
  if (level < 0)
    return "";
  std::stringstream s;
  s << (title.empty() ? m_description : title);
  if (m_symmetry >= 0)
    s << "; symmetry=" << m_symmetry + 1;
  if (rank() > 1)
    s << " parity=" << m_parity;
  if (m_diagonal)
    s << ", diagonal";
  if (level > 0) {
    if (m_buffer != nullptr)
      s << (m_managed_buffer ? "; managed" : " unmanaged") << " data at address " << &(*m_buffer)[0];
    else
      s << "; no data attached";
  }
  if (m_transposed)
    s << "; transposed data are stored";
  s << std::endl;
  if (level > 1 && m_buffer != nullptr) {
    s << "raw data:";
    for (size_t k = 0; k < size(); k++)
      s << " " << (*m_buffer)[k];
    s << std::endl;
  }
  for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
    size_t nr = m_dimensions[0][ks], nc = (rank() > 1 ? m_dimensions[1][ks ^ m_symmetry] : 1);
    if (!blockNotEmpty(ks))
      continue;
    if (rank() == 1 && m_symmetry >= 0 && (unsigned int)m_symmetry != ks)
      continue;
    if (rank() == 1 || m_diagonal)
      s << "Block (" << ks + 1 << "), dimensions (" << nr << ")" << std::endl;
    else
      s << "Block (" << ks + 1 << "," << (ks ^ m_symmetry) + 1 << "), dimensions (" << nr << "," << nc << ")"
        << std::endl;
    if (m_buffer == nullptr || level < 0)
      continue;
    if (rank() == 1 || m_diagonal)
      s << blockV(ks).transpose().format(Eigen::IOFormat(precision)) << std::endl; // use Eigen to do the printing
    else
      s << blockCopy(ks).format(Eigen::IOFormat(precision)) << std::endl; // use Eigen to do the printing
  }

  return s.str();
}

template <class T> SMat_<T>& SMat_<T>::copy(SMat_<T> const& source, dims_t sourceOffset) {
  if (sMat_debug)
    molpro::cout << "SMat::operator= symmetry()=" << symmetry() << std::endl;
  if (sMat_debug)
    molpro::cout << "SMat::operator= source.symmetry()=" << source.symmetry() << std::endl;
  if (this == &source)
    return *this;
  for (auto& o : sourceOffset)
    while (o.size() < max_symmetry_)
      o.push_back(0);

  //  molpro::cout <<
  //       (rank() == 2)
  //      << (source.rank() == 1)
  //      << (symmetry() == 0)
  //      << (source.symmetry() == -1)
  //      << (m_dimensions[0] == m_dimensions[1])
  //      << (m_dimensions[0] == source.m_dimensions[0])
  //      << std::endl;
  if (rank() == 2 and source.rank() == 2 and symmetry() == 0 and source.symmetry() == 0 and source.m_diagonal and
      !m_diagonal and m_dimensions[0] == m_dimensions[1] and
      m_dimensions[0] == source.m_dimensions[0]) { // make explicitly expanded diagonal matrix from diagonal matrix
    assign(0);
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      if (m_parity == parityNone) {
        for (size_t i = 0; i < m_dimensions[0][ks]; i++)
          blockMap(ks)(i, i) = source.block(ks)[i];
      } else if (m_parity == parityEven) {
        for (size_t i = 0; i < m_dimensions[0][ks]; i++)
          block(ks)[(i + 1) * (i + 2) / 2 - 1] = source.block(ks)[i];
      } else
        throw std::logic_error("unimplemented");
    }
    return *this;
  }
  if (rank() == 2 and source.rank() == 2 and symmetry() == 0 and source.symmetry() == 0 and !source.m_diagonal and
      m_diagonal and m_dimensions[0] == m_dimensions[1] and
      m_dimensions[0] == source.m_dimensions[0]) { // make diagonal matrix from expanded but diagonal matrix
    *this = source.diagonal();
    return *this;
  }
  if (m_diagonal and source.m_diagonal && std::accumulate(sourceOffset[0].begin(), sourceOffset[0].end(), 0) == 0 &&
      std::accumulate(sourceOffset[1].begin(), sourceOffset[1].end(), 0) == 0) {
    std::copy(source.m_buffer->begin(), source.m_buffer->end(), m_buffer->begin());
    return *this;
  }
  if (m_diagonal or source.m_diagonal)
    throw std::logic_error("incorrectly coded diagonal matrix");
  if (source.symmetry() != this->symmetry())
    ErrorExit("SMat::operator=: result matrix has wrong symmetry");
  if (rank() != source.rank())
    ErrorExit("Mismatch in ranks, SMat_::operator=");
  bool matching_dimensions = this->m_parity == source.m_parity;
  for (size_t k = 0; k < m_dimensions.size(); k++) {
    matching_dimensions &= m_dimensions[k] == source.m_dimensions[k] &&
                           std::accumulate(sourceOffset[k].begin(), sourceOffset[k].end(), 0) == 0;
  }
  if (this->m_description.empty())
    this->m_description = source.m_description;
  //  if (rank() == 2 && (source.m_parity || this->m_parity) && m_dimensions[0] != m_dimensions[1])
  //    ErrorExit("SMat::operator=: matrix must be square"); // FIXME delete this when ready

  if (rank() == 1 || matching_dimensions) { // straight copy
    this->m_transposed = source.m_transposed;
    for (size_t i = 0; i < size(); i++) {
      (*m_buffer)[i] = (*source.m_buffer)[i];
    }
    return *this;
  }

  // at this point, rank=2 and explicit conversion is needed
  bool square = m_symmetry == 0;
  for (size_t k = 0; k < sourceOffset[0].size(); k++)
    square = square and sourceOffset[0][k] == sourceOffset[1][k] and m_dimensions[0][k] == m_dimensions[1][k];
  //  molpro::cout << "m_dimensions[0]"; for (const auto& o : m_dimensions[0]) molpro::cout << " "<<o; molpro::cout <<
  //  std::endl; molpro::cout << "m_dimensions[1]"; for (const auto& o : m_dimensions[1]) molpro::cout << " "<<o;
  //  molpro::cout << std::endl; molpro::cout << "sourceOffset[0]"; for (const auto& o : sourceOffset[0]) molpro::cout
  //  << " "<<o; molpro::cout << std::endl; molpro::cout << "sourceOffset[1]"; for (const auto& o : sourceOffset[1])
  //  molpro::cout << " "<<o; molpro::cout << std::endl; molpro::cout << "square="<<square<<std::endl; molpro::cout <<
  //  "parity "<<m_parity<<source.m_parity<<std::endl; molpro::cout << "symmetry
  //  "<<m_symmetry<<source.m_symmetry<<std::endl; molpro::cout << "explicit in
  //  copy"<<this->m_parity<<source.m_parity<<std::endl;
  this->m_transposed = false;
  //  molpro::cout << "operator= source.transposed_="<<source.transposed_<<std::endl;
  if (this->m_parity != parityNone && source.m_parity != parityNone) { // copy triangle to triangle
    //    throw std::logic_error("triangle to triangle not yet done");
    //    molpro::cout << "triangle to triangle"<<std::endl;
    if (this->m_parity * source.m_parity < 0)
      throw std::logic_error("Parity flip in copy illegal");
    //    if (this->m_parity != source.m_parity)
    //      throw std::logic_error("Parity conversion not implemented"); // this would be parityOdd / parityOddPacked
    this->m_transposed = source.m_transposed; // TODO is the logic that follows working for the transposed case?
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      molpro::array<T> from = source.block(ks);
      size_t nkb = source.m_dimensions[0][ks];
      size_t nk = std::min(m_dimensions[0][ks], source.m_dimensions[0][ks] - sourceOffset[0][ks]);
      size_t nl = std::min(m_dimensions[1][ls], source.m_dimensions[1][ls] - sourceOffset[1][ls]);
      molpro::array<T> to = this->block(ks);
      if (ks == ls) {
        if (m_parity == parityOddPacked) {
          if (source.m_parity == parityOddPacked)
            for (size_t k = 0; k < nk; k++)
              for (size_t l = 0; l < k; l++)
                to[k * (k - 1) / 2 + l] =
                    from[(k + sourceOffset[0][ks]) * (k + sourceOffset[0][ks] - 1) / 2 + (l + sourceOffset[1][ks])];
          else
            for (size_t k = 0; k < nk; k++)
              for (size_t l = 0; l < k; l++)
                to[k * (k - 1) / 2 + l] =
                    from[(k + sourceOffset[0][ks]) * (k + sourceOffset[0][ks] + 1) / 2 + (l + sourceOffset[1][ks])];
        } else {
          if (source.m_parity == parityOddPacked)
            for (size_t k = 0; k < nk; k++)
              for (size_t l = 0; l <= k; l++)
                to[k * (k + 1) / 2 + l] =
                    from[(k + sourceOffset[0][ks]) * (k + sourceOffset[0][ks] - 1) / 2 + (l + sourceOffset[1][ks])];
          else
            for (size_t k = 0; k < nk; k++)
              for (size_t l = 0; l <= k; l++)
                to[k * (k + 1) / 2 + l] =
                    from[(k + sourceOffset[0][ks]) * (k + sourceOffset[0][ks] + 1) / 2 + (l + sourceOffset[1][ks])];
        }
      } else if (ks > ls) {
        for (size_t k = 0; k < nk; k++)
          for (size_t l = 0; l < nl; l++) {
            to[k + nk * l] = from[(k + sourceOffset[0][ks]) + nkb * (l + sourceOffset[1][ls])];
          }
      }
    }
  } else if (this->m_parity == parityNone && source.m_parity != parityNone) { // copy triangle to rectangle

    //    molpro::cout << "triangle to rectangle"<<std::endl;
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      size_t nkb = m_dimensions[0][ks];
      //      size_t nlb = m_dimensions[0][ls];
      size_t nk = std::min(m_dimensions[0][ks], source.m_dimensions[0][ks] - sourceOffset[0][ks]);
      size_t nl = std::min(m_dimensions[1][ls], source.m_dimensions[1][ls] - sourceOffset[1][ls]);
      //      molpro::cout << "ks="<<ks<<", ls="<<ls<<std::endl;
      //      if (ks + ls == 0) molpro::cout << "nk=" << nk << ", nl=" << nl << std::endl;
      if (nk * nl == 0)
        continue;
      molpro::array<T> to = this->block(ks);
      molpro::array<T> from = source.block(ks);
      //      molpro::cout << "from.size()="<<from.size()<<std::endl;
      if (ks == ls) { // triangular storage in triangular SMat_
        if (square) {
          //          molpro::cout << "square destination "<<nk<<nl<<nkb<<std::endl;
          if (source.m_parity == parityEven)
            for (size_t k = 0; k < nk; k++)
              for (size_t l = 0; l <= k; l++) {
                //                molpro::cout << k << l << sourceOffset[0][ks] << sourceOffset[1][ls] << std::endl;
                to[k * nkb + l] = to[l * nkb + k] =
                    from[(k + sourceOffset[0][ks]) * ((k + sourceOffset[0][ks]) + 1) / 2 + l + sourceOffset[1][ls]];
              }
          else if (source.m_parity == parityOdd && source.transposed())
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l < k; l++)
                to[l * nkb + k] =
                    -(to[k * nkb + l] = from[(k + sourceOffset[0][ks]) * ((k + sourceOffset[0][ks]) + 1) / 2 + l +
                                             sourceOffset[1][ls]]);
              to[k * nkb + k] = 0;
            }
          else if (source.m_parity == parityOdd)
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l < k; l++)
                to[k * nkb + l] =
                    -(to[l * nkb + k] = from[(k + sourceOffset[0][ks]) * ((k + sourceOffset[0][ks]) + 1) / 2 + l +
                                             sourceOffset[1][ls]]);
              to[k * nkb + k] = 0;
            }
          else if (source.m_parity == parityOddPacked && source.transposed())
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l < k; l++)
                to[l * nkb + k] =
                    -(to[k * nkb + l] = from[(k + sourceOffset[0][ks]) * ((k + sourceOffset[0][ks]) - 1) / 2 + l +
                                             sourceOffset[1][ls]]);
              to[k * nkb + k] = 0;
            }
          else if (source.m_parity == parityOddPacked)
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l < k; l++)
                to[k * nkb + l] =
                    -(to[l * nkb + k] = from[(k + sourceOffset[0][ks]) * ((k + sourceOffset[0][ks]) - 1) / 2 + l +
                                             sourceOffset[1][ls]]);
              to[k * nkb + k] = 0;
            }
        } else { // destination is rectangular
          auto ko = sourceOffset[0][ks];
          auto lo = sourceOffset[1][ls];
          //          molpro::cout << "rectangular destination "<<nk<<nl<<nkb<<ko<<lo<<std::endl;
          //          molpro::cout << "parity "<<m_parity<<source.m_parity<<std::endl;
          //          molpro::cout << "symmetry "<<m_symmetry<<source.m_symmetry<<std::endl;
          if (source.m_parity == parityEven)
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l + lo <= k + ko and l < nl; l++)
                to[k + nkb * l] = from[(k + ko) * (k + ko + 1) / 2 + l + lo];
              for (int l = nl - 1; l + lo > k + ko and l >= 0; l--)
                to[k + nkb * l] = from[(l + lo) * (l + lo + 1) / 2 + k + ko];
            }
          else if (source.m_parity == parityOdd && source.transposed())
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l + lo < k + ko and l < nl; l++)
                to[k + nkb * l] = from[(k + ko) * (k + ko + 1) / 2 + l + lo];
              for (int l = nl - 1; l + lo > k + ko and l >= 0; l--)
                to[k + nkb * l] = -from[(l + lo) * (l + lo + 1) / 2 + k + ko];
              if (k + ko - lo >= 0 and k + ko - lo < nl)
                to[k + nkb * (k + ko - lo)] = 0;
            }
          else if (source.m_parity == parityOdd)
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l + lo < k + ko and l < nl; l++)
                to[k + nkb * l] = -from[(k + ko) * (k + ko + 1) / 2 + l + lo];
              for (int l = nl - 1; l + lo > k + ko and l >= 0; l--)
                to[k + nkb * l] = from[(l + lo) * (l + lo + 1) / 2 + k + ko];
              if (k + ko - lo >= 0 and k + ko - lo < nl)
                to[k + nkb * (k + ko - lo)] = 0;
            }
          else if (source.m_parity == parityOddPacked && source.transposed())
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l + lo < k + ko and l < nl; l++)
                to[k + nkb * l] = from[(k + ko) * (k + ko - 1) / 2 + l + lo];
              for (int l = nl - 1; l + lo > k + ko and l >= 0; l--)
                to[k + nkb * l] = -from[(l + lo) * (l + lo - 1) / 2 + k + ko];
              if (k + ko - lo >= 0 and k + ko - lo < nl)
                to[k + nkb * (k + ko - lo)] = 0;
            }
          else if (source.m_parity == parityOddPacked)
            for (size_t k = 0; k < nk; k++) {
              for (size_t l = 0; l + lo < k + ko and l < nl; l++)
                to[k + nkb * l] = -from[(k + ko) * (k + ko - 1) / 2 + l + lo];
              for (int l = nl - 1; l + lo > k + ko and l >= 0; l--)
                to[k + nkb * l] = from[(l + lo) * (l + lo - 1) / 2 + k + ko];
              if (k + ko - lo >= 0 and k + ko - lo < nl)
                to[k + nkb * (k + ko - lo)] = 0;
            }
        }
      } else { // rectangular storage of unique blocks in triangular SMat_
        bool tr = source.block_transposed(ks);
        int phase = (tr != source.m_transposed && source.m_parity < 0) ? -1 : 1;
        size_t rowstride = tr ? source.m_dimensions[1][ls] : 1;
        size_t colstride = tr ? 1 : source.m_dimensions[0][ks];
        //                molpro::cout << "tr="<<tr<<", nk="<<nk<<", nl="<<nl<<", rowstride="<<rowstride<<",
        //                colstride="<<colstride<<std::endl; molpro::cout << "from=";for (size_t k=0; k<from.size();k++)
        //                molpro::cout <<from[k]<<" "; molpro::cout <<std::endl; molpro::cout << "offsets
        //                "<<sourceOffset[0][ks]<<" "<<sourceOffset[1][ls]<<std::endl;
        if (phase > 0)
          for (size_t l = 0; l < nl; l++)
            for (size_t k = 0; k < nk; k++)
              to[l * nkb + k] = from[(k + sourceOffset[0][ks]) * rowstride + (l + sourceOffset[1][ls]) * colstride];
        else
          for (size_t l = 0; l < nl; l++)
            for (size_t k = 0; k < nk; k++)
              to[l * nkb + k] = -from[(k + sourceOffset[0][ks]) * rowstride + (l + sourceOffset[1][ls]) * colstride];
      }
    }
  } else if (this->m_parity == parityNone && source.m_parity == parityNone) { // copy square to square
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      size_t nkb = m_dimensions[0][ks];
      size_t nk = std::min(m_dimensions[0][ks], source.m_dimensions[0][ks] - sourceOffset[0][ks]);
      size_t nl = std::min(m_dimensions[1][ls], source.m_dimensions[1][ls] - sourceOffset[1][ls]);
      //      if (ks + ls == 0) molpro::cout << "nk=" << nk << ", nl=" << nl << std::endl;
      //      if (ks + ls == 0) molpro::cout << "nkb=" << nkb << ", nlb=" << nlb << std::endl;
      if (nk * nl == 0)
        continue;
      molpro::array<T> to = this->block(ks);
      molpro::array<T> from = source.block(ks);
      bool tr = source.block_transposed(ks);
      size_t rowstride = tr ? source.m_dimensions[0][ls] : 1;
      size_t colstride = tr ? 1 : source.m_dimensions[0][ks];
      //      if (ks + ls == 0) molpro::cout << "rowstride=" << rowstride << ", colstride=" << colstride << std::endl;
      //      if (ks + ls == 0) molpro::cout << "sourceOffset[0][ks]=" << sourceOffset[0][ks] << ",
      //      sourceOffset[1][ls]=" << sourceOffset[1][ls] << std::endl;
      for (size_t l = 0; l < nl; l++)
        for (size_t k = 0; k < nk; k++)
          to[l * nkb + k] = from[(k + sourceOffset[0][ks]) * rowstride + (l + sourceOffset[1][ls]) * colstride];
    }
  } else { // copy square to triangle
           //    molpro::cout << "square to triangle"<<std::endl;
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      molpro::array<T> from1 = source.block(ks);
      molpro::array<T> from2 = source.block(ls);
      size_t nkb = source.m_dimensions[0][ks];
      size_t nlb = source.m_dimensions[0][ls];
      size_t nk = std::min(m_dimensions[0][ks], source.m_dimensions[0][ks] - sourceOffset[0][ks]);
      size_t nl = std::min(m_dimensions[1][ls], source.m_dimensions[1][ls] - sourceOffset[1][ls]);
      molpro::array<T> to = this->block(ks);
      scalar_type fac1 = 0.5;
      scalar_type fac2 = 0.5;
      if (this->parity() < 0) {
        if (source.m_transposed)
          fac2 = -0.5;
        else
          fac1 = -0.5;
      }
      if (ks == ls) {
        if (m_parity == parityOddPacked)
          for (size_t k = 0; k < nk; k++)
            for (size_t l = 0; l < k; l++)
              to[k * (k - 1) / 2 + l] = fac2 * from1[(k + sourceOffset[0][ks]) + nkb * (l + sourceOffset[1][ls])] +
                                        fac1 * from1[(k + sourceOffset[0][ks]) * nlb + (l + sourceOffset[1][ls])];
        else
          for (size_t k = 0; k < nk; k++)
            for (size_t l = 0; l <= k; l++)
              to[k * (k + 1) / 2 + l] = fac2 * from1[(k + sourceOffset[0][ks]) + nkb * (l + sourceOffset[1][ls])] +
                                        fac1 * from1[(k + sourceOffset[0][ks]) * nlb + (l + sourceOffset[1][ls])];
      } else if (ks < ls) {
        for (size_t k = 0; k < nk; k++)
          for (size_t l = 0; l < nl; l++) {
            to[k * nl + l] = fac1 * from1[(k + sourceOffset[0][ks]) + nkb * (l + sourceOffset[1][ls])] +
                             fac2 * from2[(k + sourceOffset[0][ks]) * nlb + (l + sourceOffset[1][ls])];
          }
      }
    }
  }
  return *this;
}

template <class T> void SMat_<T>::splice(SMat_<T> const& source, dims_t sourceOffset, dims_t offset) {
  if (parity() != source.parity())
    throw std::invalid_argument("Non-matching parities");
  if (symmetry() != source.symmetry())
    throw std::invalid_argument("Non-matching symmetries");

  if (rank() == 1 || m_diagonal) {
    size_t smin = 0, smax = max_symmetry_;
    if (!m_diagonal)
      smin = smax = m_symmetry;
    for (unsigned int ks = smin; ks < smax; ks++) {
      size_t nk = m_dimensions[0][ks];
      size_t nkfrom = source.m_dimensions[0][ks];
      if (nk * nkfrom == 0)
        continue;
      molpro::array<T> to = this->block(ks);
      molpro::array<T> from = source.block(ks);
      long rowshift = sourceOffset[0][ks] - offset[0][ks];
      for (size_t k = offset[0][ks]; k < nk && k + rowshift < nkfrom; k++)
        to[k] = from[k + rowshift];
    }
  } else if (parity() == 0) {
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned int ls = ks ^ m_symmetry;
      size_t nk = m_dimensions[0][ks];
      size_t nl = m_dimensions[1][ls];
      size_t nkfrom = source.m_dimensions[0][ks];
      size_t nlfrom = source.m_dimensions[1][ls];
      molpro::array<T> to = this->block(ks);
      molpro::array<T> from = source.block(ks);
      long rowshift = sourceOffset[0][ks] - offset[0][ks];
      long colshift = sourceOffset[1][ls] - offset[1][ls];
      for (size_t l = offset[1][ls]; l < nl && l + colshift < nlfrom; l++)
        for (size_t k = offset[0][ks]; k < nk && k + rowshift < nkfrom; k++)
          to[l * nk + k] = from[(l + colshift) * nkfrom + k + rowshift];
    }
  } else {
    throw std::logic_error("Not implemented yet");
  }
}

template <class T> SMat_<T> SMat_<T>::slice(dims_t dimensions, dims_t offset, std::string description) {
  if (rank() != 2 || m_diagonal)
    throw std::logic_error("slice works only with rank 2 matrices");
  SMat_<T> result(dimensions, (offset[0] != offset[1] || dimensions[0] != dimensions[1]) ? parityNone : this->m_parity,
                  this->m_symmetry, false, description.empty() ? this->m_description : std::move(description));
  result.copy(*this, offset);
  return result;
}

template <class T> SMat_<T>& SMat_<T>::operator*=(T a) {
  //  for (molpro::array<T>::iterator t=(*m_buffer).begin(); t!=(*m_buffer).end(); t++) *t *=a;
  for (auto& t : *m_buffer)
    t *= a;
  return *this;
}

template <class T> const SMat_<T> SMat_<T>::operator*(T a) const {
  SMat_<T> result(this);
  result = *this;
  result *= a;
  return result;
}

template <class T> SMat_<T>& SMat_<T>::operator+=(SMat_<T> const& other) {
  SMat_<T>* o = const_cast<SMat_<T>*>(&other);
  if (m_transposed and !other.m_transposed)
    eval();
  if (!m_transposed and other.m_transposed) {
    o = new SMat_<T>(other);
    *o = other;
    o->eval();
  }
  if (!compatible(*o))
    ErrorExit("Can not add incompatible SMat_ objects");
  for (size_t k = 0; k < size(); k++)
    (*m_buffer)[k] += (*o->m_buffer)[k];
  if (o != &other)
    delete o;
  return *this;
}

template <class T> const SMat_<T> SMat_<T>::operator+(SMat_<T> const& other) const {
  SMat_<T> result(this, m_parity, m_symmetry, m_dimensions.size(), m_diagonal);
  result = *this;
  result += other;
  return result;
}

template <class T> SMat_<T>& SMat_<T>::operator-=(SMat_<T> const& other) {
  if (!compatible(other))
    ErrorExit("Can not subtract incompatible SMat_ objects");
  for (size_t k = 0; k < size(); k++)
    (*m_buffer)[k] -= (*other.m_buffer)[k];
  return *this;
}

template <class T> const SMat_<T> SMat_<T>::operator-(SMat_<T> const& other) const {
  SMat_<T> result(this, m_parity, m_symmetry, m_dimensions.size(), m_diagonal);
  result = *this;
  result -= other;
  return result;
}

template <class T> const SMat_<T> SMat_<T>::operator*(SMat_<T> const& other) const {
  dims_t dims;
  size_t leftd = 0;
  if (m_transposed && rank() == 2)
    leftd = 1;
  if (rank() == 2 || other.rank() == 1)
    dims.push_back(m_dimensions[leftd]);
  size_t rightd = 1;
  if (other.m_transposed || other.rank() == 1)
    rightd = 0;
  if (rank() == 1 || other.rank() == 2)
    dims.push_back(other.m_dimensions[rightd]);
  SMat_<T> result(dims, parityNone, this->symmetry() ^ other.symmetry(), m_diagonal && other.m_diagonal,
                  (!this->m_description.empty() ? this->m_description : "SMat") + "*" +
                      (!other.m_description.empty() ? other.m_description : "SMat"));
  result.multiply(*this, other, 1, 0);
  //  molpro::cout << "constructed product: "<<result.str("",1)<<std::endl;
  return result;
}

template <class T> typename SMat_<T>::scalar_type SMat_<T>::operator&(SMat_<T> const& other) const {
  if (m_diagonal or other.m_diagonal)
    throw std::logic_error("No implementation for diagonal matrices in SMat::operator");
  bool transpose = m_transposed;
  bool othertranspose = other.m_transposed;
  typename SMat_<T>::scalar_type result = 0;
  if (rank() == 1 && other.rank() == 1) {
    if (m_dimensions != other.m_dimensions)
      ErrorExit("Mismatching dimensions, SMat_::operator&");
    for (size_t k = 0; k < size(); k++)
      result += (*m_buffer)[k] * (*other.m_buffer)[k];
  } else if (rank() == 2 && other.rank() == 2) {
    if (rank() != other.rank())
      ErrorExit("Rank mismatch in SMat_::operator&");
    if (symmetry() != other.symmetry())
      ErrorExit("Symmetry mismatch in SMat_::operator&");
    size_t d = 0;
    if (m_transposed != other.m_transposed && rank() == 2)
      d = 1;
    if (m_dimensions[d] != other.m_dimensions[1])
      ErrorExit("Mismatching dimensions, SMat_::operator&");
    if (m_dimensions[1 - d] != other.m_dimensions[0])
      ErrorExit("Mismatching dimensions, SMat_::operator&");
    auto a = const_cast<SMat_<T>*>(this);
    auto b = const_cast<SMat_<T>*>(&other);
    // non-zero parity on input causes a copy
    // it could be done more efficiently
    //      std::cout << this->str("this") << std::endl;
    if (this->parity() != parityNone) {
      transpose = false;
      a = new SMat_<T>(this, parityNone);
      *a = *this;
      //          std::cout << a->str("parity copy of a") << std::endl;
    }
    //      std::cout << other.str("other") << std::endl;
    if (other.parity() != parityNone) {
      othertranspose = false;
      b = new SMat_<T>(&other, parityNone);
      *b = other;
      //        std::cout << other.str("parity source of b") << std::endl;
      //          std::cout << b->str("parity copy of b") << std::endl;
    }

    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      unsigned ls = ks ^ m_symmetry;
      size_t cr = m_dimensions[0][ks];
      size_t cc = m_dimensions[1][ks ^ m_symmetry];
      if (cr * cc == 0)
        continue;
      molpro::array<T> ablock = a->block(m_transposed ? ls : ks);
      molpro::array<T> bblock = b->block(other.m_transposed ? ks : ls);
      for (size_t c = 0; c < cc; c++)
        for (size_t r = 0; r < cr; r++)
          result += ablock[transpose ? r * cc + c : r + cr * c] * bblock[othertranspose ? r + cr * c : r * cc + c];
    }
    if (this->parity() != parityNone)
      delete a;
    if (other.parity() != parityNone)
      delete b;
  } else
    ErrorExit("Rank mismatch in SMat_::operator&");
  return result;
}

template <class T> SMat_<T> SMat_<T>::transform(const SMat_<T>& q, bool orthogonal) const {
  if (orthogonal and m_diagonal)
    return *this;
  if (not orthogonal)
    ErrorExit("SMat::transform not yet coded for non-orthogonal transformation matrix");
  SMat_<T> result(molpro::transpose(q) * (*this) * q); // this is not the optimal algorithm
  return result;
}
template <class T> SMat_<T> SMat_<T>::diagonal() const {
  SMat_<T> result(m_dimensions, parityEven, 0, true, m_description);
  if (m_dimensions.size() != 2 or m_dimensions[0] != m_dimensions[1] or m_symmetry != 0)
    ErrorExit("SMat::diagonal must be called with a square symmetric matrix");
  for (unsigned int block = 0; block < 8; block++) {
    if (m_parity == parityNone)
      for (size_t i = 0; i < m_dimensions[0][block]; i++)
        result.block(block)[i] = blockMap(block)(i, i);
    else if (m_parity == parityEven)
      for (size_t i = 0; i < m_dimensions[0][block]; i++)
        result.block(block)[i] = this->block(block)[(i + 1) * (i + 2) / 2 - 1];
    else
      for (size_t i = 0; i < m_dimensions[0][block]; i++)
        result.block(block)[i] = 0;
  }
  return result;
}

template <class T> void SMat_<T>::orthogonalize(const molpro::SMat_<T>* metric, std::string algorithm) {
  if (algorithm == "Gram-Schmidt") {
    for (unsigned int sym = 0; sym < 8; sym++) {
      //      auto& b = this->blockMap(sym);
      if (metric == nullptr) {
        for (size_t i = 0; i < m_dimensions[1][sym]; i++) {
          for (size_t j = 0; j < i; j++) {
            auto fac = -blockMap(sym).col(j).dot(blockMap(sym).col(i)) / blockMap(sym).col(j).dot(blockMap(sym).col(j));
            blockMap(sym).col(i) += fac * blockMap(sym).col(j);
          }
          blockMap(sym).col(i) /= std::sqrt(blockMap(sym).col(i).dot(blockMap(sym).col(i)));
        }
      } else {
        auto s = metric->blockCopy(sym);
        for (size_t i = 0; i < m_dimensions[1][sym]; i++) {
          for (size_t j = 0; j < i; j++) {
            T fac;
            fac = -blockMap(sym).col(j).dot(s * blockMap(sym).col(i)) /
                  blockMap(sym).col(j).dot(s * blockMap(sym).col(j));
            blockMap(sym).col(i) += fac * blockMap(sym).col(j);
          }
          blockMap(sym).col(i) /= std::sqrt(blockMap(sym).col(i).dot(s * blockMap(sym).col(i)));
        }
      }
    }
  } else
    throw std::runtime_error("Illegal algorithm " + algorithm);
}

template <class T> void SMat_<T>::check_max_symmetry() {
  while (max_symmetry_ > 1) {
    for (auto d = m_dimensions.begin(); d != m_dimensions.end(); d++)
      if (*std::max_element(d->begin() + max_symmetry_ / 2, d->end()) > 0)
        return;
    max_symmetry_ /= 2;
  }
}

template <class T> void SMat_<T>::multiply(const SMat_<T>& aa, const SMat_<T>& bb, T alpha, T beta) {
  check_max_symmetry();
  if (size() == 0)
    return;
  if ((rank() != 2 && aa.rank() == bb.rank()) || (rank() != 1 && aa.rank() != bb.rank()))
    ErrorExit("wrong rank SMat_::multiply");
  if (m_parity != parityNone)
    ErrorExit("wrong result parity SMat_::multiply");
  if (aa.m_symmetry < 0 || aa.m_symmetry > 7)
    ErrorExit("aa.symmetry problem");
  if (bb.m_symmetry < 0 || bb.m_symmetry > 7)
    ErrorExit("bb.symmetry problem");
  if (m_symmetry != (aa.m_symmetry ^ bb.m_symmetry))
    ErrorExit("wrong result symmetry, SMat_::multiply");
  auto a = const_cast<SMat_<T>*>(&aa);
  auto b = const_cast<SMat_<T>*>(&bb);
  // non-zero parity on input causes a copy
  if (aa.parity() != parityNone) {
    a = new SMat_<T>(&aa, parityNone);
    *a = aa;
  }
  if (bb.parity() != parityNone) {
    b = new SMat_<T>(&bb, parityNone);
    *b = bb;
  }
  int transposea = (a->m_transposed) ? 1 : 0;
  int transposeb = (b->m_transposed) ? 1 : 0;

  if (a->rank() == 2 && b->rank() == 2) { // matrix * matrix
    if (a->m_dimensions[1 - transposea] != b->m_dimensions[transposeb])
      ErrorExit("Mismatch in a,b dimensions in SMat_::multiply");
    if (m_dimensions[0] != a->m_dimensions[transposea])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    if (m_dimensions[1] != b->m_dimensions[1 - transposeb])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      auto ls = ks ^ a->m_symmetry;
      auto cr = m_dimensions[0][ks];
      auto cc = m_dimensions[1][ks ^ m_symmetry];
      if (cr * cc == 0)
        continue;
      if (a->m_diagonal and b->m_diagonal) {
        if (!m_diagonal)
          throw std::logic_error("result of multiplying two diagonals must be diagonal");
        if (beta == T{0}) {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              block(ks)[i] = alpha * (a->block(ks)[i] * b->block(ks)[i]);
        } else if (beta == T{1}) {
          for (size_t i = 0; i < cr; i++)
            block(ks)[i] += alpha * (a->block(ks)[i] * b->block(ks)[i]);
        } else {
          for (size_t i = 0; i < cr; i++)
            block(ks)[i] = beta * block(ks)[i] + alpha * (a->block(ks)[i] * b->block(ks)[i]);
        }

      } else if (a->m_diagonal) {
        if (beta == T{0}) {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) = alpha * (b->blockMap(ks)(i, j) * a->block(ks)[i]);
        } else if (beta == T{1}) {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) += alpha * (b->blockMap(ks)(i, j) * a->block(ks)[i]);
        } else {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) = beta * blockMap(ks)(i, j) + alpha * (b->blockMap(ks)(i, j) * a->block(ks)[i]);
        }
      } else if (b->m_diagonal) {
        if (beta == T{0}) {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) = alpha * (a->blockMap(ks)(i, j) * b->block(ls)[j]);
        } else if (beta == T{1}) {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) += alpha * (a->blockMap(ks)(i, j) * b->block(ls)[j]);
        } else {
          for (size_t i = 0; i < cr; i++)
            for (size_t j = 0; j < cc; j++)
              blockMap(ks)(i, j) = beta * blockMap(ks)(i, j) + alpha * (a->blockMap(ks)(i, j) * b->block(ls)[j]);
        }
      } else {
#ifdef SMAT_USE_BLAS
        char transa = transposea != 0 ? 't' : 'n';
        char transb = transposeb != 0 ? 't' : 'n';
        unsigned int lst = ls;
        if (transposeb)
          lst = ls ^ b->m_symmetry;
        auto k = a->m_dimensions[1 - transposea][ls];
        auto lda = a->m_dimensions[0][transposea ? ls : ks];
        auto ldb = b->m_dimensions[0][lst];
        auto cblock = block(ks);
        auto ablock = a->block(transposea ? ls : ks);
        auto bblock = b->block(lst);
        if (sMat_debug) {
          molpro::cout << "before dgemm ablock=";
          for (auto k = ablock.begin(); k < ablock.end(); k++)
            molpro::cout << " " << *k;
          molpro::cout << std::endl;
          molpro::cout << "before dgemm bblock=";
          for (auto k = bblock.begin(); k < bblock.end(); k++)
            molpro::cout << " " << *k;
          molpro::cout << std::endl;
          molpro::cout << "transa=" << transa << std::endl;
          molpro::cout << "transb=" << transb << std::endl;
          molpro::cout << "cr=" << cr << std::endl;
          molpro::cout << "cc=" << cc << std::endl;
          molpro::cout << "k=" << k << std::endl;
          molpro::cout << "lda=" << lda << std::endl;
          molpro::cout << "ldb=" << ldb << std::endl;
        }
        dgemm(transa, transb, cr, cc, k, alpha, &ablock[0], lda, &bblock[0], ldb, beta, &cblock[0], cr);
#else // Eigen
        if (beta == T{0})
          blockMap(ks).noalias() = alpha * (a->blockMap(ks) * b->blockMap(ls));
        else if (beta == T{1})
          blockMap(ks).noalias() += alpha * (a->blockMap(ks) * b->blockMap(ls));
        else
          blockMap(ks).noalias() = beta * blockMap(ks) + alpha * (a->blockMap(ks) * b->blockMap(ls));
#endif
      }
    }
  }

  if (a->rank() == 1 && b->rank() == 1) { // vector * vector
    if (m_dimensions[0] != a->m_dimensions[0])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    if (m_dimensions[1] != b->m_dimensions[0])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    auto ablock = a->block(a->m_symmetry);
    auto bblock = b->block(b->m_symmetry);
    auto cblock = block(a->m_symmetry);
    auto cr = m_dimensions[0][a->m_symmetry];
    auto cc = m_dimensions[1][b->m_symmetry];
    for (size_t c = 0; c < cc; c++)
      for (size_t r = 0; r < cr; r++)
        if (beta == T{0})
          cblock[r + cr * c] = alpha * ablock[r] * bblock[c];
        else
          cblock[r + cr * c] = beta * cblock[r + cr * c] + alpha * ablock[r] * bblock[c];
  }

  if (a->rank() == 2 && b->rank() == 1) { // matrix * vector
    if (m_dimensions[0] != a->m_dimensions[transposea])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    if (a->m_dimensions[1 - transposea] != b->m_dimensions[0])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
#ifdef SMAT_USE_BLAS
    char transa = transposea != 0 ? 't' : 'n';
    auto ablock = a->block(transposea ? b->m_symmetry : m_symmetry);
    auto bblock = b->block(b->m_symmetry);
    auto cblock = block(m_symmetry);
    auto cr = a->m_dimensions[0][transposea ? b->m_symmetry : m_symmetry];
    auto cc = a->m_dimensions[1][transposea ? m_symmetry : b->m_symmetry];
    dgemv(transa, cr, cc, alpha, &ablock[0], cr, &bblock[0], 1, beta, &cblock[0], 1);
#else
    if (beta == T{0})
      blockV(m_symmetry).noalias() = alpha * a->blockMap(m_symmetry) * b->blockV(b->m_symmetry);
    else
      blockV(m_symmetry).noalias() =
          beta * blockV(m_symmetry) + alpha * a->blockMap(m_symmetry) * b->blockV(b->m_symmetry);
#endif
  }

  if (a->rank() == 1 && b->rank() == 2) { // vector * matrix
    if (m_dimensions[0] != b->m_dimensions[1 - transposeb])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
    if (a->m_dimensions[0] != b->m_dimensions[transposeb])
      ErrorExit("mismatch in result dimensions in SMat_::multiply");
#ifdef SMAT_USE_BLAS
    auto ablock = a->block(a->m_symmetry);
    auto bblock = b->block(transposeb ? m_symmetry : a->m_symmetry);
    auto cblock = block(m_symmetry);
    auto cr = b->m_dimensions[0][transposeb ? m_symmetry : a->m_symmetry];
    auto cc = b->m_dimensions[1][transposeb ? a->m_symmetry : m_symmetry];
    dgemv((transposeb != 0 ? 'n' : 't'), cr, cc, alpha, &bblock[0], cr, &ablock[0], 1, beta, &cblock[0], 1);
#else
    if (beta == T{0})
      blockV(m_symmetry).noalias() = alpha * a->blockV(a->m_symmetry).transpose() * b->blockMap(a->m_symmetry);
    else
      blockV(m_symmetry).noalias() =
          beta * blockV(m_symmetry) + alpha * a->blockV(a->m_symmetry).transpose() * b->blockMap(a->m_symmetry);
#endif
  }

  if (aa.parity() != parityNone)
    delete a;
  if (bb.parity() != parityNone)
    delete b;
  if (sMat_debug) {
    molpro::cout << this->str("this at end of multiply", 6);
    memory_print_status();
    molpro::cout << "at end of multiply before return result" << std::endl << std::flush;
  }
}

template <class T> void SMat_<T>::transpose() { m_transposed = not m_transposed; }

template <class T> void SMat_<T>::eval() { // FIXME incomplete implementation, incompleteness not trapped
  if (m_transposed) {
    if (rank() != 2)
      throw std::logic_error("unexpected rank");
    if (m_symmetry != 0 || m_parity != parityNone)
      throw std::logic_error("unimplemented");
    for (unsigned int sym = 0; sym < max_symmetry_; sym++) {
      auto b = blockCopy(sym);
      std::swap(m_dimensions[0][sym], m_dimensions[1][sym]);
      blockMap(sym) = b.transpose();
    }
    m_transposed = false;
  }
}

template <class T> void SMat_<T>::trim(typename SMat_<T>::scalar_type cut) {
  auto cut_ = std::abs(std::abs(cut) > 0 ? cut : Eigen::NumTraits<typename SMat_<T>::scalar_type>::epsilon());
  std::replace_if(
      m_buffer->begin(), m_buffer->end(), [cut_](T x) { return std::abs(x) < cut_; }, T{0});
}

template <class T> molpro::array<T>* SMat_<T>::data() const { return m_buffer; }

template <class T> void SMat_<T>::scal(T a, bool scaleDiagonal) {
  if (std::abs(a) == 0) {
    assign(0);
    return;
  }
  if (scaleDiagonal || m_parity == parityOdd || m_parity == parityOddPacked)
    for (size_t i = 0; i < size(); i++)
      (*m_buffer)[i] *= a;
  else { // only the off-diagonal elements are scaled
    if (rank() < 2)
      throw std::logic_error("Scaling with fixed diagonals does not make sense for a vector");
    if (m_parity == parityNone)
      throw std::logic_error("Scaling with fixed diagonals is not implemented for parity=0");
    if (m_symmetry != 0)
      throw std::logic_error("Scaling with fixed diagonals does not make sense for a non-symmetric matrix");
    auto p = &((*m_buffer)[0]);
    for (size_t sym = 0; sym < 8; sym++)
      for (size_t i = 0; i < m_dimensions[0][sym]; i++) {
        for (size_t j = 0; j < i; j++)
          *(p++) *= a;
        p++;
      }
  }
}

template <class T> void SMat_<T>::axpy(T a, const SMat_<T>& x) {
  if (!compatible(x))
    ErrorExit("Can not add incompatible SMat_ objects");
  for (size_t i = 0; i < size(); i++)
    (*m_buffer)[i] += a * (*x.m_buffer)[i];
}

template <class T>
SMat_<T>::SVD::SVD(const SMat_<T>& matrix, std::string algorithm, unsigned int computationOptions) : m_matrix(matrix) {
  if (m_matrix.symmetry() != 0)
    throw std::runtime_error("Cannot SVD non-symmetric SMat_");
  if (m_matrix.m_diagonal)
    throw std::runtime_error("Cannot SVD diagonal SMat_");
  if (m_matrix.parity() != parityNone)
    throw std::runtime_error("Cannot SVD a triangular-packed SMat_");
  dim_t nn, mm, pp;
  for (auto s = 0; s < 8; s++) {
    const auto& n = m_matrix.dimension(s, 0);
    const auto& p = m_matrix.dimension(s, 1);
    nn.push_back(n);
    pp.push_back(p);
    mm.emplace_back(std::min(n, p));
  }
  m_U = std::make_unique<SMat_<T>>((computationOptions & Eigen::ComputeFullU) ? dims_t{nn, nn} : dims_t{nn, mm});
  m_V = std::make_unique<SMat_<T>>((computationOptions & Eigen::ComputeFullU) ? dims_t{pp, pp} : dims_t{pp, mm});
  m_S = std::make_unique<SMat_<T>>(dims_t{mm, mm}, parityNone, 0, true);
  for (auto s = 0; s < 8; s++) {
    if (mm[s] == 0)
      continue;
    if (algorithm == "BDC") {
      Eigen::BDCSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(m_matrix.blockMap(s), computationOptions);
      const auto& sv = svd.singularValues();
      for (Eigen::Index i = 0; i < sv.size(); i++)
        m_S->block(s)[i] = sv[i];
      m_U->blockMap(s) = svd.matrixU();
      m_V->blockMap(s) = svd.matrixV();
    } else if (algorithm == "Jacobi") {
      Eigen::JacobiSVD<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> svd(m_matrix.blockMap(s), computationOptions);
      const auto& sv = svd.singularValues();
      for (Eigen::Index i = 0; i < sv.size(); i++)
        m_S->block(s)[i] = sv[i];
      m_U->blockMap(s) = svd.matrixU();
      m_V->blockMap(s) = svd.matrixV();
    } else
      throw std::invalid_argument("Illegal algorithm \"" + algorithm + "\" in SMat_::SVD");
  }
  m_S->m_description = matrix.m_description + " SVD S";
  m_U->m_description = matrix.m_description + " SVD U";
  m_V->m_description = matrix.m_description + " SVD V";
}

template <class T> SMat_<T> SMat_<T>::desymmetrise() const {
  dims_t newdim;
  dim_t dim(8, 0);
  for (size_t i = 0; i < rank(); i++) {
    newdim.push_back(dim);
    for (int j = 0; j < 8; j++)
      newdim[i][0] += this->m_dimensions[i][j];
  }
  SMat_<T> newmat(newdim, this->m_parity, 0, false, this->m_description);
  newmat.assign(0);
  T* to = &newmat.block(0)[0];
  size_t nrnew = newdim[0][0];
  if (m_diagonal)
    for (size_t i = 0; i < newmat.m_dimensions[0][0]; i++)
      newmat.block(0)[i] = this->block(0)[i];
  else if (rank() == 1)
    for (size_t i = 0; i < m_dimensions[0][m_symmetry]; i++)
      newmat.block(0)[i + block_offset(m_symmetry)] = this->block(m_symmetry)[i];
  else { // rank==2
    if (m_parity == 0 || m_symmetry != 0) {
      for (int sr = 0; sr < 8; sr++) {
        size_t nr = m_dimensions[0][sr];
        int sc = sr ^ m_symmetry;
        if (sc > sr)
          continue;
        size_t nc = m_dimensions[0][sc];
        for (size_t c = 0; c < nc; c++)
          for (size_t r = 0; r < nr; r++)
            to[r + c * nrnew] = block(sr)[r + c * nr];
        to += nrnew * nc + nr;
      }
    } else if (m_parity == parityOddPacked) { // symmetry==0
      size_t cabs = 0;
      for (int sr = 0; sr < 8; sr++) {
        size_t nr = m_dimensions[0][sr];
        for (size_t c = 0; c < nr; c++) {
          for (size_t r = 0; r < c; r++)
            to[r] = block(sr)[c * (c - 1) / 2 + r];
          to += cabs;
          cabs++;
        }
        to += nr;
      }
    } else { // symmetry==0 && parity !=0
      size_t cabs = 0;
      for (int sr = 0; sr < 8; sr++) {
        size_t nr = m_dimensions[0][sr];
        for (size_t c = 0; c < nr; c++) {
          for (size_t r = 0; r <= c; r++)
            to[r] = block(sr)[c * (c + 1) / 2 + r];
          cabs++;
          to += cabs;
        }
        to += nr;
      }
    }
  }
  return SMat_<T>(*this);
}

template <class T> bool SMat_<T>::operator==(SMat_<T> const& other) const {
  const T* thisbuff = this->m_buffer->data();
  const T* otherbuff = other.m_buffer->data();
  std::unique_ptr<SMat_<T>> thiscopy, othercopy;
  if (rank() != other.rank())
    return false;
  for (size_t k = 0; k < m_dimensions.size(); k++)
    if (m_dimensions[k] != other.m_dimensions[k])
      return false;
  if (symmetry() != other.symmetry())
    return false;
  //  molpro::cout << "this:\n"<<*this<<std::endl;
  //  molpro::cout << "other:\n"<<other<<std::endl;
  if (parity() != other.parity() || m_transposed != other.m_transposed or
      m_diagonal != other.m_diagonal) { // implement cross-parity checking by copying
    if (parity() != parityNone || m_transposed or (rank() == 2 and m_diagonal and not other.m_diagonal)) {
      //      molpro::cout << "copy this " << std::endl;
      thiscopy.reset(new SMat_<T>(this, parityNone));
      *thiscopy = *this;
      thisbuff = thiscopy->m_buffer->data();
      //      molpro::cout << "thiscopy\n" << *thiscopy << std::endl;
    }
    //    molpro::cout << "rank() "<<rank()<<" m_diagonal "<<m_diagonal<<" other.m_diagonal
    //    "<<other.m_diagonal<<std::endl;
    if (other.parity() != parityNone || other.m_transposed or (rank() == 2 and not m_diagonal and other.m_diagonal)) {
      //      molpro::cout << "copy other " << std::endl;
      othercopy.reset(new SMat_<T>(&other, parityNone));
      *othercopy = other;
      otherbuff = othercopy->m_buffer->data();
      //      molpro::cout << "othercopy\n" << *othercopy << std::endl;
    }
  }
  //  molpro::cout << "thisbuff, otherbuff:\n";
  for (size_t k = 0; k < size(); k++)
    //    molpro::cout << thisbuff[k] << " : " << otherbuff[k]
    //         << " : " << (thisbuff[k] == otherbuff[k]) << std::endl;
    for (size_t k = 0; k < size(); k++)
      if (thisbuff[k] != otherbuff[k])
        return false; // could implement tolerant checking
  //      if ((*m_buffer)[k] != (*other.m_buffer)[k]) return false; // could implement tolerant checking
  return true;
}

template <class T> bool SMat_<T>::operator!=(SMat_<T> const& other) const { return !(*this == other); }

template <class T> SMat_<T> SMat_<T>::exp() const {
  if (m_symmetry != 0)
    ErrorExit("SMat::exp has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::exp has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::exp has received a non-square matrix");
  SMat_<T> result(m_dimensions, parityNone, 0);
  if (m_diagonal) {
    throw std::logic_error("not implemented");
  } else {
    const SMat_<T>* source = this;
    if (parity())
      *(const_cast<SMat_<T>*>((source = new SMat_<T>(m_dimensions, parityNone, 0)))) = *this;
    for (size_t sym = 0; sym < 8; sym++)
      if (m_dimensions[0][sym] > 0)
        result.blockMap(sym) =
            (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(source->blockMap(sym))).exp();
    if (parity())
      delete source;
    return result;
  }
}

template <class T> SMat_<T> SMat_<T>::log() const {
  if (m_symmetry != 0)
    ErrorExit("SMat::log has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::log has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::log has received a non-square matrix");
  SMat_<T> result(m_dimensions, parityNone, 0);
  const SMat_<T>* source = this;
  if (parity())
    *(const_cast<SMat_<T>*>((source = new SMat_<T>(m_dimensions, parityNone, 0)))) = *this;
  for (size_t sym = 0; sym < 8; sym++)
    if (m_dimensions[0][sym] > 0)
      result.blockMap(sym) =
          (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(source->blockMap(sym))).log();
  if (parity())
    delete source;
  return result;
}

template <class T> SMat_<T> SMat_<T>::sqrt() const {
  if (m_symmetry != 0)
    ErrorExit("SMat::sqrt has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::sqrt has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::sqrt has received a non-square matrix");
  SMat_<T> result(m_dimensions, parityNone, 0);
  const SMat_<T>* source = this;
  if (parity())
    *(const_cast<SMat_<T>*>((source = new SMat_<T>(m_dimensions, parityNone, 0)))) = *this;
  for (size_t sym = 0; sym < 8; sym++)
    if (m_dimensions[0][sym] > 0)
      result.blockMap(sym) =
          (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(source->blockMap(sym))).sqrt();
  if (parity())
    delete source;
  return result;
}

template <class T> SMat_<T> SMat_<T>::pow(T p) const {
  if (m_symmetry != 0)
    ErrorExit("SMat::pow has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::pow has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::pow has received a non-square matrix");
  SMat_<T> result(m_dimensions, parityNone, 0);
  const SMat_<T>* source = this;
  if (parity())
    *(const_cast<SMat_<T>*>((source = new SMat_<T>(m_dimensions, parityNone, 0)))) = *this;
  for (size_t sym = 0; sym < 8; sym++)
    if (m_dimensions[0][sym] > 0)
      result.blockMap(sym) =
          (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(source->blockMap(sym))).pow(p);
  if (parity())
    delete source;
  return result;
}

template <class T> SMat_<T> SMat_<T>::inverse(typename SMat_<T>::scalar_type SVThresh) const {
  if (m_symmetry != 0)
    ErrorExit("SMat::inverse has received a matrix that is not of the totally symmetric irrep");
  if (rank() != 2)
    ErrorExit("SMat::inverse has not received a matrix");
  if (m_dimensions[0] != m_dimensions[1])
    ErrorExit("SMat::inverse has received a non-square matrix");
  SMat_<T> result(m_dimensions, parityNone, 0);
  const SMat_<T>* source = this;
  if (parity())
    *(const_cast<SMat_<T>*>((source = new SMat_<T>(m_dimensions, parityNone, 0)))) = *this;
  if (SVThresh != 0.0)
    result.assign(0);
  for (size_t sym = 0; sym < 8; sym++)
    if (m_dimensions[0][sym] > 0) {
      if (SVThresh == 0.0)
        result.blockMap(sym) =
            (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(source->blockMap(sym))).inverse();
      else {                                                                           // pseudo inverse
        throw std::logic_error("Pseudo-inverse unavailable because of compiler bugs"); // TODO enable the following code
                                                                                       //        auto svd =
        //            (static_cast<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(source->blockMap(sym))).jacobiSvd(
        //                Eigen::ComputeFullU | Eigen::ComputeFullV);
        //        const auto& sv = svd.singularValues();
        //        for (Eigen::Index i = 0; i < sv.size(); i++)
        //          if (std::abs(sv(i)) > std::abs(SVThresh))
        //            result.blockMap(sym) += svd.matrixV().col(i) * (1 / sv(i)) * svd.matrixU().col(i).adjoint();
      }
    }
  if (parity())
    delete source;
  return result;
}
} // namespace molpro
