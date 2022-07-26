#ifndef SYMMETRY_MATRIX_SMATMAT_IMPLEMENTATION_H
#define SYMMETRY_MATRIX_SMATMAT_IMPLEMENTATION_H
#include <algorithm>
#include <cstdint>
#include "SMatMat.h"
#include "memory.h"
#include <sstream>

using namespace molpro;

template <class T>
SMatMat_<T>::SMatMat_(const T& matrix, dims_t dimensions, typename T::value_type* buffer, parity_t parity, int symmetry,
                      std::string description)
    : m_dimensions(std::move(dimensions)), m_parity(parity), m_symmetry(symmetry),
      m_description(std::move(description)) {
  for (size_t sym = 0; sym < 8; sym++) {
    m_templates.push_back(std::make_shared<T>(matrix.dimensions(), &m_zero_buffer[0], matrix.parity(), sym));
    //      std::cout << "template parity="<<m_templates.back()->parity()<<",
    //      size="<<m_templates.back()->size()<<std::endl;
  }
  for (size_t axis = 0; axis < rank(); axis++)
    while (this->m_dimensions[axis].size() < 8)
      this->m_dimensions[axis].push_back(0);
  if ((m_managed_buffer = (buffer == nullptr))) { // no buffer provided, so store internally
    size_t n;                                     // For dumb pgc++ - bug 5086
    this->m_bufferp = std::make_shared<molpro::array<typename T::value_type>>(n = size());
  } else {    // attach to provided buffer
    size_t n; // For dumb pgc++ - bug 5086
    this->m_bufferp =
        std::make_shared<molpro::array<typename T::value_type>>(buffer, n = size()); // no checks on size of buffer!
  }
  this->m_buffer = &this->m_bufferp.get()[0];
  //  if (parity < 0) throw std::logic_error("Negative parity SMatMat not yet supported");
}

template <class T> SMatMat_<T>::SMatMat_(const char* dump, typename T::value_type* buffer) {
  //  std::cout << "construct SMatMat from bytestream dump"<<std::endl;
  //  int64_t l; std::memcpy(&l,dump,8); printf("%lX %lld\n",l,l);
  class molpro::bytestream bs(dump);
  m_symmetry = bs.ints()[0];
  m_parity = static_cast<parity_t>(bs.ints()[0]);
  for (int axis = 0; axis < 2; axis++) {
    auto ii = bs.ints();
    m_dimensions.push_back(dim_t(ii.size()));
    for (size_t k = 0; k < ii.size(); k++)
      m_dimensions.back()[k] = ii[k];
    //      for (size_t k=0; k<ii.size();k++) std::cout << "dimension "<<m_dimensions.back()[k]<<std::endl;
  }
  auto dd = bs.chars();
  std::copy(dd.begin(), dd.end(), m_description.begin());
  m_managed_buffer = true;
  for (unsigned int sym = 0; sym < 8; sym++) {
    auto bss = bs.byteStream();
    //      std::cout <<"template bytestream" <<std::endl; bss.dump();
    SMat matrix(bss);
    //      std::cout << matrix << std::endl;
    m_templates.push_back(std::make_shared<T>(matrix.dimensions(), &m_zero_buffer[0], matrix.parity(), sym));
  }
  size_t n; // For dumb pgc++ - bug 5086
  if (buffer == nullptr)
    this->m_bufferp =
        std::make_shared<molpro::array<typename T::value_type>>(n = size()); // no checks on size of buffer!
  else
    this->m_bufferp =
        std::make_shared<molpro::array<typename T::value_type>>(buffer, n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
  if (bs.size() > bs.position()) {
    auto data = bs.doubles();
    //      std::cout << "recovering floating point data "; for (auto dd=data.begin(); dd!=data.end();dd++) std::cout <<
    //      " "<<*dd; std::cout << std::endl;
    std::copy(data.begin(), data.end(), this->m_buffer->begin());
  }
}

template <class T>
SMatMat_<T>::SMatMat_(const SMatMat_<T>& source)
    : SMatMat_(*source.m_templates[0], source.m_dimensions, nullptr, source.m_parity, source.m_symmetry,
               source.m_description) {
  *this = source;
}

template <class T> SMatMat_<T>& SMatMat_<T>::operator=(const SMatMat_<T>& other) {
  checkCompatible(other);
  // pointer copy
  //  *m_buffer = *other.m_buffer;
  // deep copy
  std::copy(other.m_buffer->begin(), other.m_buffer->end(), m_buffer->begin());
  return *this;
}

template <class T> SMatMat_<T>& SMatMat_<T>::operator+=(const SMatMat_<T>& other) {
  checkCompatible(other);
  std::transform(m_buffer->begin(), m_buffer->end(), other.m_buffer->begin(), m_buffer->begin(),
                 std::plus<typename T::value_type>());
  return *this;
}

template <class T> SMatMat_<T>& SMatMat_<T>::operator-=(const SMatMat_<T>& other) {
  checkCompatible(other);
  std::transform(m_buffer->begin(), m_buffer->end(), other.m_buffer->begin(), m_buffer->begin(),
                 std::minus<typename T::value_type>());
  return *this;
}

template <class T> SMatMat_<T>& SMatMat_<T>::operator*=(typename T::value_type other) {
  for (auto& v : *m_buffer)
    v *= other;
  return *this;
}

template <class T> typename T::scalar_type SMatMat_<T>::operator&(const SMatMat_<T>& other) const {
  typename T::value_type result = 0;
  if (rank() < 2)
    throw std::logic_error("SMatMat::operator& works only with matrices not vectors");
  for (unsigned int ijsym = 0; ijsym < max_symmetry_; ijsym++)
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      size_t nr = m_dimensions[0][ks], nc = m_dimensions[1][ks ^ ijsym];
      if (nr < 1 || nc < 1)
        continue;
      //      if (rank() == 2 && m_parity != parityNone && ks < (ks^ijsym)) continue;
      //      s<<"@@ ijsym="<<ijsym<<", ks="<<ks<<std::endl;
      //      s<<"@@ nr="<<nr<<", nc="<<nc<<std::endl;
      for (size_t k = 0; k < nr; k++)
        //          if (m_parity != parityNone && ijsym==0)
        //            for (size_t l=0; l<=k; l++)
        //          else
        for (size_t l = 0; l < nc; l++) {
          result += (*smat(ijsym, ks, k, l) & *other.smat(ijsym, ks, k, l));
          //              std::cout << "matrix1:\n"<<*smat(ijsym,ks,k,l)<<std::endl;
          //              std::cout << "matrix2:\n"<<*other.smat(ijsym,ks,k,l)<<std::endl;
        }
    }
  return result;
}

template <class T> std::shared_ptr<T> SMatMat_<T>::smat(unsigned int ijsym, unsigned int isym, int i, int j) const {
  unsigned int jsym = ijsym ^ isym;
  size_t address;
  bool transpose = false;
  if (m_parity != parityNone && (jsym > isym)) { // accessing transpose, assume abij=baji symmetry
    transpose = true;
    address = block_offset(ijsym, jsym) + (j + i * m_dimensions[0][jsym]) * m_templates[ijsym ^ m_symmetry]->size();
  } else if (m_parity != parityNone && m_parity != parityOddPacked && (0 == ijsym)) {
    transpose = i < j;
    address = block_offset(ijsym, isym) +
              ((i > j) ? (i * (i + 1) / 2 + j) : (j * (j + 1) / 2 + i)) * m_templates[ijsym ^ m_symmetry]->size();
  } else if (m_parity == parityOddPacked && (0 == ijsym)) {
    transpose = i < j;
    address = block_offset(ijsym, isym) +
              ((i > j) ? (i * (i - 1) / 2 + j) : (j * (j - 1) / 2 + i)) * m_templates[ijsym ^ m_symmetry]->size();
  } else {
    address = block_offset(ijsym, isym) + (i + j * m_dimensions[0][isym]) * m_templates[ijsym ^ m_symmetry]->size();
  }
  //      std::cout << "SMatMat::smat "<<ijsym<<isym<<i<<j<<" address="<<address<<std::endl;
  auto result = std::make_shared<T>(&(*m_templates[ijsym ^ m_symmetry]), &(*m_buffer)[address],
                                    m_templates[ijsym ^ m_symmetry]->parity());
  if (transpose)
    result->transpose();
  //  std::cout <<"SMatMat::smat isym="<<isym<<" jsym="<<jsym<<" transpose="<<transpose<<std::endl;
  return result;
}

template <class T> size_t SMatMat_<T>::size() const { return block_offset(8); }

template <class T> unsigned int SMatMat_<T>::symmetry() const { return m_symmetry; }

template <class T> size_t SMatMat_<T>::block_size(unsigned int ijsym, unsigned int isym) const {
  if (rank() == 1)
    throw std::logic_error("SMatMat::block_size can not be called with two arguments for rank-1 structure");
  if (m_parity != parityNone && ijsym == 0) {
    //      std::cout << "SMatMat::block_size("<<ijsym<<","<<isym<<") " <<m_dimensions[0][isym]*(m_dimensions[0][isym]+
    //      ((m_parity==parityOddPacked)?-1:1))/2<<"; parity="<<m_parity<<std::endl;
    return m_dimensions[0][isym] * (m_dimensions[0][isym] + ((m_parity == parityOddPacked) ? -1 : 1)) / 2;
  } else
    return m_dimensions[0][isym] * m_dimensions[1][ijsym ^ isym];
}

template <class T> size_t SMatMat_<T>::block_size(unsigned int ijsym) const {
  if (rank() == 1) {
    throw std::logic_error("SMatMat::block_size not yet implemented for rank-1 structure");
  } else {
    size_t result = 0;
    for (size_t isym = 0; isym < 8; isym++)
      if (parity() == parityNone || isym >= (ijsym ^ isym))
        result += block_size(ijsym, isym);
    //      std::cout << "SMatMat::block_size("<<ijsym<<") = " <<result<<std::endl;
    return result;
  }
}

template <class T> size_t SMatMat_<T>::block_offset(unsigned int ijsym, unsigned int isym) const {
  size_t result = 0;
  if (rank() == 2) {
    unsigned int jsym = ijsym ^ isym;
    if (jsym > isym && parity() != parityNone)
      isym = jsym;
    for (unsigned int ijsym1 = 0; ijsym1 < 8; ijsym1++) {
      for (unsigned int isym1 = 0; isym1 < 8; isym1++) {
        if (ijsym1 == ijsym && isym1 == isym)
          return result;
        if (parity() == parityNone || isym1 >= (isym1 ^ ijsym1))
          result += block_size(ijsym1, isym1) * m_templates[ijsym1 ^ m_symmetry]->size();
      }
    }
  } else if (rank() == 1) {
    throw std::logic_error("SMatMat::block_offset not yet implemented for rank=1");
  }
  return result;
}

template <class T> size_t SMatMat_<T>::block_offset(unsigned int ijsym) const {
  size_t result = UINT_MAX;
  for (unsigned int isym = 0; isym < 8; isym++)
    result = std::min(result, block_offset(ijsym, isym));
  return result;
}

template <class T> std::string SMatMat_<T>::str(std::string title, int level) const {
  if (level < 0)
    return "";
  if (title.empty())
    title = m_description;
  std::stringstream s;
  s << "Supermatrix " << (title.empty() ? m_description : title);
  if (m_symmetry >= 0)
    s << "; symmetry=" << m_symmetry + 1;
  if (rank() > 1)
    s << " parity=" << m_parity;
  if (level > 0) {
    if (m_buffer != nullptr)
      s << (m_managed_buffer ? "; managed" : " unmanaged") << " data at address " << &(*m_buffer)[0];
    else
      s << "; no data attached";
  }
  s << std::endl;
  if (level > 1 && m_buffer != nullptr) {
    s << "raw data:";
    for (size_t k = 0; k < size(); k++)
      s << " " << (*m_buffer)[k];
    s << std::endl;
  }
  for (unsigned int ijsym = 0; ijsym < (rank() > 1 ? max_symmetry_ : 1); ijsym++)
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      size_t nr = m_dimensions[0][ks], nc = (rank() > 1 ? m_dimensions[1][ks ^ ijsym] : 1);
      if (nr < 1 || nc < 1)
        continue;
      if (rank() == 2 && m_parity != parityNone && ks < (ks ^ ijsym))
        continue;
      //      s<<"@@ ijsym="<<ijsym<<", ks="<<ks<<std::endl;
      //      s<<"@@ nr="<<nr<<", nc="<<nc<<std::endl;
      if (rank() == 1)
        s << "Block (" << ks + 1 << "), dimensions (" << nr << ")" << std::endl;
      else
        s << "Block (" << ks + 1 << "," << (ks ^ ijsym) + 1 << "), dimensions (" << nr << "," << nc << ")" << std::endl;
      if (m_buffer == nullptr || level < 0)
        continue;
      s << std::endl;
      if (rank() == 1)
        for (size_t k = 0; k < nr; k++)
          s << "("
            << smat(ijsym, ks, k)
                   ->str(title + " " + std::to_string(static_cast<long long>(k + 1)) + "." +
                             std::to_string(static_cast<long long>(ks + 1)) + ")",
                         level);
      else
        for (size_t k = 0; k < nr; k++) {
          size_t lmax = nc;
          if (m_parity != parityNone && ijsym == 0)
            lmax = m_parity == parityEven ? k + 1 : k;
          for (size_t l = 0; l < lmax; l++)
            s << " "
              << smat(ijsym, ks, k, l)
                     ->str(title + "(" + std::to_string(static_cast<long long>(k + 1)) + "." +
                               std::to_string(static_cast<long long>(ks + 1)) + "," +
                               std::to_string(static_cast<long long>(l + 1)) + "." +
                               std::to_string(static_cast<long long>((ks ^ ijsym) + 1)) + ")",
                           level)
              << std::endl;
        }
    }

  return s.str();
}

template <class T> void SMatMat_<T>::scal(typename T::value_type a, bool scaleDiagonal) {
  for (unsigned int ijsym = 0; ijsym < (rank() > 1 ? max_symmetry_ : 1); ijsym++)
    for (unsigned int ks = 0; ks < max_symmetry_; ks++) {
      size_t nr = m_dimensions[0][ks], nc = (rank() > 1 ? m_dimensions[1][ks ^ ijsym] : 1);
      if (nr < 1 || nc < 1)
        continue;
      if (rank() == 2 && m_parity != parityNone && ks < (ks ^ ijsym))
        continue;
      if (rank() == 1)
        for (size_t k = 0; k < nr; k++)
          smat(ijsym, ks, k)->scal(a, scaleDiagonal);
      else
        for (size_t k = 0; k < nr; k++) {
          size_t lmax = nc;
          if (m_parity != parityNone && ijsym == 0)
            lmax = m_parity == parityEven ? k + 1 : k;
          for (size_t l = 0; l < lmax; l++)
            smat(ijsym, ks, k, l)->scal(a, scaleDiagonal);
        }
    }
}

template <class T> molpro::array<typename T::value_type> SMatMat_<T>::block(unsigned int block_symmetry) const {
  size_t rows = m_templates[block_symmetry ^ m_symmetry]->size();
  return molpro::array<typename T::value_type>(&((*m_buffer)[block_offset(block_symmetry)]),
                                               block_size(block_symmetry) * rows);
}

#ifdef EIGEN_CORE_H
using namespace Eigen;
template <class T> typename SMatMat_<T>::M SMatMat_<T>::blockM(unsigned int block_symmetry) const {
  if (rank() == 1)
    throw std::logic_error("SMat::blockM call on rank=1 not yet implemented");
  size_t rows = m_templates[block_symmetry ^ m_symmetry]->size();
  Map<Matrix<typename T::value_type, Eigen::Dynamic, Eigen::Dynamic>, Unaligned, Stride<Dynamic, Dynamic>> result2(
      &((*m_buffer)[block_offset(block_symmetry)]), rows, block_size(block_symmetry),
      Stride<Dynamic, Dynamic>(rows, 1));
  return result2;
}

#endif

template <class T> molpro::bytestream SMatMat_<T>::bytestream(bool data) {
  std::vector<std::vector<char>> bs_templates;
  //  std::cout << "before template dump push"<<std::endl;
  for (unsigned int sym = 0; sym < 8; sym++) {
    auto bst = m_templates[sym]->bytestream(false);
    bs_templates.emplace_back(std::vector<char>(bst.data().begin(), bst.data().end()));
  }
  class molpro::bytestream bs;
  bs.append(m_symmetry);
  bs.append(static_cast<int>(m_parity));
  for (int i = 0; i < 2; i++)
    bs.append(&m_dimensions[i][0], m_dimensions[i].size());
  molpro::vector<char> dd(m_description.size());
  std::copy(m_description.begin(), m_description.end(), dd.begin());
  bs.append(&dd[0], dd.size());
  for (unsigned int sym = 0; sym < 8; sym++) {
    class molpro::bytestream::bytestream bss(&bs_templates[sym][0]);
    //    std::cout <<"dumping template bytestream "; bss.dump();
    bs.append(bss);
  }
  if (data)
    bs.append(&(*m_buffer)[0], m_buffer->size());
  return bs;
}

#endif // SYMMETRY_MATRIX_SMATMAT_IMPLEMENTATION_H
