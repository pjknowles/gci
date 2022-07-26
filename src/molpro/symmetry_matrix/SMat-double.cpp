#include "SMat-implementation.h"
namespace molpro {

template <>
int SMat_<double>::ev(SMat_<double>& val, SMat_<double>* vec, SMat_<double>* vali, SMat_<double>* vecl,
                      std::string algorithm, std::string sort) const {
  using T = double;
  FORTINT info = 0;
  if (this->rank() != 2 || m_dimensions[0] != m_dimensions[1])
    ErrorExit("Eigenvalues/vectors only for square matrix, SMat_::ev");
  if (this->m_symmetry != 0)
    ErrorExit("Eigenvalues/vectors only for matrix of symmetry 1, SMat_::ev");
  if (val.rank() != 2 || !val.m_diagonal || val.m_buffer == nullptr || val.m_symmetry != 0 ||
      m_dimensions[0] != val.m_dimensions[0])
    ErrorExit("Invalid val, SMat_::ev");
  val.m_description = "Eigenvalues";
  if (vec != nullptr &&
      (vec->rank() != 2 || vec->m_buffer == nullptr || vec->m_symmetry != 0 || m_dimensions != vec->m_dimensions))
    ErrorExit("Invalid vec, SMat_::ev");
  if (vali != nullptr && (vali->rank() != 2 || !vali->m_diagonal || vali->m_buffer == nullptr ||
                          vali->m_symmetry != 0 || m_dimensions != vali->m_dimensions))
    ErrorExit("Invalid vali, SMat_::ev");
  if (vecl != nullptr &&
      (vecl->rank() != 2 || vecl->m_buffer == nullptr || vecl->m_symmetry != 0 || m_dimensions != vecl->m_dimensions))
    ErrorExit("Invalid vecl, SMat_::ev");
  SMat_<T>* vvpt = vec;
  if (vec == nullptr)
    vvpt = new SMat_<T>(this, parityNone);
  *vvpt = *this;
  vvpt->m_description = "Eigenvectors";
  for (unsigned int k = 0; k < max_symmetry_; k++) {
    size_t n = this->m_dimensions[0][k];
    if (n < 1)
      continue;
    molpro::array<T> vblock = vvpt->block(k);
    molpro::array<T> valblock = val.block(k);
    if (this->m_parity > 0) {
#ifdef SMAT_USE_LAPACK
      char jobz = (vec == nullptr ? 'N' : 'V');
      T lwork;
      dsyev_X(jobz, 'L', n, &vblock[0], n, &valblock[0], &lwork, -1, info);
      {
        molpro::array<T> work((size_t)(lwork + 1));
        dsyev_X(jobz, 'L', n, &vblock[0], n, &valblock[0], &work[0], (size_t)lwork, info);
      }
      if (info != 0) {
        return info;
      }
#else // Eigen
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> es(
          vvpt->blockMap(k), vec != nullptr ? Eigen::ComputeEigenvectors : 0);
      val.blockV(k) = es.eigenvalues();
      if (vec != nullptr)
        vvpt->blockMap(k) = es.eigenvectors();
#endif
    } else if (this->m_parity == 0) {
      molpro::array<T> mptc(this->block(k)); // point
      molpro::array<T> mpt(mptc);            // copy
      molpro::array<T> vpt(n);
      molpro::array<T> mpt2(1);
      if (vecl != nullptr)
        mpt2 = vecl->block(k);
      molpro::array<T> mpt3(1);
      if (vec != nullptr)
        mpt3 = vec->block(k);
#ifdef SMAT_USE_LAPACK
      T lwork;
      char jobvl = (vecl == nullptr ? 'N' : 'V');
      dgeev_X(jobvl, jobz, n, &mpt[0], n, &valblock[0], &vpt[0], &mpt2[0], n, &mpt3[0], n, &lwork, -1, info);
      {
        molpro::array<T> work((size_t)(lwork + 1));
        if (Smat_debug)
          xout << "before second dgeev " << n << mpt << std::endl;
        dgeev_X(jobvl, jobz, n, &mpt[0], n, &valblock[0], &vpt[0], &mpt2[0], n, &mpt3[0], n, &work[0], (size_t)lwork,
                info);
      }
#else // Eigen
      EigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> es(blockMap(k), vec != nullptr);
      val.blockV(k) = es.eigenvalues().real();
      auto vmax = es.eigenvalues().real().lpNorm<Eigen::Infinity>();
      V(&vpt[0], n) = es.eigenvalues().imag();
      if (vec != nullptr)
        vec->blockMap(k) = es.eigenvectors().real(); // flaky
      if (vecl != nullptr)
        throw std::invalid_argument("SMat::ev left eigenvectors are not implemented");
#endif
      for (size_t m = 0; m < n; m++)
        if (vpt[m] != 0.0 && std::fabs(vmax / vpt[m]) < 1e14) {
          if (vali == nullptr)
            ErrorExit("Matrix has complex eigenvalues, but no array to receive imaginary part was passed, SMat_::ev");
          for (size_t l = 0; l < n; l++)
            vali->block(m)[l] = vpt[l];
          goto validone;
        }
    validone:;
    } else
      ErrorExit(" Eigenvalues/vectors cannot be computed for antisymmetric matrix, SMat_::ev");
    // sort eigensolutions
    std::vector<long> map(n, -1);
    {
      molpro::array<T> vpt(n);
      if (sort[0] == 'A' || sort[0] == 'a')
        for (size_t l = 0; l < n; l++)
          vpt[l] = -valblock[l];
      else if (sort[0] == 'D' || sort[0] == 'd')
        for (size_t l = 0; l < n; l++)
          vpt[l] = valblock[l];
      else if (sort[0] == 'O' || sort[0] == 'o') {
        if (vec == nullptr)
          ErrorExit("Sorting of eigensolutions by overlap requested, but eigenvectors not calculated, SMat_::ev");
        for (size_t l = 0; l < n; l++) {
          typename SMat_<T>::scalar_type tester = 0;
          for (size_t m = 0; l < n; l++)
            if (std::abs(vec->block(k)[m + l * n]) > tester) {
              tester = std::fabs(vec->block(k)[m + l * n]);
              vpt[l] = (T)m;
            }
        }
      } else
        ErrorExit("Unknown sorting algorithm, SMat_::ev");

      // naive sort
      {
        std::vector<bool> chosen(n, false);
        for (size_t l = 0; l < n; l++) {
          typename SMat_<T>::scalar_type tester = -(std::numeric_limits<double>::max()) / 4;
          for (size_t m = 0; m < n; m++) {
            if (chosen[m] || vpt[m] < tester)
              continue;
            tester = vpt[m];
            map[l] = m;
          }
          chosen[map[l]] = true;
        }
      }
    }

    // sort eigenvalues
    {
      molpro::array<T> vpt(n);
      for (size_t m = 0; m < n; m++)
        vpt[m] = valblock[map[m]];
      for (size_t m = 0; m < n; m++)
        valblock[m] = vpt[m];
    }
    if (vali != nullptr) {
      molpro::array<T> valiblock = vali->block(k);
      molpro::array<T> vpt(n);
      for (size_t m = 0; m < n; m++)
        vpt[m] = valiblock[map[m]];
      for (size_t m = 0; m < n; m++)
        valiblock[m] = vpt[m];
    }
    // sort eigenvectors
    if (vec != nullptr) {
      molpro::array<T> mpt = vec->block(k);
      molpro::array<T> mpt2(n * n);
      for (size_t l = 0; l < n; l++) {
        size_t iphase = 0;
        for (size_t m = 0; m < n; m++)
          if (std::fabs(mpt[m + map[l] * n]) > std::fabs(mpt[iphase + map[l] * n]))
            iphase = m;
        if (mpt[iphase + map[l] * n] > 0)
          for (size_t m = 0; m < n; m++)
            mpt2[m + l * n] = mpt[m + map[l] * n];
        else
          for (size_t m = 0; m < n; m++)
            mpt2[m + l * n] = -mpt[m + map[l] * n];
      }
      for (size_t m = 0; m < n * n; m++)
        mpt[m] = mpt2[m];
    }
    if (vecl != nullptr) {
      molpro::array<T> mpt = vecl->block(k);
      molpro::array<T> mpt2(n * n);
      for (size_t l = 0; l < n; l++) {
        size_t iphase = 0;
        for (size_t m = 0; m < n; m++)
          if (std::fabs(mpt[m + map[l] * n]) > std::fabs(mpt[iphase + map[l] * n]))
            iphase = m;
        if (mpt[iphase + map[l] * n] > 0)
          for (size_t m = 0; m < n; m++)
            mpt2[m + l * n] = mpt[m + map[l] * n];
        else
          for (size_t m = 0; m < n; m++)
            mpt2[m + l * n] = -mpt[m + map[l] * n];
      }
      for (size_t m = 0; m < n * n; m++)
        mpt[m] = mpt2[m];
    }
  }
  if (vec == nullptr)
    delete vvpt;
  return info;
}

template <class T> SMat_<T> SMat_<T>::solve(const SMat_<T>& rhs, std::string algorithm) const {
  auto mat = const_cast<SMat_<T>*>(this);
  // non-zero parity on input causes a copy
  if (parity() != 0) {
    mat = new SMat_<T>(this, parityNone);
    *mat = *this;
  }
  if (mat->rank() != 2)
    ErrorExit("SMat::solve can be called only on a matrix, not a vector");
  int itranspose = mat->m_transposed ? 1 : 0;
  int itransposer = rhs.m_transposed ? 1 : 0;
  if (mat->m_dimensions[itranspose] != rhs.m_dimensions[itransposer])
    ErrorExit("Dimension mismatch, SMat_::solve");
  SMat_<T> result(dims_t{mat->m_dimensions[1 - itranspose], rhs.m_dimensions[1 - itransposer]}, parityNone,
                  m_symmetry ^ rhs.m_symmetry, false, std::string());
  for (unsigned int k = 0; k < 8; k++) {
    unsigned int l = k ^ mat->m_symmetry;
    auto A = mat->blockMap(k);
    if (A.size() == 0)
      continue;
    auto b = rhs.blockMap(k);
    if (b.size() == 0)
      continue;
    auto x = result.blockMap(l);
    if (x.size() == 0)
      continue;
    if (algorithm == "partialPivLU")
      x = A.partialPivLu().solve(b);
    else if (algorithm == "FullPivLU")
      x = A.fullPivLu().solve(b);
    else if (algorithm == "HouseholderQR")
      x = A.householderQr().solve(b);
    else if (algorithm == "ColPivHouseholderQR")
      x = A.colPivHouseholderQr().solve(b);
    else if (algorithm == "FullPivHouseholderQR")
      x = A.fullPivHouseholderQr().solve(b);
    else if (algorithm == "LLT")
      x = A.llt().solve(b);
    else if (algorithm == "LDLT")
      x = A.ldlt().solve(b);
    else {
      std::cout << algorithm << std::endl;
      ErrorExit("Unknown algorithm in SMat_::solve");
    }
  }
  if (parity() != 0)
    delete mat;
  return result;
}

template <> molpro::bytestream SMat_<double>::bytestream(bool data) {
  class molpro::bytestream bs;
  //  xout << "bs constructed, size="<<bs.size()<<std::endl;
  bs.append(m_symmetry);
  //  xout << "appended symmetry_, size="<<bs.size()<<std::endl;
  //  memory_print_status();
  bs.append(m_parity);
  bs.append((m_transposed ? 1 : 0));
  bs.append(rank());
  bs.append(&m_dimensions[0][0], m_dimensions[0].size());
  if (rank() > 1)
    bs.append(&m_dimensions[1][0], m_dimensions[1].size());
  molpro::vector<char> dd(m_description.size());
  std::copy(m_description.begin(), m_description.end(), dd.begin());
  bs.append(&dd[0], dd.size());
  //  xout << "SMat::dump before data append size()="<<bs.size()<<std::endl;
  if (data)
    bs.append(&(*m_buffer)[0], m_buffer->size());
  //  xout << "SMat::dump after data append size()="<<bs.size()<<std::endl;
  return bs;
}

template <> SMat_<double>::SMat_(const char* dump, double* buffer) : m_diagonal(false) {
  using T = double;
  //  xout << "construct smat from bytestream T*"<<std::endl;
  class molpro::bytestream bs(dump);
  //  printf("input bytestream hash=%lX\n",bs.hash());
  m_symmetry = bs.ints()[0];
  //  xout << "recovered symmetry="<<m_symmetry<<std::endl;
  m_parity = static_cast<parity_t>(bs.ints()[0]);
  //  xout << "m_parity="<<m_parity<<std::endl;
  m_transposed = bs.ints()[0] != 0;
  //  xout << "m_transposed="<<m_transposed<<std::endl;
  int rank = bs.ints()[0];
  //  xout << "rank="<<rank<<std::endl;
  for (int axis = 0; axis < rank; axis++) {
    //      xout << "axis="<<axis<<std::endl;
    auto ii = bs.ints();
    m_dimensions.push_back(dim_t(ii.size()));
    for (size_t k = 0; k < ii.size(); k++)
      m_dimensions.back()[k] = ii[k];
    //      for (size_t k=0; k<ii.size();k++) xout << "dimension "<<k<<" "<<m_dimensions.back()[k]<<std::endl;
  }
  auto dd = bs.chars();
  //  xout <<"back from bs.chars()"<<dd<<std::endl;
  m_description.resize(dd.size());
  std::copy(dd.begin(), dd.end(), m_description.begin());
  //  xout << "m_description="<<m_description<<std::endl;
  m_managed_buffer = true;
  size_t n; // For dumb pgc++ - bug 5086
  if (buffer == nullptr)
    this->m_bufferp = std::make_shared<molpro::array<T>>(n = size()); // no checks on size of buffer!
  else
    this->m_bufferp = std::make_shared<molpro::array<T>>(buffer, n = size()); // no checks on size of buffer!
  this->m_buffer = &this->m_bufferp.get()[0];
  if (bs.size() > bs.position()) {
    //      xout << "bs.size()="<<bs.size()<<", bs.position()="<<bs.position()<<std::endl;
    auto data = bs.doubles(); // FIXME for non-double
                              //      xout << "importing data of size "<<data.size()<<std::endl;
    std::copy(data.begin(), data.end(), this->m_buffer->begin());
    //      xout << "bs.size()="<<bs.size()<<", bs.position()="<<bs.position()<<std::endl;
  }
}
} // namespace molpro

template class molpro::SMat_<double>;
