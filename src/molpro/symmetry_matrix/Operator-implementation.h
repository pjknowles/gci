#ifndef SYMMETRY_MATRIX_OPERATOR_IMPLEMENTATION_H
#define SYMMETRY_MATRIX_OPERATOR_IMPLEMENTATION_H
#include <algorithm>
#include <molpro/symmetry_matrix/Operator.h>

using namespace molpro;

inline parity_t parity(int hermiticity, parity_t odd = parityOddPacked) {
  if (hermiticity == 0)
    return parityNone;
  if (hermiticity > 0)
    return parityEven;
  return odd;
}

template <class T>
Operator_<T>::Operator_(std::array<dims_t, 2> dimensions, int rank, bool uhf, std::vector<int> hermiticity,
                        std::vector<int> exchange, unsigned int symmetry, bool covariant, bool diagonal,
                        std::string description)
    : m_dimensions(std::move(dimensions)), m_rank(std::move(rank)), m_uhf(std::move(uhf)),
      m_hermiticity(std::move(hermiticity)), m_exchange(std::move(exchange)), m_symmetry(std::move(symmetry)),
      m_description(std::move(description)), m_O0(0), m_covariant(covariant), m_diagonal(diagonal),
      m_dirac_out_of_date(true) {
  if (m_dimensions.size() != 2)
    throw std::logic_error("dimensions must be dimensioned 2");
  if (int(m_dimensions[0].size()) < m_rank * 2)
    throw std::logic_error("dimensions must be dimensioned at least rank*2");
  if (int(m_dimensions[1].size()) < m_rank * 2)
    throw std::logic_error("dimensions must be dimensioned at least rank*2");
  if (m_rank > 1 && m_hermiticity.size() != 2)
    throw std::logic_error("hermiticity must be dimensioned 2");
  if (int(m_exchange.size()) < (m_rank - 1) * 2)
    throw std::logic_error("exchange must be dimensioned at least 2");
  T* dummyBuffer = nullptr;
  if (m_rank == 1 || (m_rank == 2 && m_hermiticity[0] == m_hermiticity[1])) {
    m_O1a = std::make_shared<SMat_<T>>(dims_t{m_dimensions[0][0], m_dimensions[0][1]},
                                       parity(m_hermiticity[0] * (rank < 2 ? 1 : m_hermiticity[1]), parityOdd),
                                       m_symmetry, m_diagonal, m_description + (m_uhf ? ":a" : ""));
    m_O1b = m_uhf ? std::make_shared<SMat_<T>>(dims_t{m_dimensions[1][0], m_dimensions[1][1]},
                                               parity(m_hermiticity[0] * (rank < 2 ? 1 : m_hermiticity[1]), parityOdd),
                                               m_symmetry, false, m_description + ":b")
                  : m_O1a;
  }
  if (m_rank == 2) {
    {
      SMat_<T> innera(dims_t{m_dimensions[0][0], m_dimensions[0][1]}, dummyBuffer, parity(m_hermiticity[0]), 0, false,
                      m_description + (m_uhf ? ":maa" : "m"));
      SMat_<T> innerb(dims_t{m_dimensions[1][0], m_dimensions[1][1]}, dummyBuffer, parity(m_hermiticity[0]), 0, false,
                      m_description + (m_uhf ? ":mbb" : "m"));
      //      std::cout << "memory_remaining()="<<memory_remaining()<<std::endl;
      m_O2aa_mulliken = std::make_shared<SMatMat_<SMat_<T>>>(innera, dims_t{m_dimensions[0][2], m_dimensions[0][3]},
                                                             nullptr, parity(m_hermiticity[1]), m_symmetry,
                                                             m_description + (m_uhf ? ":maa" : ":m"));
      //      std::cout << "O2aa_mulliken created at "<<&m_O2aa_mulliken->block(0)[0]<<std::endl;
      //      std::cout << "memory_remaining()="<<memory_remaining()<<std::endl;
      if (m_uhf) {
        m_O2ab_mulliken =
            std::make_shared<SMatMat_<SMat_<T>>>(innera, dims_t{m_dimensions[1][2], m_dimensions[1][3]}, nullptr,
                                                 parity(m_hermiticity[1]), m_symmetry, m_description + ":mab");
        if (m_exchange[0] * m_exchange[1] == 0)
          m_O2ba_mulliken =
              std::make_shared<SMatMat_<SMat_<T>>>(innerb, dims_t{m_dimensions[0][2], m_dimensions[0][3]}, nullptr,
                                                   parity(m_hermiticity[1]), m_symmetry, m_description + ":mba");
        m_O2bb_mulliken =
            std::make_shared<SMatMat_<SMat_<T>>>(innerb, dims_t{m_dimensions[1][2], m_dimensions[1][3]}, nullptr,
                                                 parity(m_hermiticity[1]), m_symmetry, m_description + ":mbb");
        //            std::cout << "UHF Operator_<T> load "<<std::endl;
      } else {
        m_O2ab_mulliken = m_O2ba_mulliken = m_O2bb_mulliken = m_O2aa_mulliken;
        //            std::cout << "RHF Operator_<T> load "<<std::endl;
        //            std::cout << m_O2aa_mulliken << std::endl;
        //            std::cout << m_O2ab_mulliken << std::endl;
      }
    }
    {
      SMat_<T> innera(dims_t{m_dimensions[0][0], m_dimensions[0][2]}, dummyBuffer, parity(m_exchange[0]), 0, false,
                      m_description + ":daa");
      SMat_<T> innerb(dims_t{m_dimensions[1][0], m_dimensions[1][2]}, dummyBuffer, parity(m_exchange[0]), 0, false,
                      m_description + ":dbb");
      //        std::cout << "inner parity="<<inner.parity()<<std::endl;
      //        std::cout << "O2aa_mulliken at "<<&m_O2aa_mulliken->block(0)[0]<<std::endl;
      //        std::cout << "memory_remaining()="<<memory_remaining()<<std::endl;
      m_O2aa_dirac = std::make_shared<SMatMat_<SMat_<T>>>(
          innera, dims_t{m_dimensions[0][1], m_dimensions[0][3]}, nullptr,
          (m_exchange[1] ? (m_exchange[1] > 0 ? parityEven : parityOddPacked) : parityNone), m_symmetry,
          m_description + ":daa");
      //        std::cout << "O2aa_dirac created at "<<&m_O2aa_dirac->block(0)[0]<<std::endl;
      //        std::cout << "memory_remaining()="<<memory_remaining()<<std::endl;
      SMat_<T> inner2(dims_t{m_dimensions[0][0], m_dimensions[1][2]}, dummyBuffer, parityNone, 0, false,
                      m_description + ":dab");
      //        std::cout << "inner2 parity="<<inner2.parity()<<std::endl;
      m_O2ab_dirac = std::make_shared<SMatMat_<SMat_<T>>>(inner2, dims_t{m_dimensions[0][1], m_dimensions[1][3]},
                                                          nullptr, parityNone, m_symmetry, m_description + ":dab");
      //        std::cout << "O2ab_dirac created at "<<&m_O2ab_dirac->block(0)[0]<<std::endl;
      //        std::cout << "memory_remaining()="<<memory_remaining()<<std::endl;
      if (m_uhf) {
        m_O2bb_dirac = std::make_shared<SMatMat_<SMat_<T>>>(
            innerb, dims_t{m_dimensions[1][1], m_dimensions[1][3]}, nullptr,
            (m_exchange[1] ? (m_exchange[1] > 0 ? parityEven : parityOddPacked) : parityNone), m_symmetry,
            m_description + ":dbb");
        SMat_<T> inner2(dims_t{m_dimensions[1][0], m_dimensions[0][2]}, dummyBuffer, parityNone, 0, false,
                        m_description + ":dba");
        if (m_exchange[0] * m_exchange[1] == 0)
          m_O2ba_dirac = std::make_shared<SMatMat_<SMat_<T>>>(inner2, dims_t{m_dimensions[1][1], m_dimensions[0][3]},
                                                              nullptr, parityNone, m_symmetry, m_description + ":dba");
      } else {
        m_O2bb_dirac = m_O2aa_dirac;
        if (m_exchange[0] * m_exchange[1] == 0)
          m_O2ba_dirac = m_O2ab_dirac;
      }
      m_dirac_out_of_date = true;
    }
  }
  //  std::cout << "at end of Operator_<T>::Operator_<T>, memory_remaining()="<<memory_remaining()<<std::endl;
}

template <class T>
Operator_<T>::Operator_(const Operator_<T>& source)
    : Operator_<T>::Operator_(source.m_dimensions, source.m_rank, source.m_uhf, source.m_hermiticity, source.m_exchange,
                              source.m_symmetry, source.m_covariant, source.m_diagonal, source.m_description) {
  *this = source;
}

template <class T> Operator_<T>& Operator_<T>::operator=(const Operator_<T>& source) {
  checkCompatible(source);
  m_O0 = source.m_O0;
  *m_O1a = *source.m_O1a;
  if (m_rank > 1)
    *m_O2aa_mulliken = *source.m_O2aa_mulliken;
  if (m_uhf) {
    *m_O1b = *source.m_O1b;
    if (m_rank > 1)
      *m_O2ab_mulliken = *source.m_O2ab_mulliken;
    if (m_rank > 1)
      *m_O2bb_mulliken = *source.m_O2bb_mulliken;
  }
  m_dirac_out_of_date = true;
  return *this;
}

template <class T> Operator_<T> Operator_<T>::construct(const char* dump) {
  class molpro::bytestream bs(dump);
  return Operator_<T>::construct(bs);
}

template <class T> Operator_<T> Operator_<T>::construct(class molpro::bytestream& bs) {
  std::array<dims_t, 2> dimensions;
  for (int spin = 0; spin < 2; spin++)
    for (int axis = 0; axis < 4; axis++) {
      auto ii = bs.ints();
      dimensions[spin].push_back(dim_t(ii.size()));
      for (size_t k = 0; k < ii.size(); k++)
        dimensions[spin].back()[k] = ii[k];
    }
  int rank;                     //! the rank of the operator
  bool uhf;                     //!< Whether the underlying 1-particle spaces are different for alpha and beta spin.
  std::vector<int> hermiticity; //!< Conjugation symmetry, [ai] -> [ia], [bj] -> [jb]: 0=none, 1=symmetric,
                                //!< -1=antisymmetric, default 0.
  std::vector<int> exchange;
  //!< For each of bra, ket, whether there is an interchange symmetry between two
  //! indices in the two-particle part of the operator: [abij] -> [baij]: 0=none, 1=symmetric, -1=antisymmetric, default
  //! -1.
  unsigned int symmetry;   //!< The spatial symmetry of the operator.
  std::string description; //!< A string describing the object.
  rank = bs.ints()[0];
  uhf = bs.ints()[0];
  {
    auto ii = bs.ints();
    hermiticity.push_back(ii[0]);
    hermiticity.push_back(ii[1]);
  }
  {
    auto ii = bs.ints();
    exchange.push_back(ii[0]);
    exchange.push_back(ii[1]);
  }
  symmetry = bs.ints()[0];
  bool covariant = bs.ints()[0] != 0;
  {
    auto dd = bs.chars();
    description.resize(dd.size());
    std::copy(dd.begin(), dd.end(), description.begin());
  }
  bool dirac_out_of_date = bs.ints()[0] != 0;
  Operator_<T> result(dimensions, rank, uhf, hermiticity, exchange, symmetry, covariant, false, description);
  result.m_O0 = bs.doubles()[0];
  if (result.m_rank > 0) {
    result.m_O1a = std::make_shared<SMat_<T>>(bs.byteStream());
    result.m_O1b = result.m_uhf ? std::make_shared<SMat_<T>>(bs.byteStream()) : result.m_O1a;
  }
  if (result.m_rank > 1) {
    result.m_O2aa_mulliken = std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream());
    result.m_O2ab_mulliken =
        result.m_uhf ? std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream()) : result.m_O2aa_mulliken;
    result.m_O2bb_mulliken =
        result.m_uhf ? std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream()) : result.m_O2aa_mulliken;
    if (!dirac_out_of_date) {
      result.m_O2aa_dirac = std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream());
      result.m_O2ab_dirac = result.m_uhf ? std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream()) : result.m_O2aa_dirac;
      result.m_O2bb_dirac = result.m_uhf ? std::make_shared<SMatMat_<SMat_<T>>>(bs.byteStream()) : result.m_O2aa_dirac;
      result.m_dirac_out_of_date = false;
    }
  }
  return result;
}

template <class T>
class Operator_<T> Operator_<T>::slice(std::array<dims_t, 2> dimensions, std::array<dims_t, 2> offset,
                                       std::string description) const {
  auto hermiticity = m_hermiticity;
  for (auto spin = 0; spin < 2; spin++)
    for (auto particle = 0; particle < m_rank; particle++)
      if (dimensions[spin][particle * 2] != dimensions[spin][particle * 2 + 1])
        hermiticity[particle] = 0;
  auto exchange = m_exchange;
  if (m_rank > 1) {
    for (auto spin = 0; spin < 2; spin++)
      for (auto braket = 0; braket < 2; braket++)
        if (dimensions[spin][braket] != dimensions[spin][braket + 2])
          exchange[braket] = 0;
  }
  Operator_ result(dimensions, this->m_rank, this->m_uhf, hermiticity, exchange, this->m_symmetry, this->m_covariant,
                   this->m_diagonal, description == "" ? this->m_description : description);
  *result.m_O1a = m_O1a->slice(dimensions[0], offset[0]);
  if (m_uhf)
    *result.m_O1b = m_O1b->slice(dimensions[1], offset[1]);
  if (m_rank > 1) {
    throw std::logic_error("no implementation");
  }
  return result;
}

template <class T>
void Operator_<T>::ensure_dirac() const {
  if (m_rank > 1)
    if (m_dirac_out_of_date && m_covariant) {
      m_dirac_out_of_date = false;
      if (false)
        std::cout << O2(true, true, true) << std::endl;
      if (false) {
        for (unsigned int symik = 0; symik < 8; symik += 4)
          for (unsigned int symi = 0; symi < 8; symi += 4)
            for (unsigned int i = 0; i < m_dimensions[0][0][symi]; i++)
              for (unsigned int k = 0; k < m_dimensions[0][0][symik ^ symi]; k++) {
                std::cout << "@@ Mulliken integrals for i=" << i << "." << symi << ", k=" << k << "." << (symi ^ symik)
                          << std::endl;
                std::cout << *(O2(true, true, true).smat(symik, symi, i, k)) << std::endl;
              }
      }
      if (false) {
        unsigned int i = 0, symi = 0, k = 1, symk = 4, symik = 4;
        std::cout << "@@ Mulliken integrals for i=" << i << "." << symi << ", k=" << k << "." << (symi ^ symik)
                  << std::endl;
        std::cout << (O2(true, true, true).smat(symik, symi, i, k))->str("", 3) << std::endl;
        //        std::swap(i,k);
        std::swap(symi, symk);
        std::cout << "@@ Mulliken integrals for i=" << i << "." << symi << ", k=" << k << "." << (symi ^ symik)
                  << std::endl;
        std::cout << (O2(true, true, true).smat(symik, symi, i, k))->str("", 3) << std::endl;
      }
      if (false) {
        T el1 = element(1, 0, 0, 4, 0, 0, 1, 4);
        T el2 = element(1, 4, 0, 0, 0, 4, 1, 0);
        std::cout << "@@@ el1=" << el1 << std::endl;
        std::cout << "@@@ el2=" << el2 << std::endl;
        //      exit(0);
      }
      if (m_hermiticity[0] * m_hermiticity[1] != 1 || m_exchange[0] * m_exchange[1] != 1 || m_symmetry != 0)
        throw std::logic_error(
            "Dirac ordering of integrals not yet implemented for this combination of hermiticity/exchange symmetry");
      for (unsigned int symik = 0; symik < 8; symik++) {
        if (m_O2ab_dirac->block_size(symik) <= 0)
          continue;
        //          std::cout << "@@symik="<<symik<<std::endl;
        for (unsigned int symi = 0; symi < 8; symi++) {
          for (unsigned int symj = 0; symj < 8; symj++) {
            unsigned int symk = symi ^ symik;
            unsigned int syml = symj ^ symik;
            //                  std::cout << "symi="<<symi<<" symj="<<symj<< " symk="<<symk<<" syml="<<syml<<std::endl;
            //                  if (symi <= symk && symj <= syml)
            {
              for (size_t j = 0; j < m_dimensions[0][2][symj]; j++)
                for (size_t l = 0; l < m_dimensions[1][3][syml]; l++) {
                  size_t ki = 0;
                  for (size_t i = 0; i < m_dimensions[0][0][symi]; i++)
                    for (size_t k = 0; k < m_dimensions[1][1][symk]; k++) {
                      //                            std::cout << "opposite spin
                      //                            <"<<i<<"."<<symi<<","<<k<<"||"<<j<<l<<">="<<
                      //                            "("<<i<<j<<"|"<<k<<l<<") = " <<
                      //                            element(i,symi,j,symj,k,symk,l,syml,true,false)<<std::endl;
                      O2(true, false, false, false).smat(symik, syml, l, j)->block(symk)[ki] =
                          element(i, symi, j, symj, k, symk, l, syml, true, false);
                      ki++;
                    }
                }
              for (size_t j = 0; j < m_dimensions[0][2][symj]; j++)
                for (size_t l = 0; l < (symik ? m_dimensions[0][3][syml] : j); l++) {
                  size_t ki = 0;
                  for (size_t k = 0; k < m_dimensions[0][1][symk]; k++)
                    for (size_t i = 0; i < (symik ? m_dimensions[0][0][symi] : k); i++) {
                      O2(true, true, false, false).smat(symik, syml, l, j)->block(symk)[ki] =
                          symik ? (element(i, symi, j, symj, k, symk, l, syml, true, true) -
                                   element(i, symi, l, syml, k, symk, j, symj, true, true))
                                : (element(i, symi, l, syml, k, symk, j, symj, true, true) -
                                   element(i, symi, j, symj, k, symk, l, syml, true, true));
                      ki++;
                    }
                }
              for (size_t j = 0; j < m_dimensions[1][2][symj]; j++)
                for (size_t l = 0; l < (symik ? m_dimensions[1][3][syml] : j); l++) {
                  size_t ki = 0;
                  for (size_t k = 0; k < m_dimensions[1][1][symk]; k++)
                    for (size_t i = 0; i < (symik ? m_dimensions[1][0][symi] : k); i++) {
                      O2(true, true, false, false).smat(symik, syml, l, j)->block(symk)[ki] =
                          symik ? (element(i, symi, j, symj, k, symk, l, syml, true, true) -
                                   element(i, symi, l, syml, k, symk, j, symj, true, true))
                                : (element(i, symi, l, syml, k, symk, j, symj, true, true) -
                                   element(i, symi, j, symj, k, symk, l, syml, true, true));
                      O2(false, false, false, false).smat(symik, syml, l, j)->block(symk)[ki] =
                          symik ? (element(i, symi, j, symj, k, symk, l, syml, false, false) -
                                   element(i, symi, l, syml, k, symk, j, symj, false, false))
                                : (element(i, symi, l, syml, k, symk, j, symj, false, false) -
                                   element(i, symi, j, symj, k, symk, l, syml, false, false));
                      ki++;
                    }
                }
            }
          }
        }
        //          std::cout << " block_offsets:"; for (auto ii=0; ii<8; ii++) std::cout <<"
        //          "<<m_O2aa_dirac->block_offset(symik,ii); std::cout << std::endl; if
        //          (m_O2aa_dirac->blockM(symik).rows()>0)
        //            std::cout << "Dirac aa symmetry block
        //            "<<&m_O2aa_dirac->blockM(symik)(0,0)-&m_O2aa_dirac->blockM(0)(0,0)<<std::endl<<m_O2aa_dirac->blockM(symik)<<std::endl;
        //          std::cout << "Dirac ab symmetry block "<<std::endl<<m_O2ab_dirac->blockM(symik)<<std::endl;
        //                    std::cout <<
        //                    "m_O2aa_dirac-element="<<&m_O2aa_dirac->block(0)[0]-&m_O2aa_mulliken->block(0)[0]<<std::endl;
      }
    }
}
template <class T> void Operator_<T>::mulliken_from_dirac() {
  if (m_rank > 1) {
    //      std::cout << "@@ mulliken_from_dirac"<<std::endl;
    //          std::cout << "Dirac ab symmetry block "<<std::endl<<m_O2ab_dirac->blockM(0)<<std::endl;
    //          std::cout << "Dirac ab "<<std::endl<<O2(true,false,false)<<std::endl;
    //      std::cout <<str("at start of from_dirac",2)<<std::endl;
    if (false)
      std::cout << O2(true, true, true) << std::endl;
    if (m_hermiticity[0] != m_hermiticity[1] || m_exchange[0] * m_exchange[1] != 1 || m_symmetry != 0)
      throw std::logic_error(
          "Dirac ordering of integrals not yet implemented for this combination of hermiticity/exchange symmetry");
    for (unsigned int symik = 0; symik < 8; symik++) {
      if (m_O2ab_dirac->block_size(symik) <= 0)
        continue;
      for (unsigned int symi = 0; symi < 8; symi++) {
        for (unsigned int symj = 0; symj < 8; symj++) {
          unsigned int symk = symi ^ symik;
          unsigned int syml = symj ^ symik;
          //          size_t ni = m_dimensions[0][symi];
          //          size_t nj = m_dimensions[0][symj];
          //          size_t nk = m_dimensions[0][symk];
          //          size_t nl = m_dimensions[0][syml];
          //          if (ni * ni * nk * nl == 0) continue;
          {
            if (m_uhf) {
              for (size_t j = 0; j < m_dimensions[0][2][symj]; j++)
                for (size_t l = 0; l < m_dimensions[1][3][syml]; l++) {
                  size_t ki = 0;
                  for (size_t i = 0; i < m_dimensions[0][0][symi]; i++)
                    for (size_t k = 0; k < m_dimensions[1][1][symk]; k++) {
                      element(i, symi, j, symj, k, symk, l, syml, true, false) =
                          O2(true, false, false)
                              .smat(symik, syml, l,
                                    j)
                              ->block(symk)[ki]; // TODO do we really have lj or jl?
                      //                                std::cout << "opposite spin
                      //                                <"<<i<<"."<<symi<<","<<k<<"||"<<j<<l<<">="<<
                      //                                "("<<i<<j<<"|"<<k<<l<<") = " <<
                      //                                element(i,symi,j,symj,k,symk,l,syml,true,false)<<std::endl;
                      ki++;
                    }
                }
              for (size_t j = 0; j < m_dimensions[0][2][symj]; j++)
                for (size_t l = 0; l < (symik ? m_dimensions[0][3][syml] : j); l++) {
                  //                      for (size_t l=(symik?0:j+1); l<nl; l++) {
                  size_t ki = 0;
                  for (size_t k = 0; k < m_dimensions[0][1][symk]; k++)
                    //                            for (size_t i=(symik ? 0 : k+1); i<ni; i++) {
                    for (size_t i = 0; i < (symik ? m_dimensions[0][0][symi] : k); i++) {
                      if (symik) {
                        element(i, symi, l, syml, k, symk, j, symj, true, true) =
                            -(element(i, symi, j, symj, k, symk, l, syml, true, true) =
                                  O2(true, true, false).smat(symik, syml, l, j)->block(symk)[ki]);
                      } else {
                        element(i, symi, j, symj, k, symk, l, syml, true, true) =
                            -(element(i, symi, l, syml, k, symk, j, symj, true, true) =
                                  O2(true, true, false).smat(symik, syml, l, j)->block(symk)[ki]);
                      }
                      ki++;
                    }
                }
              for (size_t j = 0; j < m_dimensions[1][2][symj]; j++)
                for (size_t l = 0; l < (symik ? m_dimensions[1][3][syml] : j); l++) {
                  //                      for (size_t l=(symik?0:j+1); l<nl; l++) {
                  size_t ki = 0;
                  for (size_t k = 0; k < m_dimensions[1][1][symk]; k++)
                    //                            for (size_t i=(symik ? 0 : k+1); i<ni; i++) {
                    for (size_t i = 0; i < (symik ? m_dimensions[1][0][symi] : k); i++) {
                      if (symik) {
                        element(i, symi, l, syml, k, symk, j, symj, false, false) =
                            -(element(i, symi, j, symj, k, symk, l, syml, false, false) =
                                  O2(false, false, false).smat(symik, syml, l, j)->block(symk)[ki]);
                      } else {
                        element(i, symi, j, symj, k, symk, l, syml, false, false) =
                            -(element(i, symi, l, syml, k, symk, j, symj, false, false) =
                                  O2(false, false, false).smat(symik, syml, l, j)->block(symk)[ki]);
                      }
                      ki++;
                    }
                }
            } else { // ! m_uhf
              if (m_symmetry != 0)
                throw std::runtime_error("Non-symmetric operator not yet implemented");
              //                        std::cout <<symi<<symj<<symk<<syml<<std::endl;
              //                        std::cout << "element(1,0,1,0,0,0,0,0,true,true,false)  " <<
              //                        element(1,0,1,0,0,0,0,0,true,true,false) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,1,0,0,0,true,true,false)  " <<
              //                        element(1,0,0,0,1,0,0,0,true,true,false) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,0,0,1,0,true,true,false)  " <<
              //                        element(1,0,0,0,0,0,1,0,true,true,false) <<std::endl;
              for (size_t i = 0; i < m_dimensions[0][0][symi]; i++)
                for (size_t j = 0; j < ((m_hermiticity[0] != 0 && symi == symj) ? i + (m_hermiticity[0] == -1 ? 0 : 1)
                                                                                : m_dimensions[0][2][symj]);
                     j++)
                  for (size_t k = 0; k < m_dimensions[0][1][symk]; k++) {
                    for (size_t l = 0;
                         l < ((m_hermiticity[1] != 0 && symk == syml) ? k + (m_hermiticity[0] == -1 ? 0 : 1)
                                                                      : m_dimensions[0][3][syml]);
                         l++) {
                      m_dirac_out_of_date = false;
                      //                                  size_t ki=k+i*nk;
                      //                                    std::cout << "spin summed <"
                      //                                             << "("<<i<<j<<"|"<<k<<l<<") = ";
                      element(i, symi, j, symj, k, symk, l, syml, true, true) =
                          element(i, symi, k, symk, l, syml, j, symj, true, false, false) +
                          element(k, symk, i, symi, j, symj, l, syml, true, false, false) +
                          element(i, symi, k, symk, l, syml, j, symj, true, true, false) +
                          element(i, symi, l, syml, k, symk, j, symj, true, true, false);
                    }
                  }
              //                        std::cout << "element(1,0,1,0,0,0,0,0,true,false,true)  " <<
              //                        element(1,0,1,0,0,0,0,0,true,false,true) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,1,0,0,0,true,false,true)  " <<
              //                        element(1,0,0,0,1,0,0,0,true,false,true) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,0,0,1,0,true,false,true)  " <<
              //                        element(1,0,0,0,0,0,1,0,true,false,true) <<std::endl; std::cout <<
              //                        "element(1,0,1,0,0,0,0,0,true,true,true)  " <<
              //                        element(1,0,1,0,0,0,0,0,true,true,true) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,1,0,0,0,true,true,true)  " <<
              //                        element(1,0,0,0,1,0,0,0,true,true,true) <<std::endl; std::cout <<
              //                        "element(1,0,0,0,0,0,1,0,true,true,true)  " <<
              //                        element(1,0,0,0,0,0,1,0,true,true,true) <<std::endl;
            }
          }
        }
      }
    }
    m_dirac_out_of_date = false;
    //      std::cout <<str("at end of from_dirac",2)<<std::endl;
  }
}

template <class T> const SMat_<T>& Operator_<T>::O1(bool spinUp) const {
  if (spinUp)
    return *m_O1a;
  else
    return *m_O1b;
}

template <class T> const SMatMat_<SMat_<T>>& Operator_<T>::O2(bool spinUp1, bool spinUp2, bool mulliken, bool ensure) const {
  if (mulliken)
    if (spinUp1)
      if (spinUp2)
        return *m_O2aa_mulliken;
      else
        return *m_O2ab_mulliken;
    else if (spinUp2)
      return *m_O2ba_mulliken;
    else
      return *m_O2bb_mulliken;
  else {
    if (ensure)
      ensure_dirac();
    if (spinUp1)
      if (spinUp2)
        return *m_O2aa_dirac;
      else
        return *m_O2ab_dirac;
    else if (spinUp2)
      return *m_O2ba_dirac;
    else
      return *m_O2bb_dirac;
  }
}

template <class T>
const T& Operator_<T>::element(int i, int isym, int j, int jsym, int k, int ksym, int l, int lsym, bool spinUp1,
                               bool spinUp2, bool mulliken) const {
  //  std::cout <<"=====\nelement
  //  ("<<i<<"."<<isym<<","<<j<<"."<<jsym<<","<<k<<"."<<ksym<<","<<l<<"."<<lsym<<")"<<spinUp1<<spinUp2<<mulliken<<std::endl;
  if (k < 0) { // 1-particle
    if (m_hermiticity[0])
      return O1(spinUp1).block(isym)[i > j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i];
    else
      return const_cast<Operator_<T>*>(this)->O1(spinUp1).blockMap(isym)(i, j);
  } else { // 2-particle
    if (((mulliken && m_hermiticity[0]) || (!mulliken && m_exchange[0] && spinUp1 == spinUp2)) && (ksym == lsym)) {
      auto m = O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l);
      //          std::cout << "smat\n"<<*m<<std::endl;
      if ((m->parity() == parityOdd || m->parity() == parityOddPacked) && (k == l || i == j))
        return constants.s_zero; // FIXME not thread-safe
      T phase = (i > j ? 1 : m_hermiticity[0]) * (k > l ? 1 : m_hermiticity[1]);
      if (m->parity() == parityOddPacked)
        return constants.return_value = phase * m->block(isym)[i > j ? i * (i - 1) / 2 + j : j * (j - 1) / 2 + i];
      else
        return constants.return_value = phase * m->block(isym)[i > j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i];
    } else if (isym > jsym) {
      if (false) {
        std::cout << "element rectangular block 1 (" << i << "," << j << ")" << std::endl
                  << *O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l) << std::endl
                  << O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(isym).rows() << std::endl
                  << O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(isym).cols() << std::endl
                  << O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(isym) << std::endl;
      }
      auto m = O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l);
      if (m->parity() == parityOddPacked) {
        T phase = (ksym > lsym ? 1 : m_hermiticity[1]);
        return constants.return_value = phase * m->blockMap(isym)(i, j);
      } else
        return m->blockMap(isym)(i, j);
    } else {
      if (false) {
        auto m = O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(jsym);
        std::cout << "element rectangular block 2 (" << j << "," << i << ")" << std::endl
                  << m.rows() << m.cols() << std::endl
                  << "isym=" << isym << std::endl
                  << "jsym=" << jsym << std::endl
                  << "m=blockMap(jsym):\n"
                  << m << std::endl
                  << "blockMap(isym):\n"
                  << O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(isym) << std::endl;
        std::cout << "m(0,) " << m(0, 0) << m(0, 1) << std::endl;
        std::cout << "m(1,) " << m(1, 0) << m(1, 1) << std::endl;
        std::cout << m(j, i) << std::endl;
      }
      auto m = O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l);
      if (m->parity() == parityOddPacked) {
        T phase = (ksym > lsym ? 1 : m_hermiticity[1]) * m_hermiticity[0];
        return constants.return_value = phase * m->blockMap(jsym)(j, i);
      } else
        return O2(spinUp1, spinUp2, mulliken).smat(ksym ^ lsym, ksym, k, l)->blockMap(jsym)(j, i);
    }
  }
}
template <class T>
T& Operator_<T>::element(int i, int isym, int j, int jsym, int k, int ksym, int l, int lsym, bool spinUp1,
                         bool spinUp2) {
  m_dirac_out_of_date = true;
  return const_cast<T&>(
      static_cast<const Operator_<T>*>(this)->element(i, isym, j, jsym, k, ksym, l, lsym, spinUp1, spinUp2, true));
}

template <class T> std::string Operator_<T>::str(std::string title, int level) const {
  if (level < 0)
    return "";
  if (title.empty())
    title = m_description;
  std::stringstream s;
  s << "Operator_<T> " << (title.empty() ? m_description : title);
  s << "; rank=" << m_rank;
  s << "; symmetry=" << m_symmetry + 1;
  s << "; hermiticity=" << m_hermiticity[0] << "," << m_hermiticity[0];
  if (m_rank > 1)
    s << "; exchange=" << m_exchange[0] << "," << m_exchange[0];
  s << std::endl;
  if (level > 0) {
    s << "Constant part: " << m_O0 << std::endl;
    if (m_rank > 0)
      s << O1(true).str(title + ": 1-particle, " + (m_uhf ? "alpha" : "spin-free"), level) << std::endl;
    if (m_rank > 0 && m_uhf)
      s << O1(false).str(title + ": 1-particle, beta", level) << std::endl;
    ensure_dirac();
    if (m_O2aa_mulliken != nullptr)
      s << O2(true, true, true)
               .str(title + ": 2-particle, " + (m_uhf ? "alpha-alpha" : "spin-free") + ", Mulliken", level)
        << std::endl;
    if (m_O2ab_mulliken != nullptr && m_uhf)
      s << O2(true, false, true).str(title + ": 2-particle, alpha-beta, Mulliken", level) << std::endl;
    if (m_O2bb_mulliken != nullptr && m_uhf)
      s << O2(false, false, true).str(title + ": 2-particle, beta-beta, Mulliken", level) << std::endl;
    if (m_O2aa_dirac != nullptr)
      s << O2(true, true, false).str(title + ": 2-particle, alpha-alpha, Dirac", level) << std::endl;
    if (m_O2ab_dirac != nullptr)
      s << O2(true, false, false).str(title + ": 2-particle, alpha-beta, Dirac", level) << std::endl;
    if (m_O2bb_dirac != nullptr && m_uhf)
      s << O2(false, false, false).str(title + ": 2-particle, beta-beta, Dirac", level) << std::endl;
  }
  return s.str();
}

template <class T>
Operator_<T> Operator_<T>::fock(const Operator_<T>& density, bool oneElectron, std::string description) const {
  if (m_exchange[0] >= 0 || m_exchange[1] >= 0 || m_hermiticity[0] * m_hermiticity[1] <= 0)
    throw std::logic_error("Operator_<T>::fock can work only with hermitian fermionic hamiltonian");
  if (m_dimensions[0] != m_dimensions[1])
    throw std::logic_error("Operator_<T>::fock can work only with alpha and beta spaces of same dimension");
  Operator_<T> result(m_dimensions, 1, m_uhf || density.m_uhf, std::vector<int>{1, 1}, std::vector<int>{1, 1},
                      m_symmetry ^ density.m_symmetry, true, false, description.empty() ? "Fock" : description);
  result.O1() = O1();
  if (result.m_uhf)
    result.O1(false) = O1(false);
  if (not oneElectron)
    result.O1() *= 0;
  if (not oneElectron && result.m_uhf)
    result.O1(false) *= 0;
  if (density.m_rank == 1) {
    for (unsigned int symk = 0; symk < 8; symk++) {
      unsigned int syml = symk ^ result.m_symmetry;
      if (m_hermiticity[1] && syml > symk)
        continue;
      for (size_t k = 0; k < m_dimensions[0][2][symk]; k++) {
        for (size_t l = (m_hermiticity[1] && symk == syml) ? k : 0; l < m_dimensions[0][3][syml]; l++) {
          if (result.m_uhf) {
            result.element(k, symk, l, syml, true) +=
                ((density.O1(true)) & (*O2(true, true, true).smat(result.m_symmetry, symk, k, l))) +
                ((density.O1(false)) & (*O2(false, true, true).smat(result.m_symmetry, symk, k, l)));
            result.element(k, symk, l, syml, false) +=
                ((density.O1(true)) & (*O2(true, false, true).smat(result.m_symmetry, symk, k, l))) +
                ((density.O1(false)) & (*O2(false, false, true).smat(result.m_symmetry, symk, k, l)));
          } else
            result.element(k, symk, l, syml, true) +=
                ((density.O1(true)) & (*O2(true, true, true).smat(result.m_symmetry, symk, k, l)));

          for (unsigned int symi = 0; symi < 8; symi++) {
            for (size_t i = 0; i < m_dimensions[0][0][symi]; i++) { // N^5. If it matters, this could be speeded up.
              if (result.m_uhf) {
                result.element(k, symk, l, syml, true) -=
                    (*O2(true, true).smat(symi ^ symk, symi, i, k) * density.O1(true)).blockMap(syml)(l, i);
                result.element(k, symk, l, syml, false) -=
                    (*O2(false, false).smat(symi ^ symk, symi, i, k) * density.O1(false)).blockMap(syml)(l, i);
              } else
                result.element(k, symk, l, syml, true) -=
                    0.5 * (*O2(true, true).smat(symi ^ symk, symi, i, k) * density.O1(true)).blockMap(syml)(l, i);
            }
          }
        }
      }
    }
  } else
    throw std::logic_error("Unimplemented fock of 2-particle density");
  return result;
}

template <class T> molpro::bytestream Operator_<T>::bytestream() const {
  class molpro::bytestream bs;
  for (int spin = 0; spin < 2; spin++)
    for (int i = 0; i < 4; i++)
      bs.append(&m_dimensions[spin][i][0], m_dimensions[spin][i].size());
  bs.append(m_rank);
  bs.append(m_uhf);
  bs.append(&m_hermiticity[0], 2);
  bs.append(&m_exchange[0], 2);
  bs.append(m_symmetry);
  bs.append(m_covariant ? 1 : 0);
  bs.append(&m_description[0], m_description.size());
  bs.append(m_dirac_out_of_date);
  bs.append(&m_O0, 1);
  if (m_rank > 0)
    bs.append(m_O1a->bytestream());
  if (m_rank > 0 && m_uhf)
    bs.append(m_O1b->bytestream());
  if (m_rank > 1)
    bs.append(m_O2aa_mulliken->bytestream());
  if (m_rank > 1 && m_uhf)
    bs.append(m_O2ab_mulliken->bytestream());
  if (m_rank > 1 && m_uhf)
    bs.append(m_O2bb_mulliken->bytestream());
  if (m_rank > 1 && !m_dirac_out_of_date)
    bs.append(m_O2aa_dirac->bytestream());
  if (m_rank > 1 && m_uhf && !m_dirac_out_of_date)
    bs.append(m_O2ab_dirac->bytestream());
  if (m_rank > 1 && m_uhf && !m_dirac_out_of_date)
    bs.append(m_O2bb_dirac->bytestream());
  return bs;
}

template <class T> void Operator_<T>::zero() {
  m_O0 = 0;
  m_O1a->data()->assign(0);
  if (m_rank > 1)
    m_O2aa_mulliken->data()->assign(0);
  if (m_uhf) {
    m_O1b->data()->assign(0);
    if (m_rank > 1)
      m_O2ab_mulliken->data()->assign(0);
    if (m_rank > 1)
      m_O2bb_mulliken->data()->assign(0);
  }
  m_dirac_out_of_date = true;
  ensure_dirac();
}

template <class T> void Operator_<T>::set_dirty() const { m_dirac_out_of_date = true; }

template <class T> Operator_<T>& Operator_<T>::operator-=(const Operator_<T>& other) {
  checkCompatible(other);
  if (other.m_rank > m_rank)
    throw std::logic_error("Incompatible operator ranks");
  m_O0 -= other.m_O0;
  *m_O1a -= *other.m_O1a;
  if (other.m_rank > 1)
    *m_O2aa_mulliken -= *other.m_O2aa_mulliken;
  if (m_uhf) {
    *m_O1b -= *other.m_O1b;
    if (other.m_rank > 1)
      *m_O2ab_mulliken -= *other.m_O2ab_mulliken;
    if (other.m_rank > 1)
      *m_O2bb_mulliken -= *other.m_O2bb_mulliken;
  }
  m_dirac_out_of_date = true;
  return *this;
}

template <class T> Operator_<T>& Operator_<T>::operator+=(const Operator_<T>& other) {
  checkCompatible(other);
  if (other.m_rank > m_rank)
    throw std::logic_error("Incompatible operator ranks");
  m_O0 += other.m_O0;
  *m_O1a += *other.m_O1a;
  if (other.m_rank > 1)
    *m_O2aa_mulliken += *other.m_O2aa_mulliken;
  if (m_uhf) {
    *m_O1b += *other.m_O1b;
    if (other.m_rank > 1)
      *m_O2ab_mulliken += *other.m_O2ab_mulliken;
    if (other.m_rank > 1)
      *m_O2bb_mulliken += *other.m_O2bb_mulliken;
  }
  m_dirac_out_of_date = true;
  return *this;
}

template <class T> Operator_<T>& Operator_<T>::operator*=(T other) {
  m_O0 *= other;
  *m_O1a *= other;
  if (m_rank > 1)
    *m_O2aa_mulliken *= other;
  if (m_uhf) {
    *m_O1b *= other;
    if (m_rank > 1)
      *m_O2ab_mulliken *= other;
    if (m_rank > 1)
      *m_O2bb_mulliken *= other;
  }
  m_dirac_out_of_date = true;
  return *this;
}

template <class T> typename SMat_<T>::scalar_type Operator_<T>::operator&(const Operator_<T>& other) const {
  if (other.m_rank != m_rank)
    throw std::logic_error("Incompatible operator ranks");
  T result = m_O0 * other.m_O0;
  //  std::cout << "@@ Operator_<T>::operator&"<<std::endl;
  //  std::cout << "1a: " <<this->O1(true).data()[0]<<other.O1(true).data()[0]<<(this->O1(true) &
  //  other.O1(true))<<std::endl; std::cout << "1b: "
  //  <<this->O1(true).data()[0]<<other.O1(true).data()[0]<<(this->O1(false) & other.O1(false))<<std::endl; std::cout <<
  //  "2aa: " <<this->O2(true,true).data()[0]<<other.O2(true,true).data()[0]<<
  //          (this->O2(true,true) & other.O2(true,true))<<std::endl;
  //  std::cout << "2ab: " <<this->O2(true,false).data()[0]<<other.O2(true,false).data()[0]<<
  //          (this->O2(true,false) & other.O2(true,false))<<std::endl;
  //  std::cout << "2bb: " <<this->O2(false,false).data()[0]<<other.O2(false,false).data()[0]<<
  //          (this->O2(false,false) & other.O2(false,false))<<std::endl;
  if (m_uhf) {
    if (m_rank > 0)
      result += (this->O1(true) & other.O1(true)) + (this->O1(false) & other.O1(false));
    if (m_rank > 1)
      result += (this->O2(true, true) & other.O2(true, true)) + T(2) * (this->O2(true, false) & other.O2(true, false)) +
                (this->O2(false, false) & other.O2(false, false));
  } else {
    if (m_rank > 0)
      result += T(1) * (this->O1(true) & other.O1(true));
    //      std::cout << "result after 1-elec "<<result;
    //      std::cout << this->O2(true,false).str( "this->O2(true,false)",2)<<std::endl;
    //      std::cout << other.O2(true,false).str( "other.O2(true,false)",2)<<std::endl;
    if (m_rank > 1)
      result += T(0.5) * (this->O2(true, false) & other.O2(true, false));
  }
  //      std::cout << "result after 2-elec "<<result<<std::endl;;
  return result;
}

/*!
 * \brief Construct M(pq,rs)=this.mulliken(q,p,r,s)+delta(p,r)*this(q,s).
 * The composite indices pq, rs correspond to the memory layout of this(p,q)
 * Implementation defined only when p,q,r,s are in the same space.
 * Implementation defined only for a spin-summed operator ie m_uhf=false
 * \return
 */
template <class T> SMat_<T> Operator_<T>::metric() {
  // FIXME some policing of restrictions
  dim_t pqlen;
  const auto& Gamma = O2(true, false, true);
  //  std::cout << "Gamma\n"<<Gamma<<std::endl;
  const auto& gamma = O1(true);
  //  std::cout << "gamma\n"<<gamma<<std::endl;
  for (auto pqsym = 0; pqsym < 8; pqsym++)
    pqlen.emplace_back(Gamma.block_size(pqsym));
  //  for (auto pqsym=0; pqsym<8; pqsym++) {
  //    std::cout << "2pdm ab block\n" <<Gamma.blockM(pqsym)<<std::endl;
  //    }
  //  std::cout << "metric: this\n"<<*this<<std::endl;
  //  std::cout <<"metric dimensions:";for (auto&p : pqlen) std::cout << " "<<p; std::cout <<std::endl;
  SMat_<T> result(dims_t{pqlen, pqlen}, parityNone, 0, false, m_description + " metric");
  for (auto pqsym = 0; pqsym < 8; pqsym++) {
    std::copy(Gamma.block(pqsym).begin(), Gamma.block(pqsym).end(), result.block(pqsym).begin());
    size_t offset = 0;
    for (unsigned int psym = 0; psym < 8; psym++) {
      unsigned int qsym = psym ^ pqsym;
      for (size_t p = 0; p < m_dimensions[0][0][psym]; p++)
        for (size_t q = 0; q < m_dimensions[0][0][qsym]; q++)
          for (unsigned int rsym = 0; rsym < 8; rsym++) {
            unsigned int ssym = rsym ^ pqsym;
            for (size_t r = 0; r < m_dimensions[0][0][rsym]; r++)
              for (size_t s = 0; s < m_dimensions[0][0][ssym]; s++) {
                result.block(pqsym)[offset] = Gamma.smat(pqsym, qsym, p, q)->blockMap(rsym)(r, s);
                //                        if (rsym==psym && r==p) result.block(pqsym)[offset] +=
                //                        gamma.blockMap(qsym)(q,s);
                if (ssym == qsym && s == q)
                  result.block(pqsym)[offset] += gamma.blockMap(psym)(p, r);
                if (p == q && psym == qsym)
                  result.block(pqsym)[offset] = 0;
                if (r == s && rsym == ssym)
                  result.block(pqsym)[offset] = 0;
                //                        std::cout << "p="<<p+1<<",q="<<q+1<<",r="<<r+1<<",s="<<s+1<<",
                //                        Gamma(q,p,r,s)="<<
                //                        Gamma.smat(pqsym,qsym,q,p)->blockMap(rsym)(r,s)<<"="<<element(q,qsym,p,psym,r,rsym,s,ssym,true,false)<<",
                //                        S="<<result.block(pqsym)[offset]<<std::endl;
                ++offset;
              }
          }
    }
  }
  //  std::cout <<"constructed\n"<<result<<std::endl;
  return result;
}

#endif // SYMMETRY_MATRIX_OPERATOR_IMPLEMENTATION_H
