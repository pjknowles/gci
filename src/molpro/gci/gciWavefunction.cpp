#include "gciWavefunction.h"

#include <molpro/cblas.h>

#include "gci.h"
#include "gciOrbitals.h"
#include "gciRun.h"
#include "gciStringSet.h"
#include "gciTransitionDensity.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <molpro/Profiler.h>
#include <molpro/linalg/array/util/Distribution.h>
#include <set>
#include <vector>

namespace molpro {
namespace gci {
using molpro::linalg::array::DistrArrayMPI3;

int _mpi_rank(MPI_Comm comm) {
  int r;
  MPI_Comm_rank(comm, &r);
  return r;
}

int _mpi_size(MPI_Comm comm) {
  int s;
  MPI_Comm_size(comm, &s);
  return s;
}

Wavefunction::Wavefunction(OrbitalSpace h, int nelec, int symmetry, int ms2, MPI_Comm communicator)
    : State(h, nelec, symmetry, ms2), m_sparse(false), m_communicator(communicator),
      m_parallel_rank(_mpi_rank(m_communicator)), m_parallel_size(_mpi_size(m_communicator)), dimension(0) {
  buildStrings();
}

Wavefunction::Wavefunction(OrbitalSpace *h, int n, int s, int m2, MPI_Comm communicator)
    : State(h, n, s, m2), m_sparse(false), m_communicator(communicator), m_parallel_rank(_mpi_rank(m_communicator)),
      m_parallel_size(_mpi_size(m_communicator)), dimension(0) {
  buildStrings();
}

Wavefunction::Wavefunction(const State &state, MPI_Comm communicator)
    : State(state), m_sparse(false), m_communicator(communicator), m_parallel_rank(_mpi_rank(m_communicator)),
      m_parallel_size(_mpi_size(m_communicator)), dimension(0) {
  buildStrings();
}

Wavefunction::Wavefunction(const Wavefunction &source, int option, MPI_Comm communicator)
    : m_sparse(source.m_sparse), m_communicator(communicator), m_parallel_rank(-1), m_parallel_size(0),
      dimension(source.dimension) {
  *this = source;
  if (!buffer.empty()) {
    int ranks;
    MPI_Comm_size(m_communicator,&ranks);
    auto distribution = molpro::linalg::array::util::make_distribution_spread_remainder<DistrArrayMPI3::index_type>(dimension,ranks);
    DistrArrayMPI3::index_type start, end;
    std::tie(start, end) = distribution.range(m_parallel_rank);
    distr_buffer.reset(new DistrArrayMPI3(std::make_unique<molpro::linalg::array::util::Distribution<DistrArrayMPI3::index_type>>(distribution),
                                          m_communicator,
                                          molpro::linalg::array::span::Span<double>(&buffer[start], end - start)));
  }
  if (m_communicator == MPI_COMM_NULL)
    m_communicator = source.m_communicator;
  m_parallel_rank = _mpi_rank(m_communicator);
  m_parallel_size = _mpi_size(m_communicator);
}

Wavefunction::Wavefunction(const Wavefunction &source) : Wavefunction(source, 0, source.m_communicator) {}

void Wavefunction::buildStrings() {
  profiler->start("buildStrings");
  alphaStrings.resize(8);
  betaStrings.resize(8);
  dimension = 0;
  _blockOffset.resize(8);
  for (unsigned int syma = 0; syma < 8; syma++) {
    unsigned int symb = syma ^ symmetry;
    String stringa(this, 1);
    stringa.first((nelec + ms2) / 2);
    alphaStrings[syma] = StringSet(stringa, true, syma);
    String stringb(this, -1);
    stringb.first((nelec - ms2) / 2);
    betaStrings[symb] = StringSet(stringb, true, symb);
    _blockOffset[syma] = dimension;
    dimension += alphaStrings[syma].size() * betaStrings[symb].size();
  }
  profiler->stop("buildStrings");
}

void Wavefunction::allocate_buffer() {
  if (m_sparse)
    buffer_sparse.clear();
  else {
    buffer.resize(dimension, (double)0);
    int ranks;
    MPI_Comm_size(m_communicator,&ranks);
    auto distribution = molpro::linalg::array::util::make_distribution_spread_remainder<DistrArrayMPI3::index_type>(dimension,ranks);
    DistrArrayMPI3::index_type start, end;
    std::tie(start, end) = distribution.range(m_parallel_rank);
    distr_buffer.reset(new DistrArrayMPI3(std::make_unique<molpro::linalg::array::util::Distribution<DistrArrayMPI3::index_type>>(
        distribution),
                                          m_communicator,
                                          molpro::linalg::array::span::Span<double>(&buffer[start], end - start)));
//    std::tie(start, end) = distr_buffer.distribution().range(m_parallel_rank);
//    distr_buffer.allocate_buffer({&buffer[start], end - start});
  }
}

void Wavefunction::set(size_t offset, const double val) {
  if (m_sparse)
    buffer_sparse[offset] = val;
  else
    buffer.at(offset) = val;
}

void Wavefunction::set(const double value) {
  allocate_buffer();
  for (auto &b : buffer)
    b = value;
}

void Wavefunction::diagonalOperator(const molpro::Operator &op) {
  auto p = profiler->push("diagonalOperator");
  auto ha = int1(op, true);
  auto hbb = int1(op, false);
  Eigen::MatrixXd Jaa;
  Eigen::MatrixXd Jab;
  Eigen::MatrixXd Jbb;
  Eigen::MatrixXd Kaa;
  Eigen::MatrixXd Kbb;
  if (op.m_rank > 1) {
    Jaa = intJ(op, true, true);
    Jab = intJ(op, true, false);
    Jbb = intJ(op, false, false);
    Kaa = intK(op, true);
    Kbb = intK(op, false);
  }
  size_t offset = 0;
  set(op.m_O0);
  for (unsigned int syma = 0; syma < 8; syma++) {
    unsigned int symb = syma ^ symmetry;
    size_t nsa = alphaStrings[syma].size();
    size_t nsb = betaStrings[symb].size();
    size_t nact = orbitalSpace->total();
    if (!nsa || !nsb)
      continue;
    std::vector<double> ona = alphaStrings[syma].occupationNumbers();
    std::vector<double> onb = betaStrings[symb].occupationNumbers();
    if (false && orbitalSpace->spinUnrestricted) { // UHF
    } else {                                       // RHF
      // static task distribution
      size_t chunk = (nsa - 1) / m_parallel_size + 1;
      for (size_t ia = chunk * m_parallel_rank; ia < chunk * (m_parallel_rank + 1) && ia < nsa; ia++) {
        std::vector<double> on(onb);
        for (size_t i = 0; i < nact; i++) {
          for (size_t ib = 0; ib < nsb; ib++) {
            on[ib + i * nsb] += ona[ia + i * nsa];
          }
        }
        for (size_t i = 0; i < nact; i++) {
          for (size_t ib = 0; ib < nsb; ib++) {
            buffer[offset + ia * nsb + ib] += on[ib + i * nsb] * ha[i];
          }
        }
        if (op.m_rank > 1) {
          for (size_t i = 0; i < nact; i++) {
            for (size_t j = 0; j <= i; j++) {
              double zz = Jaa(j, i) - (double)0.5 * Kaa(j, i);
              if (i == j)
                zz *= (double)0.5;
              for (size_t ib = 0; ib < nsb; ib++) {
                buffer[offset + ia * nsb + ib] += on[ib + i * nsb] * on[ib + j * nsb] * zz;
              }
            }
          }
          double vv = (double)ms2 * ms2;
          std::vector<double> f(nsb, (double)0);
          for (size_t i = 0; i < nact; i++) {
            for (size_t ib = 0; ib < nsb; ib++) {
              on[ib + i * nsb] *= ((double)2 - on[ib + i * nsb]); // on becomes a mask for singly-occupied orbitals
              f[ib] += on[ib + i * nsb];
            }
          }
          for (size_t ib = 0; ib < nsb; ib++) {
            f[ib] = f[ib] < (double)1.1 ? (double)1.1 : f[ib];
          } // mask off singularities (this code does nothing for 0 or 1 open-shell orbitals)
          for (size_t ib = 0; ib < nsb; ib++) {
            f[ib] = (vv - f[ib]) / (f[ib] * (f[ib] - (double)1));
          }
          for (size_t i = 0; i < nact; i++) {
            for (size_t j = 0; j < i; j++) {
              double zz = -(double)0.5 * Kaa(i, j);
              for (size_t ib = 0; ib < nsb; ib++) {
                buffer[offset + ia * nsb + ib] += f[ib] * on[ib + i * nsb] * on[ib + j * nsb] * zz;
              }
            }
            double zz = -(double)0.25 * Kaa(i, i);
            for (size_t ib = 0; ib < nsb; ib++) {
              buffer[offset + ia * nsb + ib] += on[ib + i * nsb] * zz;
            }
          }
        }
      }
      gather_chunks(&buffer[offset], nsa * nsb, chunk * nsb, m_communicator);
    }
    offset += nsa * nsb;
  }
  //  cout << "diagonal elements"<<std::endl; for (size_t i=0; i < buffer.size(); i++) cout <<" "<<buffer[i]; cout
  //  <<std::endl;
}

void Wavefunction::axpy(double a, const Wavefunction &xx) {
  //  cout << "Wavefunction::axpy initial=";
  //  for (size_t i=0; i<buffer.size(); i++) cout<<" "<<buffer[i]; cout << std::endl;
  //  cout << "Wavefunction::axpy x=";
  //  for (size_t i=0; i<buffer.size(); i++) cout<<" "<<xx.buffer[i]; cout << std::endl;
  for (size_t i = 0; i < buffer.size(); i++)
    buffer[i] += xx.buffer[i] * a;
  //  cout << "Wavefunction::axpy result=";
  //  for (size_t i=0; i<buffer.size(); i++) cout<<" "<<buffer[i]; cout << std::endl;
}

void Wavefunction::axpy(double a, const std::map<size_t, double> &x) {
  for (const auto &xx : x) {
    buffer[xx.first] += xx.second * a;
  }
}

std::tuple<std::vector<size_t>, std::vector<double>>
Wavefunction::select(const Wavefunction &measure, const size_t maximumNumber, const double threshold) const {
  std::multimap<double, size_t, std::greater<double>> sortlist;
  assert(buffer.size() == measure.size());
  for (size_t i = 0; i < buffer.size(); i++) {
    auto test = buffer[i] * measure.buffer[i];
    //    cout << "select "<<buffer[i]<<" "<<measure.buffer[i]<<std::endl;
    if (test > threshold) {
      sortlist.insert(std::make_pair(test, i));
      if (sortlist.size() > maximumNumber)
        sortlist.erase(std::prev(sortlist.end()));
    }
  }
  std::vector<size_t> indices;
  indices.reserve(sortlist.size());
  std::vector<double> values;
  values.reserve(sortlist.size());
  for (const auto &p : sortlist) {
    indices.push_back(p.second);
    values.push_back(p.first);
  }
  return std::make_tuple(indices, values);
}

void Wavefunction::scal(double a) {
  for (auto &b : buffer)
    b *= a;
}

void Wavefunction::fill(double a) { std::fill(begin(), end(), a); }

Wavefunction &Wavefunction::operator*=(const double &value) {
  for (auto &b : buffer)
    b *= value;
  return *this;
}

// Wavefunction& Wavefunction::operator = ( const Wavefunction &other)
//{
//    if (this != &other) {
//    if (! compatible(other)) throw std::domain_error("attempt to copy between incompatible Wavefunction objects");
//        *this = other;
//    }
//    return *this;
//}

Wavefunction &Wavefunction::operator+=(const Wavefunction &other) {
  if (!compatible(other))
    throw std::domain_error("attempt to add incompatible Wavefunction objects");
  //  cout << "Wavefunction::operator += &this=" << this <<" &other="<<&other <<std::endl;
  for (size_t i = 0; i < buffer.size(); i++)
    buffer[i] += other.buffer[i];
  return *this;
}

Wavefunction &Wavefunction::operator-=(const Wavefunction &other) {
  if (!compatible(other))
    throw std::domain_error("attempt to add incompatible Wavefunction objects");
  for (size_t i = 0; i < buffer.size(); i++) {
    buffer[i] -= other.buffer[i];
  }
  return *this;
}

Wavefunction &Wavefunction::operator+=(const double other) {
  for (auto &b : buffer)
    b += other;
  return *this;
}

Wavefunction &Wavefunction::operator-=(const double other) {
  for (auto &b : buffer)
    b -= other;
  return *this;
}

Wavefunction &Wavefunction::operator/=(const Wavefunction &other) {
  if (!compatible(other))
    throw std::domain_error("attempt to add incompatible Wavefunction objects");
  for (size_t i = 0; i < buffer.size(); i++) {
    buffer[i] /= other.buffer[i];
  }
  return *this;
}

double Wavefunction::update(const Wavefunction &diagonalH, double &eTruncated, double const dEmax) {
  if (!compatible(diagonalH))
    throw std::domain_error("attempt to combine incompatible Wavefunction objects");
  //  cout << "Wavefunction::update this="; for (size_t i=0; i<buffer.size(); i++) cout<<buffer[i]<<" ";cout
  //  <<std::endl; cout << "Wavefunction::update diagonalH="; for (size_t i=0; i<diagonalH.buffer.size(); i++)
  //  cout<<diagonalH.buffer[i]<<" ";cout <<std::endl;
  size_t imin = 0;
  size_t imax = buffer.size();
  eTruncated = 0;
  double ePredicted = 0;
  for (size_t i = imin; i < imax; i++) {
    double dE = buffer[i] * buffer[i] / diagonalH.buffer[i];
    if (diagonalH.buffer[i] > 1e-5)
      ePredicted -= dE;
    if (dE > dEmax)
      buffer[i] = -buffer[i] / diagonalH.buffer[i];
    else {
      buffer[i] = (double)0;
      eTruncated += dE;
    }
  }
  return ePredicted;
}

Wavefunction operator+(const Wavefunction &w1, const Wavefunction &w2) {
  Wavefunction result = w1;
  return result += w2;
}

Wavefunction operator-(const Wavefunction &w1, const Wavefunction &w2) {
  Wavefunction result = w1;
  return result -= w2;
}

Wavefunction operator/(const Wavefunction &w1, const Wavefunction &w2) {
  Wavefunction result = w1;
  return result /= w2;
}

Wavefunction operator*(const Wavefunction &w1, const double &value) {
  Wavefunction result = w1;
  return result *= value;
}

Wavefunction operator*(const double &value, const Wavefunction &w1) {
  Wavefunction result = w1;
  return result *= value;
}

double Wavefunction::dot(const Wavefunction &other) const {
  return (*this) * ((dynamic_cast<const Wavefunction &>(other)));
}

double Wavefunction::dot(const std::map<size_t, double> &other) const {
  double result = 0;
  if (m_sparse)
    for (const auto &o : other) {
      if (buffer_sparse.count(o.first))
        result += buffer_sparse.at(o.first) * o.second;
    }
  else
    for (const auto &o : other) {
      result += buffer[o.first] * o.second;
    }
  return result;
}

void Wavefunction::times(const Wavefunction *a, const Wavefunction *b) {
  const auto wa = dynamic_cast<const Wavefunction *>(a);
  const auto wb = dynamic_cast<const Wavefunction *>(b);
  for (size_t i = 0; i < buffer.size(); i++) {
    buffer[i] = wa->buffer[i] * wb->buffer[i];
  }
}

void Wavefunction::divide(const Wavefunction *a, const Wavefunction *b, double shift, bool append, bool negative) {
  const auto wa = dynamic_cast<const Wavefunction *>(a);
  const auto wb = dynamic_cast<const Wavefunction *>(b);
  if (negative) {
    if (append) {
      for (size_t i = 0; i < buffer.size(); i++) {
        buffer[i] -= wa->buffer[i] / (wb->buffer[i] + shift);
      }
    } else {
      for (size_t i = 0; i < buffer.size(); i++) {
        buffer[i] = -wa->buffer[i] / (wb->buffer[i] + shift);
      }
    }
  } else {
    if (append) {
      for (size_t i = 0; i < buffer.size(); i++) {
        buffer[i] += wa->buffer[i] / (wb->buffer[i] + shift);
      }
    } else {
      for (size_t i = 0; i < buffer.size(); i++) {
        buffer[i] = wa->buffer[i] / (wb->buffer[i] + shift);
      }
    }
  }
}

void Wavefunction::zero() {
  for (auto &b : buffer) {
    b = 0;
  }
}

double operator*(const Wavefunction &w1, const Wavefunction &w2) {
  if (!w1.compatible(w2))
    throw std::domain_error("attempt to form scalar product between incompatible Wavefunction objects");
  auto result = (double)0;
  for (size_t i = 0; i < w1.buffer.size(); i++) {
    result += w1.buffer[i] * w2.buffer[i];
  }
  return result;
}

Wavefunction &Wavefunction::operator-() {
  for (auto &b : buffer)
    b = -b;
  return *this;
}

std::string Wavefunction::values() const {
  std::ostringstream s;
  for (auto &v : buffer)
    s << " " << v;
  return s.str();
}

std::string Wavefunction::str(int verbosity, unsigned int columns) const {
  std::ostringstream s;
  if (verbosity >= 0) {
    if (verbosity >= 1)
      s << std::endl << "values at address " << &buffer << " and of length " << buffer.size();
    size_t address = 0;
    for (unsigned int syma = 0; syma < 8; syma++) {
      unsigned int symb = syma ^ symmetry;
      if (!alphaStrings[syma].empty() && !betaStrings[symb].empty()) {
        if (verbosity >= 1)
          s << std::endl << "Alpha strings of symmetry " << syma + 1 << ":";
        if (verbosity >= 2) {
          for (const auto &i : alphaStrings[syma]) {
            s << std::endl << i.str();
          }
          for (const auto &i : betaStrings[symb]) {
            s << std::endl << i.str();
          }
        }
        if (!m_sparse && buffer.size() == dimension && verbosity >= 0) {
          s << std::endl << "Values:";
          for (size_t i = 0; i < alphaStrings[syma].size(); i++) {
            s << std::endl;
            for (size_t j = 0; j < betaStrings[symb].size(); j++) {
              s << buffer[address++] << " ";
            }
          }
        }
      }
    }
    if (m_sparse && !buffer_sparse.empty() && verbosity >= 0) {
      s << std::endl << "Sparse storage, values:";
      for (const auto &b : buffer_sparse) {
        s << std::endl << determinantAt(b.first) << " : " << b.second;
      }
    }
  }
  return "Wavefunction object\n" + this->State::str(verbosity) + s.str();
}

bool Wavefunction::compatible(const Wavefunction &other) const {
  return dimension == other.dimension && buffer.size() == other.buffer.size();
}

std::vector<size_t> Wavefunction::minlocN(size_t n) const {
  std::vector<size_t> results;
  for (size_t k = 0; k < n; k++) {
    auto m = cbegin();
    while (std::count(results.begin(), results.end(), m - cbegin()) != 0)
      m++;
    for (auto i = cbegin(); i != cend(); i++) {
      if (*i < *m && std::count(results.begin(), results.end(), i - cbegin()) == 0)
        m = i;
    }
    results.push_back(m - cbegin());
  }
  return results;
}

size_t Wavefunction::minloc(size_t n) const { return minlocN(n).back(); }

double Wavefunction::at(size_t offset) const { return buffer.at(offset); }

Determinant Wavefunction::determinantAt(size_t offset) const {
  size_t address = 0;
  for (unsigned int syma = 0; syma < 8; syma++) {
    unsigned int symb = syma ^ symmetry;
    size_t newaddress = address + alphaStrings[syma].size() * betaStrings[symb].size();
    if (offset >= address && offset < newaddress) {
      size_t a = (offset - address) / betaStrings[symb].size();
      size_t b = offset - address - a * betaStrings[symb].size();
      return Determinant(this, &alphaStrings[syma][a], &betaStrings[symb][b]);
    }
    address = newaddress;
  }
  throw std::runtime_error("Wavefunction::determinantAt cannot find");
}

size_t Wavefunction::stringAddress(size_t offset, unsigned int axis) const {
  size_t address = 0;
  for (unsigned int syma = 0; syma < 8; syma++) {
    unsigned int symb = syma ^ symmetry;
    size_t newaddress = address + alphaStrings[syma].size() * betaStrings[symb].size();
    if (offset >= address && offset < newaddress) {
      size_t a = (offset - address) / betaStrings[symb].size();
      size_t b = offset - address - a * betaStrings[symb].size();
      if (axis == 0)
        return b;
      if (axis == 1)
        return a;
      throw std::runtime_error("Wavefunction::stringAddress bad axis");
    }
    address = newaddress;
  }
  throw std::runtime_error("Wavefunction::stringAddress cannot find");
}

size_t Wavefunction::stringSymmetry(size_t offset, unsigned int axis) const {
  size_t address = 0;
  for (unsigned int syma = 0; syma < 8; syma++) {
    unsigned int symb = syma ^ symmetry;
    size_t newaddress = address + alphaStrings[syma].size() * betaStrings[symb].size();
    if (offset >= address && offset < newaddress) {
      if (axis == 0)
        return symb;
      if (axis == 1)
        return syma;
      throw std::runtime_error("Wavefunction::stringSymmetry bad axis");
    }
    address = newaddress;
  }
  throw std::runtime_error("Wavefunction::stringAddress cannot find");
}

size_t Wavefunction::blockOffset(unsigned int syma) const { return _blockOffset.at(syma); }

using uint = unsigned int;
void MXM(double *Out, const double *A, const double *B, uint nRows, uint nLink, uint nCols, bool AddToDest,
         int nStrideLink = -1) {
  const bool debug = false;
  auto prof = profiler->push("MXM");
  prof += 2 * nLink * size_t(nCols) * nRows;
  if (nStrideLink < 0) {
    if (debug && AddToDest)
      cout << "MXM initial Out\n" << Eigen::Map<Eigen::MatrixXd>(Out, nRows, nCols) << std::endl;
    if (debug)
      cout << "MXM A\n" << Eigen::Map<const Eigen::MatrixXd>(A, nRows, nLink) << std::endl;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nRows, nCols, nLink, 1, A, nRows, B, nLink,
                AddToDest ? 1 : 0, Out, nRows);
  } else {
    //      if (debug) // how to make const and Stride work together?
    //        cout << "MXM A\n"<<Eigen::Map<const Eigen::MatrixXd>(A,nRows,nLink,
    //        Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1, -nStrideLink))<<std::endl;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, nRows, nCols, nLink, 1, A, static_cast<uint>(nStrideLink), B,
                nLink, AddToDest ? 1 : 0, Out, nRows);
  }
  if (debug) {
    cout << "MXM B\n" << Eigen::Map<const Eigen::MatrixXd>(B, nLink, nCols) << std::endl;
    cout << "MXM Out\n" << Eigen::Map<Eigen::MatrixXd>(Out, nRows, nCols) << std::endl;
  }
}

void Wavefunction::operatorOnWavefunction(
    const molpro::Operator &h, const Wavefunction &w,
    bool parallel_stringset) { // FIXME not really thoroughly checked if the symmetry of h is not zero.
  auto prof = profiler->push("operatorOnWavefunction" + std::string(m_sparse ? "/sparse" : "") +
                             std::string(w.m_sparse ? "/sparse" : ""));
  if (m_sparse)
    buffer_sparse.clear();
  if (m_parallel_rank == 0) {
    if (m_sparse) {
      if (!w.m_sparse)
        throw std::runtime_error("Cannot make sparse residual from full vector");
      for (const auto &b : w.buffer_sparse) {
        buffer_sparse[b.first] += h.m_O0 * b.second;
      }
    } else if (w.m_sparse) {
      for (const auto &b : w.buffer_sparse) {
        buffer[b.first] += h.m_O0 * b.second;
      }
    } else {
      for (size_t i = 0; i < buffer.size(); i++) {
        buffer[i] += h.m_O0 * w.buffer[i];
      }
    }
  } else if (!m_sparse)
    for (auto &b : buffer)
      b = 0;
  //  cout <<"residual after 0-electron:"<<std::endl<<str(2)<<std::endl;

  //  cout <<std::endl<<"w in operatorOnWavefunction="<<w.str(2)<<std::endl;
  DivideTasks(99999999, 1, 1, m_communicator);
  const auto alphaActiveStrings = w.activeStrings(true);
  const auto betaActiveStrings = w.activeStrings(false);
  //  cout << "betaActiveStrings"<<std::endl;
  //  for (const auto& s : betaActiveStrings) for (const auto& ss : s) cout <<ss<<std::endl;

  if (true) {
    auto p = profiler->push("1-electron RI");
    size_t nsaaMax = 1000000000;
    size_t nsbbMax = 1000000000;
    std::vector<StringSet> bbs;
    //    cout << "before bbs emplace "<<std::endl;
    for (unsigned int symb = 0; symb < 8; symb++) {
      //        bbs.emplace_back(w.betaStrings,1,0,symb,parallel_stringset);
      bbs.emplace_back(betaActiveStrings, 1, 0, symb, parallel_stringset);
    }
    //    cout << "after bbs emplace "<<std::endl;
    for (unsigned int syma = 0; syma < 8; syma++) {
      //          StringSet aa(w.alphaStrings,1,0,syma,parallel_stringset);
      StringSet aa(alphaActiveStrings, 1, 0, syma, parallel_stringset);
      //      cout << "after making StringSet aa "<<std::endl;
      for (unsigned int symb = 0; symb < 8; symb++) { // symmetry of N-1 electron state
        unsigned int symexc = w.symmetry ^ syma ^ symb;
        //        cout << "syma="<<syma<<" symb="<<symb<<" symexc="<<symexc<<std::endl;
        if (m_tilesize > 0) {
          nsaaMax = m_tilesize / double(betaStrings[symb].size()) + 1;
          nsbbMax = m_tilesize / double(alphaStrings[symb].size()) + 1;
        }
        if (m_alphatilesize > 0 && m_betatilesize > 0) {
          nsaaMax = m_alphatilesize;
          nsbbMax = m_betatilesize;
        }
        if (!aa.empty() && !betaStrings[symb].empty()) {
          auto ham = h.O1(true).blockCopy(symexc);
          if (ham.cols() < 1)
            continue;
          //        for (const auto& aaa: aa) cout <<"N-1 alpha member "<<aaa<<std::endl;
          for (StringSet::const_iterator aa1, aa0 = aa.begin();
               aa1 = aa0 + nsaaMax > aa.end() ? aa.end() : aa0 + nsaaMax, aa0 < aa.end();
               aa0 = aa1) { // loop over alpha batches
            if (!NextTask(m_communicator))
              continue;
            TransitionDensity d(w, aa0, aa1, w.betaStrings[symb].begin(), w.betaStrings[symb].end(), molpro::parityEven,
                                true, false);
            //            cout << "alpha transition density" << d << "\n" << w.betaStrings[symb].size() << aa1 - aa0 <<
            //            " "
            //                 << d.size() << std::endl;
            TransitionDensity e(d, false);
            //            cout << "hamiltonian block\n" << ham << std::endl;
            //                        cout << "ham dimensions "<<ham.rows()<<" "<<ham.cols()<<std::endl;
            MXM(&e[0], &d[0], &ham(0, 0), std::distance(aa0, aa1) * w.betaStrings[symb].size(), ham.rows(), ham.cols(),
                false);
            //            cout << "alpha e" << e << std::endl;
            e.action(*this);
            //            cout << "this after action " << values() << std::endl;
          }
        }
        const auto &bb = bbs[symb];
        if (!bb.empty() && !alphaStrings[syma].empty()) {
          auto ham = h.O1(false).blockCopy(symexc);
          if (ham.cols() < 1)
            continue;
          for (StringSet::const_iterator bb1, bb0 = bb.begin();
               bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
               bb0 = bb1) { // loop over beta batches
            if (!NextTask(m_communicator))
              continue;
            TransitionDensity d(w, w.alphaStrings[syma].begin(), w.alphaStrings[syma].end(), bb0, bb1,
                                molpro::parityEven, false, true);
            //            cout << "beta transition density"<<d<<std::endl;
            TransitionDensity e(d, false);
            MXM(&e[0], &d[0], &ham(0, 0), std::distance(bb0, bb1) * w.alphaStrings[syma].size(), ham.rows(), ham.cols(),
                false);
            //            cout << "beta e"<<e<<std::endl;
            e.action(*this);
            //            cout << "this after action "<<values()<<std::endl;
          }
        }
      }
    }

  } else {
    auto p = profiler->push("1-electron");
    auto tilesize = m_tilesize;
    auto alphatilesize = m_alphatilesize;
    auto betatilesize = m_betatilesize;
    //    cout << "object tile sizes "<< m_tilesize<<","<<m_alphatilesize<<","<<m_betatilesize<<std::endl;
    size_t nsaaMax = 1000000000;
    size_t nsbbMax = 1000000000;
    size_t offset = 0, nsa = 0, nsb = 0;
    for (unsigned int syma = 0; syma < 8; syma++) {
      unsigned int symb = w.symmetry ^ syma;
      if (alphaStrings[syma].empty() || betaStrings[symb].empty())
        continue;
      auto &aa = w.alphaStrings[syma];
      auto &bb = w.betaStrings[symb];
      if (tilesize > 0)
        nsaaMax = tilesize / double(bb.size()) + 1;
      if (alphatilesize > 0 && betatilesize > 0) {
        nsaaMax = alphatilesize;
        nsbbMax = betatilesize;
      }
      for (StringSet::const_iterator aa1, aa0 = aa.begin();
           aa1 = aa0 + nsaaMax > aa.end() ? aa.end() : aa0 + nsaaMax, aa0 < aa.end();
           aa0 = aa1) { // loop over alpha batches
        for (StringSet::const_iterator bb1, bb0 = bb.begin();
             bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
             bb0 = bb1) { // loop over beta batches
          nsa = aa1 - aa0;
          nsb = bb1 - bb0;
          if (NextTask(m_communicator)) {
            //                    cout << "offset="<<offset <<", nsa="<<nsa <<", nsb="<<nsb<<std::endl;
            //        cout << bb0-bb.begin()<<std::endl;
            {
              TransitionDensity d(w, aa0, aa1, bb0, bb1, molpro::parityEven, true, !h.m_uhf);
              if (nsb == bb.size())
                MXM(&buffer[offset], &d[0], &(*h.O1(true).data())[0], nsa * nsb, w.orbitalSpace->total(0, 1), 1, true);
              else {
                molpro::vector<double> result(nsa * nsb);
                MXM(result.data(), &d[0], &(*h.O1(true).data())[0], nsa * nsb, w.orbitalSpace->total(0, 1), 1, false);
                for (size_t ia = 0; ia < nsa; ia++) {
                  for (size_t ib = 0; ib < nsb; ib++) {
                    buffer[offset + ia * bb.size() + bb0 - bb.begin() + ib] += result[ia * nsb + ib];
                  }
                }
              }
            }
            if (h.m_uhf) {
              TransitionDensity d(w, aa0, aa1, bb0, bb1, molpro::parityEven, false, true);
              MXM(&buffer[offset], &d[0], &(*h.O1(false).data())[0], nsa * nsb, w.orbitalSpace->total(0, 1), 1, true);
            }
          }
        }
        offset += nsa * bb.size();
      }
    }
  }
  //  cout <<"residual after 1-electron:"<<values()<<std::endl;
  //  cout <<"residual after 1-electron:"<<std::endl<<str(2)<<std::endl;

  if (h.m_rank > 1) { // two-electron contribution, alpha-alpha
    auto p = profiler->push("aa integrals");
    size_t nsbbMax = 64; // temporary static
    for (unsigned int syma = 0; syma < 8; syma++) {
      profiler->start("StringSet aa");
      StringSet aa(alphaActiveStrings, 2, 0, syma, parallel_stringset);
      profiler->stop("StringSet aa");
      //      cout << "number of alpha-alpha-excited strings=" << aa.size() << std::endl;
      if (aa.empty())
        continue;
      for (unsigned int symb = 0; symb < 8; symb++) {
        if (!NextTask(m_communicator))
          continue;
        auto praa = profiler->push("aa1 loop");
        unsigned int symexc = syma ^ symb ^ w.symmetry;
        // size_t nexc = h.pairSpace.find(-1)->second[symexc];
        size_t nexc = h.O2(true, true, false).block_size(symexc);
        size_t nsb = betaStrings[symb].size();
        if (nsb == 0)
          continue;
        for (StringSet::iterator aa1, aa0 = aa.begin();
             aa1 = aa0 + nsbbMax > aa.end() ? aa.end() : aa0 + nsbbMax, aa0 < aa.end();
             aa0 = aa1) { // loop over alpha batches
          size_t nsa = aa1 - aa0;
          TransitionDensity d(w, aa0, aa1, w.betaStrings[symb].begin(), w.betaStrings[symb].end(),
                              molpro::parityOddPacked, false, false);
          TransitionDensity e(d, false);
          MXM(&e[0], &d[0], &h.O2(true, true, false).block(symexc)[0], nsa * nsb, nexc, nexc, false);
          e.action(*this);
        }
      }
    }
  }
  //  cout <<"residual after alpha-alpha on process "<<m_parallel_rank<<" "<<buffer[0]<<std::endl<<str(2)<<std::endl;

  if (h.m_rank > 1)
    if (true || h.m_uhf) { // two-electron contribution, beta-beta
      auto p = profiler->push("bb integrals");
      size_t nsbbMax = 64; // temporary static
      for (unsigned int symb = 0; symb < 8; symb++) {
        StringSet bb(betaActiveStrings, 2, 0, symb, parallel_stringset);
        if (bb.empty())
          continue;
        for (unsigned int syma = 0; syma < 8; syma++) {
          if (!NextTask(m_communicator))
            continue;
          unsigned int symexc = symb ^ syma ^ w.symmetry;
          size_t nexc = h.O2(false, false, false).block_size(symexc);
          size_t nsa = alphaStrings[syma].size();
          if (nsa == 0)
            continue;
          for (StringSet::iterator bb1, bb0 = bb.begin();
               bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
               bb0 = bb1) { // loop over beta batches
            size_t nsb = bb1 - bb0;
            TransitionDensity d(w, w.alphaStrings[syma].begin(), w.alphaStrings[syma].end(), bb0, bb1,
                                molpro::parityOddPacked, false, false);
            TransitionDensity e(d, false);
            MXM(&e[0], &d[0], &h.O2(false, false, false).block(symexc)[0], nsa * nsb, nexc, nexc, false);
            e.action(*this);
          }
        }
      }
      //      cout <<"residual after beta-beta on process "<<m_parallel_rank<<"
      //      "<<buffer[0]<<std::endl<<str(2)<<std::endl;
    }

  if (h.m_rank > 1) { // two-electron contribution, alpha-beta
    auto p = profiler->push("ab integrals");
    for (unsigned int symb = 0; symb < 8; symb++) {
      StringSet bb(betaActiveStrings, 1, 0, symb, parallel_stringset);
      if (bb.empty())
        continue;
      size_t nsbbMax = std::min((size_t)256, bb.size() / m_parallel_size); // temporary static
      nsbbMax = std::max(nsbbMax, (size_t)16);
      for (unsigned int syma = 0; syma < 8; syma++) {
        StringSet aa(alphaActiveStrings, 1, 0, syma, parallel_stringset);
        if (aa.empty())
          continue;
        size_t nsaaMax = std::min((size_t)256, aa.size() / m_parallel_size); // temporary static
        nsaaMax = std::max(nsaaMax, (size_t)16);
        unsigned int symexc = symb ^ syma ^ w.symmetry;
        size_t nexc = h.O2(true, false, false).block_size(symexc);
        {
          //                        cout << "syma="<<syma<<", symb="<<symb<<", symexc="<<symexc<<std::endl;
          auto pro = profiler->push("StringSet iterator loops");
          for (StringSet::iterator aa1, aa0 = aa.begin();
               aa1 = aa0 + nsaaMax > aa.end() ? aa.end() : aa0 + nsaaMax, aa0 < aa.end();
               aa0 = aa1) { // loop over alpha batches
            size_t nsa = aa1 - aa0;
            for (StringSet::iterator bb1, bb0 = bb.begin();
                 bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
                 bb0 = bb1) { // loop over beta batches
              size_t nsb = bb1 - bb0;
              if (!NextTask(m_communicator))
                continue;
              TransitionDensity d(w, aa0, aa1, bb0, bb1, molpro::parityNone, false, false);
              TransitionDensity e(d, false);
              MXM(&e[0], &d[0], &h.O2(true, false, false).block(symexc)[0], nsa * nsb, nexc, nexc, false);
              e.action(*this);
            }
          }
        }
      }
    }
    //      cout <<"residual after alpha-beta on process "<<m_parallel_rank<<"
    //      "<<buffer[0]<<std::endl<<str(2)<<std::endl;
  }

  EndTasks(m_communicator);

  if (m_sparse) {
    gsum(buffer_sparse, m_communicator);
  } else {
    gsum(&buffer[0], buffer.size(), m_communicator);
  }
}

molpro::Operator Wavefunction::density(int rank, bool uhf, bool hermitian, const Wavefunction *bra,
                                       std::string description, bool parallel_stringset) const {
  if (bra == nullptr)
    bra = this;
  auto prof = profiler->push("density");

  molpro::dim_t dim;
  for (const auto s : *orbitalSpace)
    dim.push_back(s);
  molpro::Operator result(molpro::dims_t{dim, dim, dim, dim}, rank, uhf,
                          (hermitian ? std::vector<int>{1, 1} : std::vector<int>{0, 0}), {-1, -1},
                          symmetry ^ bra->symmetry, false, false, description);
  //  std::cout << "result\n"<<result.str("result",3)<<std::endl;
  result.zero();
  result.m_O0 = (*this) * (*bra);
  //  std::cout << "@@@ density after zero before construct\n"<<result.str("result",3)<<std::endl;

  DivideTasks(99999999, 1, 1, m_communicator);

  {
    auto p = profiler->push("1-electron");
    size_t offset = 0, nsa = 0, nsb = 0;
    for (unsigned int syma = 0; syma < 8; syma++) {
      offset += nsa * nsb;
      unsigned int symb = symmetry ^ syma;
      nsa = alphaStrings[syma].size();
      nsb = betaStrings[symb].size();
      if (!NextTask(m_communicator) || nsa * nsb == 0)
        continue;
      {
        TransitionDensity d(*this, alphaStrings[syma].begin(), alphaStrings[syma].end(), betaStrings[symb].begin(),
                            betaStrings[symb].end(), (hermitian ? molpro::parityEven : molpro::parityNone), true,
                            !result.m_uhf);
        //      cout << "Transition density\n"<<d<<std::endl;
        MXM(&((*result.O1(true).data())[0]), &d[0], &(bra->buffer[offset]),
            orbitalSpace->total(0, (hermitian ? molpro::parityEven : molpro::parityNone)), nsa * nsb, 1, true,
            nsa * nsb);
      }
      if (result.m_uhf) {
        TransitionDensity d(*this, alphaStrings[syma].begin(), alphaStrings[syma].end(), betaStrings[symb].begin(),
                            betaStrings[symb].end(), (hermitian ? molpro::parityEven : molpro::parityNone), false,
                            true);
        MXM(&((*result.O1(false).data())[0]), &d[0], &(bra->buffer[offset]),
            orbitalSpace->total(0, (hermitian ? molpro::parityEven : molpro::parityNone)), nsa * nsb, 1, true,
            nsa * nsb);
      }
    }
    std::vector<bool> spincases(1, true);
    if (result.m_uhf)
      spincases.push_back(false);
    for (auto spincase : spincases) {
      if (result.O1(spincase).parity() == molpro::parityEven) {
        if (result.O1(spincase).symmetry() == 0) {
          //          result.O1(spincase) *=.5;
          //          for (uint isym=0; isym<8; isym++)
          //            for (uint i=0; i<result.O1(spincase).dimension(isym); i++)
          //              result.O1(spincase).block(isym)[(i+2)*(i+1)/2-1]*=2;
          result.O1(spincase).scal(0.5, false);
        } else {
          throw std::runtime_error("non-symmetric density not yet supported");
        }
      }
    }
  }

  if (rank > 1) { // two-electron contribution, alpha-alpha
    //  std::cout << "@@@ density before construct 2 elec\n"<<result.str("result before construct 2e",3)<<std::endl;
    auto p = profiler->push("aa density");
    size_t nsbbMax = 64; // temporary static
    for (unsigned int syma = 0; syma < 8; syma++) {
      profiler->start("StringSet aa");
      StringSet aa(alphaStrings, 2, 0, syma, parallel_stringset);
      profiler->stop("StringSet aa");
      if (aa.empty())
        continue;
      for (unsigned int symb = 0; symb < 8; symb++) {
        if (!NextTask(m_communicator))
          continue;
        auto praa = profiler->push("aa1 loop");
        unsigned int symexc = syma ^ symb ^ symmetry;
        // size_t nexc = h.pairSpace.find(-1)->second[symexc];
        size_t nexc = result.O2(true, true, false).block_size(symexc);
        size_t nsb = betaStrings[symb].size();
        if (nsb == 0)
          continue;
        for (StringSet::iterator aa1, aa0 = aa.begin();
             aa1 = aa0 + nsbbMax > aa.end() ? aa.end() : aa0 + nsbbMax, aa0 < aa.end();
             aa0 = aa1) { // loop over alpha batches
          size_t nsa = aa1 - aa0;
          TransitionDensity d(*this, aa0, aa1, betaStrings[symb].begin(), betaStrings[symb].end(),
                              molpro::parityOddPacked, false, false);
          TransitionDensity e(*bra, aa0, aa1, betaStrings[symb].begin(), betaStrings[symb].end(),
                              molpro::parityOddPacked, false, false);
          //                      cout <<"Daa: "<<d.str(2)<<std::endl;
          //                      cout <<"Eaa: "<<e.str(2)<<std::endl;
          //                cout << "result goes to "<<&result.O2(true,true).block(symexc)[0]<<std::endl;
          //                cout << "result goes to "<<&result.O2(true,true).blockM(symexc)(0,0)<<std::endl;
          //                cout << "initial result "<<result.O2(true,true).blockM(symexc)(0,0)<<std::endl;
          MXM(&result.O2(true, true, false).block(symexc)[0], &d[0], &e[0], nexc, nsa * nsb, nexc, true,
              nsa * nsb); // not yet right for non-symmetric tdm
          //                cout << "final result "<<result.O2(true,true,false).blockM(symexc)(0,0)<<std::endl;
          //                      cout << "result"<<result.O2(true,true,false).blockM(symexc)<<std::endl;
        }
      }
    }
  }

  if (rank > 1 && result.m_uhf) { // two-electron contribution, beta-beta
    auto p = profiler->push("bb density");
    size_t nsbbMax = 64; // temporary static
    for (unsigned int symb = 0; symb < 8; symb++) {
      StringSet bb(betaStrings, 2, 0, symb, parallel_stringset);
      if (bb.empty())
        continue;
      for (unsigned int syma = 0; syma < 8; syma++) {
        if (!NextTask(m_communicator))
          continue;
        unsigned int symexc = symb ^ syma ^ symmetry;
        size_t nexc = result.O2(false, false, false).block_size(symexc);
        size_t nsa = alphaStrings[syma].size();
        if (nsa == 0)
          continue;
        for (StringSet::iterator bb1, bb0 = bb.begin();
             bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
             bb0 = bb1) { // loop over beta batches
          size_t nsb = bb1 - bb0;
          TransitionDensity d(*this, alphaStrings[syma].begin(), alphaStrings[syma].end(), bb0, bb1,
                              molpro::parityOddPacked, false, false);
          TransitionDensity e(*bra, alphaStrings[syma].begin(), alphaStrings[syma].end(), bb0, bb1,
                              molpro::parityOddPacked, false, false);
          MXM(&result.O2(false, false, false).block(symexc)[0], &d[0], &e[0], nexc, nsa * nsb, nexc, true,
              nsa * nsb); // not yet right for non-symmetric tdm
        }
      }
    }
  }

  if (rank > 1) { // two-electron contribution, alpha-beta
    auto p = profiler->push("ab density");
    size_t nsaaMax = 128; // temporary static
    size_t nsbbMax = 128; // temporary static
    for (unsigned int symb = 0; symb < 8; symb++) {
      StringSet bb(betaStrings, 1, 0, symb, parallel_stringset);
      if (bb.empty())
        continue;
      for (unsigned int syma = 0; syma < 8; syma++) {
        StringSet aa(alphaStrings, 1, 0, syma, parallel_stringset);
        if (aa.empty())
          continue;
        unsigned int symexc = symb ^ syma ^ symmetry;
        size_t nexc = result.O2(true, false, false).block_size(symexc);
        {
          auto pro = profiler->push("StringSet iterator loops");
          for (StringSet::iterator aa1, aa0 = aa.begin();
               aa1 = aa0 + nsaaMax > aa.end() ? aa.end() : aa0 + nsaaMax, aa0 < aa.end();
               aa0 = aa1) { // loop over alpha batches
            size_t nsa = aa1 - aa0;
            for (StringSet::iterator bb1, bb0 = bb.begin();
                 bb1 = bb0 + nsbbMax > bb.end() ? bb.end() : bb0 + nsbbMax, bb0 < bb.end();
                 bb0 = bb1) { // loop over beta batches
              size_t nsb = bb1 - bb0;
              if (!NextTask(m_communicator))
                continue;
              TransitionDensity d(*this, aa0, aa1, bb0, bb1, molpro::parityNone, false, false);
              TransitionDensity e(*bra, aa0, aa1, bb0, bb1, molpro::parityNone, false, false);
              //                      cout <<"D: "<<d.str(2)<<std::endl;
              //                      cout <<"E: "<<e.str(2)<<std::endl;
              //                      cout << "ab result before
              //                      MXM"<<result.O2(true,false,false).str("result",2)<<std::endl; cout << "dimensions
              //                      "<<nexc<<nsa*nsb<<std::endl;
              MXM(&result.O2(true, false, false).block(symexc)[0], &d[0], &e[0], nexc, nsa * nsb, nexc, true,
                  nsa * nsb); // not yet right for non-symmetric tdm
              //                      cout << "ab result"<<result.O2(true,false,false).str("result",2)<<std::endl;
            }
          }
        }
      }
    }
  }

  EndTasks(m_communicator);

  gsum(result);
  //  std::cout << "Density before from_dirac:\n"<<result<<std::endl;
  result.mulliken_from_dirac();
  if (result.m_rank > 1 && result.m_hermiticity == std::vector<int>{1, 1}) {
    result.O2(true, true).scal(0.5, false);
    if (result.m_uhf) {
      result.O2(true, false).scal(0.5, false);
      result.O2(false, false).scal(0.5, false);
    }
  }
  //  std::cout << "Density:\n"<<result<<std::endl;
  return result;
}

Orbitals Wavefunction::naturalOrbitals() {
  Orbitals orb(*orbitalSpace);
  molpro::SMat dens1 = density(1).O1();
  dens1.ev(orb.m_occupations, &orb.m_orbitals, nullptr, nullptr, "lapack", "descending");
  orb.m_orbitals.m_description = "Natural orbitals";
  orb.m_occupations.m_description = "Natural orbital occupation numbers";
  return orb;
}

void Wavefunction::putw(File &f, int index) {
  auto p = profiler->push("Wavefunction::putw");
  size_t chunk = (buffer.size() - 1) / m_parallel_size + 1;
  size_t offset = chunk * m_parallel_rank;
  if (chunk + offset > buffer.size())
    chunk = buffer.size() - offset;
  if (offset < buffer.size())
    f.write(&buffer[offset], chunk, index * chunk);
  p += chunk;
}

void Wavefunction::getw(File &f, int index) {
  auto p = profiler->push("Wavefunction::getw");
  size_t chunk = (buffer.size() - 1) / m_parallel_size + 1;
  size_t offset = chunk * m_parallel_rank;
  if (chunk + offset > buffer.size())
    chunk = buffer.size() - offset;
  if (offset < buffer.size())
    f.read(&buffer[offset], chunk, index * chunk);
  p += chunk;
}

void Wavefunction::getAll(File &f, int index) {
  auto p = profiler->push("Wavefunction::getAll");
  getw(f, index);
  size_t chunk = (buffer.size() - 1) / m_parallel_size + 1;
  gather_chunks(&buffer[0], buffer.size(), chunk, m_communicator);
  p += buffer.size();
}

void Wavefunction::replicate() {
  size_t chunk = (buffer.size() - 1) / m_parallel_size + 1;
  gather_chunks(&buffer[0], buffer.size(), chunk, m_communicator);
}

std::vector<std::size_t> Wavefunction::histogram(const std::vector<double> &edges, bool parallel, std::size_t start,
                                                 std::size_t stop) {
  if (stop > edges.size())
    stop = edges.size();
  std::vector<std::size_t> cumulative(edges.size(), 0);
  if (parallel)
    DivideTasks(edges.size(), 0, 0, m_communicator);
  for (size_t i = start; i < stop; i++) {
    if (parallel && NextTask(m_communicator))
      for (double j : buffer) {
        if (std::fabs(j) > edges[i])
          cumulative[i]++;
      }
  }
  if (parallel) {
    EndTasks(m_communicator);
    std::vector<double> dcumulative(cumulative.size());
    for (size_t i = 0; i < cumulative.size(); i++)
      dcumulative[i] = (double)cumulative[i];
    gsum(&dcumulative[0], cumulative.size(), m_communicator);
    for (size_t i = 0; i < cumulative.size(); i++)
      cumulative[i] = (std::size_t)dcumulative[i];
  }
  return cumulative;
}

double Wavefunction::norm(const double k) {
  auto result = (double)0;
  size_t imin = 0;
  size_t imax = buffer.size();
  for (size_t i = imin; i < imax; i++) {
    if (buffer[i])
      result += pow(fabs(buffer[i]), k);
  }
  return result;
}

Wavefunction &Wavefunction::addAbsPower(const Wavefunction &c, const double k, const double factor) {

  if (!compatible(c))
    throw std::domain_error("attempt to add incompatible Wavefunction objects");
  //  cout <<"addAbsPower initial=";
  //  for (size_t i=0; i<buffer.size(); i++)
  //    cout <<" "<<buffer[i];
  //  cout <<std::endl;
  size_t imin = 0;
  size_t imax = buffer.size();
  for (size_t i = imin; i < imax; i++) {
    if (k == -1)
      buffer[i] += c.buffer[i] < 0 ? -factor : factor;
    else if (c.buffer[i])
      buffer[i] += factor * pow(fabs(c.buffer[i]), k) * c.buffer[i];
  }
  //  cout <<"addAbsPower final=";
  //  for (size_t i=0; i<buffer.size(); i++)
  //    cout <<" "<<buffer[i];
  //  cout <<std::endl;
  return *this;
}

molpro::vector<double>::iterator Wavefunction::begin() { return buffer.begin(); }

molpro::vector<double>::iterator Wavefunction::end() { return buffer.end(); }

molpro::vector<double>::const_iterator Wavefunction::cbegin() const { return buffer.cbegin(); }

molpro::vector<double>::const_iterator Wavefunction::cend() const { return buffer.cend(); }

std::vector<StringSet> Wavefunction::activeStrings(bool spinUp) const {
  auto p = profiler->push("activeStrings");
  const std::vector<StringSet> &sources = spinUp ? alphaStrings : betaStrings;
  //  return sources;
  const std::vector<StringSet> &complements = spinUp ? betaStrings : alphaStrings;
  std::vector<StringSet> results(8);
  if (m_sparse) {
    auto axis = spinUp ? 1 : 0;
    String proto;
    for (unsigned int sym = 0; sym < 8; sym++) {
      if (!sources[sym].empty())
        proto = sources[sym].front();
    }
    for (unsigned int sym = 0; sym < 8; sym++) {
      results[sym] = StringSet(proto, false, static_cast<int>(sym));
    }
    for (const auto &ww : buffer_sparse) {
      const auto det = determinantAt(ww.first);
      //      cout << "activeStrings " << ww.first << " : " << ww.second << std::endl;
      //      cout <<" symmetry calculated "<<stringSymmetry(ww.first,axis)<<std::endl;
      //      cout <<" det "<<det<<std::endl;
      //      cout <<" strings "<<det.stringAlpha<<", "<<det.stringBeta<<std::endl;
      results[stringSymmetry(ww.first, axis)].push_back((spinUp ? det.stringAlpha : det.stringBeta));
    }
  } else {
    for (unsigned int sym = 0; sym < 8; sym++) {
      const StringSet &source = sources[sym];
      if (source.empty())
        continue;
      const auto syma = spinUp ? sym : symmetry ^ sym;
      const auto symb = symmetry ^ syma;
      results[sym] = StringSet(source.front(), false, static_cast<int>(sym));
      size_t off = _blockOffset[syma];
      size_t dout = spinUp ? betaStrings[symb].size() : 1;
      size_t din = spinUp ? 1 : betaStrings[symb].size();
      size_t nc = complements[sym ^ symmetry].size();
      for (const auto &s : source) {
        for (size_t c = 0; c < nc; c++) {
          if (std::abs(buffer[off + c * din]) > m_activeStringTolerance) {
            results[sym].push_back(s);
            break;
          }
        }
        off += dout;
      }
    }
  }
  //  for (unsigned int sym = 0; sym < 8; sym++) {
  //    cout << "symmetry "<<sym<<std::endl;
  //    for (const auto &r : results[sym])
  //      cout << r << std::endl;
  //  }
  return results;
}
} // namespace gci
} // namespace molpro
