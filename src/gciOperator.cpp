#include "gci.h"
#include "gciOperator.h"
#include <numeric>
#include <algorithm>

gci::Operator gci::Operator::construct(const FCIdump &dump) {
  std::vector<char> portableByteStream;
  int lPortableByteStream;
  int rank = 0;
#ifdef HAVE_MPI_H
  MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
#endif
  if (rank == 0) {
    int verbosity = 0;
    std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
    dim_t dim(8);
    for (const auto &s : orbital_symmetries)
      dim.at(s - 1)++;

    gci::Operator result(dim, 2, dump.parameter("IUHF")[0] > 0, 0, true, true, "Hamiltonian");
//  for (auto i=0; i<orbital_symmetries.size(); i++) xout << "i="<<i+1<<", symmetry="<<orbital_symmetries[i]<<", offset="<<result.offset(i+1)<<std::endl;;

    dump.rewind();
    double value;
    FCIdump::integralType type;
    int i, j, k, l;
    auto &integrals_a = result.O1(true);
    integrals_a.assign(0);
    auto &integrals_b = result.O1(false);
    integrals_b.assign(0);
    auto &integrals_aa = result.O2(true, true);
    auto &integrals_ab = result.O2(true, false);
    auto &integrals_bb = result.O2(false, false);
    if (verbosity > 0) {
      xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
      xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
      xout << "integral addresses " << &integrals_aa << " " << &integrals_ab << " " << &integrals_bb << std::endl;
      xout << "integral sizes " << integrals_aa.size() << " " << integrals_ab.size() << " " << integrals_bb.size()
           << std::endl;
    }
    off_t si,sj,sk,sl,oi,oj,ok,ol;
    while ((type = dump.nextIntegral(si, oi, sj, oj, sk, ok, sl, ol, value)) != FCIdump::endOfFile) {
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
      if (si < sj || (si == sj && oi < oj)) {
        std::swap(oi, oj);
        std::swap(si, sj);
      }
      if (sk < sl || (sk == sl && ok < ol)) {
        std::swap(ok, ol);
        std::swap(sk, sl);
      }
      unsigned int sij = si ^sj;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

      if (type == FCIdump::I2aa) {
        if (verbosity > 2) xout << "aa(" << i << j << "|" << k << l << ") = " << value << std::endl;
        if (verbosity > 2) xout << "aa(" << k << l << "|" << i << j << ") = " << value << std::endl;
        (sij ? integrals_aa.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_aa.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
        (sij ? integrals_aa.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) : integrals_aa.smat(sij, sk, ok, ol)->block(si)[
            oi * (oi + 1) / 2 + oj]) = value;
//          if (sij)
//           xout << "aa("<< i << j <<"|"<< k << l <<") = " << value
//                <<" "<< integrals_aa.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                <<" "<< &integrals_aa.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                 -&integrals_aa.smat(0,0,0,0)->block(0)[0]
//                <<std::endl;
      } else if (type == FCIdump::I2ab) {
        if (verbosity > 2) xout << "ab(" << i << j << "|" << k << l << ") = " << value << std::endl;
        (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_ab.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
//          if (sij)
//           xout << "ab("<< i << j <<"|"<< k << l <<") = " << value
//                <<" "<< integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                <<" "<< &integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                 -&integrals_ab.smat(0,0,0,0)->block(0)[0]
//                <<std::endl;
      } else if (type == FCIdump::I2bb) {
        if (verbosity > 2) xout << "bb(" << i << j << "|" << k << l << ") = " << value << std::endl;
        if (verbosity > 2) xout << "bb(" << k << l << "|" << i << j << ") = " << value << std::endl;
        (sij ? integrals_bb.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_bb.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
        (sij ? integrals_bb.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) : integrals_bb.smat(sij, sk, ok, ol)->block(si)[
            oi * (oi + 1) / 2 + oj]) = value;
      } else if (type == FCIdump::I1a) {
        if (verbosity > 1) xout << "ha(" << i << "," << j << ") = " << value << std::endl;
        integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
      } else if (type == FCIdump::I1b) {
        if (verbosity > 1) xout << "hb(" << i << "," << j << ") = " << value << std::endl;
        integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
      } else if (type == FCIdump::I0)
        result.m_O0 = value;
    }
    if (verbosity > 0) xout << result << std::endl;
    if (verbosity > 1) xout << "int1:\n" << result.int1(1) << std::endl;
    if (verbosity > 1) xout << "intJ:\n" << result.intJ(1, 1) << std::endl;
    if (verbosity > 1) xout << "intK:\n" << result.intK(1) << std::endl;
    portableByteStream = result.bytestream().data();
    lPortableByteStream = portableByteStream.size();
  }
#ifdef HAVE_MPI_H
  MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, MPI_COMM_COMPUTE);
#endif
  char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
  MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, MPI_COMM_COMPUTE);
#endif
  auto result = Operator::construct(buf);
  if (rank != 0) free(buf);
  return result;
}

gci::Operator gci::Operator::construct(const char *dump) {
  class bytestream bs(dump);
  auto so = SymmetryMatrix::Operator::construct(bs);
  return gci::Operator(so);
}

void gci::Operator::FCIDump(const std::string filename, std::vector<int> orbital_symmetries) const {
  FCIdump dump;
  int verbosity = 0;
  if (orbital_symmetries.empty())
    for (auto sym = 0; sym < 8; sym++)
      for (auto i = 0; i < m_dimensions[0][sym]; i++)
        orbital_symmetries.push_back(sym);
  size_t n = orbital_symmetries.size();
  dump.addParameter("IUHF", m_uhf ? 1 : 0);
  dump.addParameter("ORBSYM", orbital_symmetries);
  dump.addParameter("NORB", int(orbital_symmetries.size()));

  dump.write(filename);
  dump.rewind();
  size_t i, j, k, l;
  const auto &integrals_a = O1(true);
  const auto &integrals_b = O1(false);
  const auto &integrals_aa = O2(true, true);
  const auto &integrals_ab = O2(true, false);
  const auto &integrals_bb = O2(false, false);
  if (verbosity > 0) {
    xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
    xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
    xout << "integral addresses " << &integrals_aa << " " << &integrals_ab << " " << &integrals_bb << std::endl;
//      xout << "integral sizes "<< integrals_aa.size()<<" "<<integrals_ab.size()<<" "<<integrals_bb.size()<<std::endl;
    xout << "n=" << n << std::endl;
    xout << "m_rank=" << m_rank << std::endl;
  }
  if (m_uhf) throw std::logic_error("UHF not supported");
  if (m_rank > 1) { // xout << "2 electron"<<std::endl;
    for (i = 1; i <= n; i++)
      for (j = 1; j <= i; j++) { // xout << "j="<<j<<std::endl;
        for (k = 1; k <= i; k++) { // xout << "k="<<k<<std::endl;
          for (l = 1; l <= (i == k ? j : k); l++) {
            auto oi = dump.orbital_offset(i);
            auto oj = dump.orbital_offset(j);
            auto ok = dump.orbital_offset(k);
            auto ol = dump.orbital_offset(l);
            auto si = dump.orbital_symmetry(i);
            auto sj = dump.orbital_symmetry(j);
            auto sk = dump.orbital_symmetry(k);
            auto sl = dump.orbital_symmetry(l);
//      xout << "ijkl "<<i<<j<<k<<l<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
            if (si < sj || (si == sj && oi < oj)) {
              std::swap(oi, oj);
              std::swap(si, sj);
            }
            if (sk < sl || (sk == sl && ok < ol)) {
              std::swap(ok, ol);
              std::swap(sk, sl);
            }
            unsigned int sij = si ^sj;
            unsigned int skl = sk ^sl;
            if ((sij ^ skl) != m_symmetry) continue;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
//      xout << (sij ? integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) <<" "<< i << " "<<j<<" "<<k<<" "<<l <<std::endl;;
//      xout << (sij ? &integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : &integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) <<" "<< i << " "<<j<<" "<<k<<" "<<l <<std::endl;;
            auto value = (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_ab.smat(sij,
                                                                                                             si,
                                                                                                             oi,
                                                                                                             oj)->block(
                sk)[ok * (ok + 1) / 2 + ol]);
            if (std::abs(value) > 1e-15)
              dump.writeIntegral(i, j, k, l, value);
          }
        }
      }
  }
//  xout << "n="<<n<<std::endl;
  if (m_rank > 0)
    for (size_t ii = 1; ii <= n; ii++)
      for (size_t jj = 1; jj <= ii; jj++) {
//        xout << "ii="<<ii<<", jj="<<jj<<std::endl;
        auto oi = dump.orbital_offset(ii);
        auto oj = dump.orbital_offset(jj);
        auto si = dump.orbital_symmetry(ii);
        auto sj = dump.orbital_symmetry(jj);
        if ((si ^ sj) == m_symmetry) {
          auto value = integrals_a.block(si).at(oi * (oi + 1) / 2 + oj);
          if (std::abs(value) > 1e-15)
            dump.writeIntegral(ii, jj, 0, 0, value);
        }
      }
  dump.writeIntegral(0, 0, 0, 0, m_O0);

}

gci::Operator *gci::Operator::projector(const std::string special, const bool forceSpinUnrestricted) const {
  auto result = new gci::Operator(m_dimensions[0],
                                  1,
                                  m_uhf > 0 || forceSpinUnrestricted,
                                  0,
                                  true,
                                  true,
                                  special + " projector");
  result->m_O0 = 0;

  if (special == "P" || special == "Q") {
    if (special == "P")
      result->O1(true).setIdentity();
    else
      result->O1(true).assign(0);
    if (m_uhf)
      result->O1(false) = result->O1(true);
    // determine the uncoupled orbital
    size_t uncoupled_orbital = 0;
    unsigned int uncoupled_orbital_symmetry = 0;
    double min_rowsum = 1e50;
    for (unsigned int sym = 0; sym < 8; sym++) {
      for (size_t k = 0; k < m_dimensions[0][sym]; k++) {
        double rowsum = 0;
        for (size_t l = 0; l < m_dimensions[0][sym]; l++)
          rowsum += O1(true).block(sym)[k > l ? k * (k + 1) / 2 + l : l * (l + 1) / 2 + k];
        if (std::fabs(rowsum) < min_rowsum) {
          min_rowsum = std::fabs(rowsum);
          uncoupled_orbital = k;
          uncoupled_orbital_symmetry = sym;
        }
      }
    }
    xout << "non-interacting orbital is " << uncoupled_orbital + 1 << "." << uncoupled_orbital_symmetry + 1
         << std::endl;
    result->O1(false).block(uncoupled_orbital_symmetry)[(uncoupled_orbital + 2) * (uncoupled_orbital + 1) / 2 - 1] =
        (special == "P" ? 0 : 1);
  }
  return result;
}

Eigen::VectorXd gci::Operator::int1(int spin) const {
  size_t basisSize = 0;
  for (auto si=0; si<8; si++)
    basisSize += O1(spin>0).dimension(si);
  Eigen::VectorXd result(basisSize);
  size_t off=0;
  for (auto si=0; si<8; si++)
    for (auto oi=0; oi<m_dimensions[0][si]; oi++)
      result[off++] = O1(spin>0).block(si)[(oi+1)*(oi+2)/2-1];
  return result;
}

Eigen::MatrixXd gci::Operator::intJ(int spini, int spinj) const {
  if (spinj > spini) return intJ(spinj, spini).transpose();
  size_t basisSize = 0;
  for (auto si=0; si<8; si++)
    basisSize += O1().dimension(si);
  Eigen::MatrixXd result(basisSize, basisSize);
  size_t i = 0;
  for (auto si = 0; si < 8; si++)
    for (auto oi = 0; oi < dimension(si, 0, spini > 0); oi++) {
      size_t j = 0;
      for (auto sj = 0; sj < 8; sj++)
        for (auto oj = 0; oj < dimension(sj, 0, spinj > 0); oj++)
          result(j++, i) = O2(spini > 0, spinj > 0).smat(0, si, oi, oi)->block(sj)[(oj + 2) * (oj + 1) / 2 - 1];
      i++;
    }
  return result;
}

Eigen::MatrixXd gci::Operator::intK(int spin) const {
  size_t basisSize = 0;
  for (auto si=0; si<8; si++)
    basisSize += O1().dimension(si);
  Eigen::MatrixXd result(basisSize, basisSize);
  size_t i = 0;
  for (auto si = 0; si < 8; si++)
    for (auto oi = 0; oi < dimension(si, 0, spin > 0); oi++) {
      size_t j = 0;
      for (auto sj = 0; sj < 8; sj++)
        for (auto oj = 0; oj < dimension(sj, 0, spin > 0); oj++)
          result(j++, i) =
              ((si < sj) ?
               O2(spin > 0, spin > 0).smat(si ^ sj, si, oi, oj)->blockMap(si)(oi, oj)
                         :
               ((si > sj) ?
                O2(spin > 0, spin > 0).smat(si ^ sj, sj, oj, oi)->blockMap(sj)(oj, oi)
                          : ((i > j) ?
                             O2(spin > 0, spin > 0).smat(0, si, oi, oj)->block(si)[(oi * (oi + 1)) / 2 + oj]
                                     :
                             O2(spin > 0, spin > 0).smat(0, sj, oj, oi)->block(sj)[(oj * (oj + 1)) / 2 + oi]
                )));
      i++;
    }
  return result;
}


SymmetryMatrix::Operator gci::fockOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference, const std::string description) {
  Operator f(hamiltonian.m_dimensions[0],
             1,
             hamiltonian.m_uhf,
             hamiltonian.m_symmetry,
             true,
             true,
             description);
  // xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout << "gci::Operator::fockOperator Reference beta: "<<reference.stringBeta<<std::endl;
  //  bool closed = reference.stringAlpha==reference.stringBeta;
  //  xout << "gci::Operator::fockOperator Reference alpha=beta: "<<closed<<std::endl;
  auto refAlphaOrbitals = reference.stringAlpha.orbitals();
  auto refBetaOrbitals = reference.stringBeta.orbitals();
  f.m_O0 = hamiltonian.m_O0;
  f.O1(true) = hamiltonian.O1(true);
  if (hamiltonian.m_uhf) f.O1(false) = hamiltonian.O1(false);
  size_t basisSize = 0;
  for (auto si=0; si<8; si++)
    basisSize += hamiltonian.O1().dimension(si);
  // xout <<"reference.stringAlpha.orbitals ";for (size_t i=0; i < reference.stringAlpha.orbitals().size(); i++) xout <<reference.stringAlpha.orbitals()[i]<<" ";xout <<std::endl;
  for (const auto &o : refAlphaOrbitals) {
//       xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
//       xout<< "f alpha, alpha occ: " <<*o << std::endl;
    unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
    unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
//    unsigned int oo = offset(o);
    for (unsigned int i = 1; i <= basisSize; i++)
      for (unsigned int j = 1; j <= i; j++) {
        if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
        f.element(reference.orbitalSpace->orbitalIndex(i),
                  reference.orbitalSpace->orbital_symmetries[i - 1],
                  reference.orbitalSpace->orbitalIndex(j),
                  reference.orbitalSpace->orbital_symmetries[j - 1],
                  true) +=
            hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    oo,
                    os,
                    oo,
                    os,
                    true,
                    true) -
                hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                        reference.orbitalSpace->orbital_symmetries[i - 1],
                        oo,
                        os,
                        oo,
                        os,
                        reference.orbitalSpace->orbitalIndex(j),
                        reference.orbitalSpace->orbital_symmetries[j - 1],
                        true,
                        true);
      }
  }
  for (const auto &o : refBetaOrbitals) {
    // xout<< "f alpha, beta occ: " <<*o << std::endl;
    unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
    unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
    for (unsigned int i = 1; i <= basisSize; i++)
      for (unsigned int j = 1; j <= i; j++) {
        if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
        f.element(reference.orbitalSpace->orbitalIndex(i),
                  reference.orbitalSpace->orbital_symmetries[i - 1],
                  reference.orbitalSpace->orbitalIndex(j),
                  reference.orbitalSpace->orbital_symmetries[j - 1],
                  true) +=
            hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    oo,
                    os,
                    oo,
                    os,
                    true,
                    false);
      }
  }
  if (f.m_uhf) {
    for (const auto &o : refBetaOrbitals) {
      // xout<< "f beta, beta occ: " <<*o << std::endl;
      unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
      unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
      for (unsigned int i = 1; i <= basisSize; i++)
        for (unsigned int j = 1; j <= i; j++) {
          if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
          f.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    false) +=
              hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                      reference.orbitalSpace->orbital_symmetries[i - 1],
                      reference.orbitalSpace->orbitalIndex(j),
                      reference.orbitalSpace->orbital_symmetries[j - 1],
                      oo,
                      os,
                      oo,
                      os,
                      false,
                      false) -
                  hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                          reference.orbitalSpace->orbital_symmetries[i - 1],
                          oo,
                          os,
                          oo,
                          os,
                          reference.orbitalSpace->orbitalIndex(j),
                          reference.orbitalSpace->orbital_symmetries[j - 1],
                          false,
                          false);
        }
    }
    for (const auto &o : refAlphaOrbitals) {
      // xout<< "f beta, alpha occ: " <<*o << std::endl;
      unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
      unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
      for (unsigned int i = 1; i <= basisSize; i++)
        for (unsigned int j = 1; j <= i; j++) {
          if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
          f.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    false) +=
              hamiltonian.element(oo,
                      os,
                      oo,
                      os,
                      reference.orbitalSpace->orbitalIndex(i),
                      reference.orbitalSpace->orbital_symmetries[i - 1],
                      reference.orbitalSpace->orbitalIndex(j),
                      reference.orbitalSpace->orbital_symmetries[j - 1],
                      true,
                      false);
        }
    }
  }
  return f;
}

SymmetryMatrix::Operator gci::sameSpinOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference, const std::string description) {
  Operator result(hamiltonian.m_dimensions[0],
                  hamiltonian.m_rank,
                  true,
                  hamiltonian.m_symmetry,
                  true,
                  true,
                  description);
  {
    Determinant ra = reference;
    ra.stringBeta.nullify();
    auto f = fockOperator(hamiltonian,ra);
    for (size_t i = 0; i < hamiltonian.O1(true).size(); i++)
      (*result.O1(true).data())[i] = (*hamiltonian.O1(true).data())[i] - (*f.O1(true).data())[i];
  }
  {
    Determinant ra = reference;
    ra.stringAlpha.nullify();
    auto f = fockOperator(hamiltonian,ra);
    for (size_t i = 0; i < hamiltonian.O1(false).size(); i++)
      (*result.O1(false).data())[i] = (*hamiltonian.O1(false).data())[i] - (*f.O1(false).data())[i];
  }

  *result.O2(true, true).data() = *hamiltonian.O2(true, true).data();
  *result.O2(false, false).data() = *hamiltonian.O2(false, false).data();
  for (auto &s : *result.O2(true, false).data()) s = 0;
  return result;
}

void gci::gsum(SymmetryMatrix::Operator& op) {
  std::vector<double> O0;
  O0.push_back(op.m_O0);
  ::gci::gsum(O0);
  op.m_O0 = O0[0];
  if (op.m_rank > 0) {
    ::gci::gsum(*op.O1(true).data());
    if (op.m_uhf)
      ::gci::gsum(*op.O1(false).data());
  }
  if (op.m_rank > 1) {
    ::gci::gsum(*op.O2(true, true).data());
    if (op.m_uhf) {
      ::gci::gsum(*op.O2(true, false).data());
      ::gci::gsum(*op.O2(false, false).data());
    }
  }
}
