#include "gci.h"
#include "gciOperator.h"

gci::Operator gci::Operator::construct(const FCIdump &dump)
{
  int verbosity=0;
  std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
  dim_t dim(8);
  for (auto& s : orbital_symmetries)
    dim.at(s-1)++;

  bool uhf = dump.parameter("IUHF")[0]>0;
  gci::Operator result(dims_t(4,dim),2,uhf,{1,1},{-1,-1},0,"Hamiltonian");
  for (auto i=0; i<4; i++) result.m_orbitalSpaces.push_back(OrbitalSpace(orbital_symmetries,uhf)); // in this implementation, all four orbital spaces are the same
  for (auto& s : orbital_symmetries ) result.m_orbital_symmetries.push_back(--s);
  result.m_fcidump = &dump;
//  for (auto i=0; i<orbital_symmetries.size(); i++) xout << "i="<<i+1<<", symmetry="<<orbital_symmetries[i]<<", offset="<<result.offset(i+1)<<std::endl;;
//  for (auto i=0; i<orbital_symmetries.size(); i++) xout << "i="<<i+1<<", symmetry="<<result.m_orbital_symmetries[i]<<std::endl;;

  dump.rewind();
  double value;
  FCIdump::integralType type;
  int i,j,k,l;
  auto& integrals_a = result.O1(true); integrals_a.assign(0);
  auto& integrals_b = result.O1(false);integrals_b.assign(0);
  auto& integrals_aa = result.O2(true,true);
  auto& integrals_ab = result.O2(true,false);
  auto& integrals_bb = result.O2(false,false);
  if (verbosity>0) {
      xout << "integral addresses "<< &integrals_a<<" "<<&integrals_b<<std::endl;
      xout << "integral addresses "<< &integrals_a.block(0)[0]<<" "<<&integrals_b.block(0)[0]<<std::endl;
      xout << "integral addresses "<< &integrals_aa<<" "<<&integrals_ab<<" "<<&integrals_bb<<std::endl;
      xout << "integral sizes "<< integrals_aa.size()<<" "<<integrals_ab.size()<<" "<<integrals_bb.size()<<std::endl;
    }
  while ((type=dump.nextIntegral(i,j,k,l,value))!=FCIdump::endOfFile) {
      auto oi = result.offset(i);
      auto oj = result.offset(j);
      auto ok = result.offset(k);
      auto ol = result.offset(l);
      auto si = i<1 ? 0 : orbital_symmetries[i-1];
      auto sj = j<1 ? 0 : orbital_symmetries[j-1];
      auto sk = k<1 ? 0 : orbital_symmetries[k-1];
      auto sl = l<1 ? 0 : orbital_symmetries[l-1];
//      xout << "ijkl "<<i<<j<<k<<l<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
      if (si<sj || (si==sj && oi < oj)) { std::swap(oi,oj); std::swap(si,sj);}
      if (sk<sl || (sk==sl && ok < ol)) { std::swap(ok,ol); std::swap(sk,sl);}
      unsigned int sij = si^sj;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

      if (type == FCIdump::I2aa) {
          if (verbosity>2) xout << "aa("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          if (verbosity>2) xout << "aa("<< k << l <<"|"<< i << j <<") = " << value <<std::endl;
          (sij ? integrals_aa.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : integrals_aa.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) = value;
          (sij ? integrals_aa.smat(sij,sk,ok,ol)->blockMap(si)(oi,oj) : integrals_aa.smat(sij,sk,ok,ol)->block(si)[oi*(oi+1)/2+oj]) = value;
//          if (sij)
//           xout << "aa("<< i << j <<"|"<< k << l <<") = " << value
//                <<" "<< integrals_aa.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                <<" "<< &integrals_aa.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                 -&integrals_aa.smat(0,0,0,0)->block(0)[0]
//                <<std::endl;
        } else if (type == FCIdump::I2ab) {
          if (verbosity>2) xout << "ab("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          (sij ? integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) = value;
//          if (sij)
//           xout << "ab("<< i << j <<"|"<< k << l <<") = " << value
//                <<" "<< integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                <<" "<< &integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol)
//                 -&integrals_ab.smat(0,0,0,0)->block(0)[0]
//                <<std::endl;
        } else if (type == FCIdump::I2bb) {
          if (verbosity>2) xout << "bb("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          if (verbosity>2) xout << "bb("<< k << l <<"|"<< i << j <<") = " << value <<std::endl;
          (sij ? integrals_bb.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : integrals_bb.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) = value;
          (sij ? integrals_bb.smat(sij,sk,ok,ol)->blockMap(si)(oi,oj) : integrals_bb.smat(sij,sk,ok,ol)->block(si)[oi*(oi+1)/2+oj]) = value;
        } else if (type == FCIdump::I1a) {
          if (verbosity>1) xout << "ha("<< i <<","<< j <<") = " << value <<std::endl;
          integrals_a.block(si).at(oi*(oi+1)/2+oj)=value;
        } else if (type == FCIdump::I1b) {
          if (verbosity>1) xout << "hb("<< i <<","<< j <<") = " << value <<std::endl;
          integrals_b.block(si).at(oi*(oi+1)/2+oj)=value;
        } else if (type == FCIdump::I0)
        result.m_O0 = value;
    }
//  exit(0);
  if (verbosity>0) xout << result << std::endl;
  if (verbosity>1) xout << "int1:\n" <<result.int1(1)<< std::endl;
  if (verbosity>1) xout << "intJ:\n" <<result.intJ(1,1)<< std::endl;
  if (verbosity>1) xout << "intK:\n" <<result.intK(1)<< std::endl;

  return result;
}

gci::Operator* gci::Operator::projector(const std::string special, const bool forceSpinUnrestricted) const
{
  auto result = new gci::Operator({m_dimensions[0]},1,m_uhf>0||forceSpinUnrestricted,m_hermiticity,m_exchange,0,special+" projector");
  result->m_orbital_symmetries = m_orbital_symmetries;
  result->m_O0=0;

  if (special=="P" || special=="Q") {
      if (special=="P")
        result->O1(true).setIdentity();
      else
        result->O1(true).assign(0);
      if (m_uhf)
        result->O1(false) = result->O1(true);
      // determine the uncoupled orbital
      size_t uncoupled_orbital=0;
      unsigned int uncoupled_orbital_symmetry=0;
      double min_rowsum=1e50;
      for (unsigned int sym=0; sym<8; sym++) {
          for (size_t k=0; k<m_dimensions[0][sym]; k++) {
              double rowsum=0;
              for (size_t l=0; l<m_dimensions[0][sym]; l++)
                rowsum+=O1(true).block(sym)[k>l ? k*(k+1)/2+l : l*(l+1)/2+k];
              if (std::fabs(rowsum) < min_rowsum) {
                  min_rowsum=std::fabs(rowsum);
                  uncoupled_orbital=k;
                  uncoupled_orbital_symmetry=sym;
                }
            }
        }
      xout << "non-interacting orbital is "<<uncoupled_orbital<<"."<<uncoupled_orbital_symmetry<<std::endl;
      result->O1(false).block(uncoupled_orbital_symmetry)[(uncoupled_orbital+2)*(uncoupled_orbital+1)/2-1]=(special=="P" ? 1:0);
    }
  return result;
}

Eigen::VectorXd gci::Operator::int1(int spin) const
{
  size_t basisSize = m_orbital_symmetries.size();
  Eigen::VectorXd result(basisSize);
  for (unsigned int i=0; i < basisSize; i++) {
    result[i] = O1(spin>0).block(m_orbital_symmetries[i])[(offset(i+1)+1)*(offset(i+1)+2)/2-1];
  }
  return result;
}


Eigen::MatrixXd gci::Operator::intJ(int spini, int spinj) const
{
  if (spinj > spini) return intJ(spinj,spini).transpose();
  size_t basisSize = m_orbital_symmetries.size();
  Eigen::MatrixXd result(basisSize,basisSize);
  for (unsigned int j=0; j < basisSize; j++) {
    for (unsigned int i=0; i < basisSize; i++) {
      auto si = m_orbital_symmetries[i];
      auto sj = m_orbital_symmetries[j];
      auto oi = offset(i+1);
      auto oj = offset(j+1);
      result(i,j) = O2(spini>0,spinj>0).smat(0,si,oi,oi)->block(sj)[(oj+2)*(oj+1)/2-1];
    }
  }
  return result;
}

Eigen::MatrixXd gci::Operator::intK(int spin) const
{
  size_t basisSize = m_orbital_symmetries.size();
  Eigen::MatrixXd result(basisSize,basisSize);
  for (unsigned int j=0; j < basisSize; j++) {
    for (unsigned int i=0; i < basisSize; i++) {
      auto si = m_orbital_symmetries[i];
      auto sj = m_orbital_symmetries[j];
      auto oi = offset(i+1);
      auto oj = offset(j+1);
      result(i,j) =
          ((si < sj) ?
             O2(spin>0,spin>0).smat(si^sj,si,oi,oj)->blockMap(si)(oi,oj)
           :
             ((si > sj) ?
                O2(spin>0,spin>0).smat(si^sj,sj,oj,oi)->blockMap(sj)(oj,oi)
              : ((i>j) ?
                   O2(spin>0,spin>0).smat(0,si,oi,oj)->block(si)[(oi*(oi+1))/2+oj]
                 :
                 O2(spin>0,spin>0).smat(0,sj,oj,oi)->block(sj)[(oj*(oj+1))/2+oi]
             )));
      }
  }
  return result;
}

gci::Operator gci::Operator::fockOperator(const Determinant &reference, const std::string description) const
{
  Operator f({m_dimensions[0]},
             1,
             m_uhf,
             m_hermiticity,
             m_exchange,
             m_symmetry,
             description);
  // xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout << "gci::Operator::fockOperator Reference beta: "<<reference.stringBeta<<std::endl;
  //  bool closed = reference.stringAlpha==reference.stringBeta;
  //  xout << "gci::Operator::fockOperator Reference alpha=beta: "<<closed<<std::endl;
  auto refAlphaOrbitals=reference.stringAlpha.orbitals();
  auto refBetaOrbitals=reference.stringBeta.orbitals();
  f.m_O0 = m_O0;
  f.O1(true) = O1(true);
  if (m_uhf) f.O1(false) = O1(false);
  f.m_orbital_symmetries = m_orbital_symmetries;
  unsigned int basisSize=m_orbital_symmetries.size();
  // xout <<"reference.stringAlpha.orbitals ";for (size_t i=0; i < reference.stringAlpha.orbitals().size(); i++) xout <<reference.stringAlpha.orbitals()[i]<<" ";xout <<std::endl;
  // for (std::vector<unsigned int>::const_iterator o=reference.stringAlpha.orbitals().begin(); o != reference.stringAlpha.orbitals().end(); o++)
  for (auto o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
    {
      // xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
      // xout<< "f alpha, alpha occ: " <<*o << std::endl;
      unsigned int os=m_orbital_symmetries[*o];
      unsigned int oo=offset(*o);
      for (unsigned int i=1; i<=basisSize; i++)
        for (unsigned int j=1; j<=i; j++) {
            if (m_orbital_symmetries[i-1]!=m_orbital_symmetries[j-1]) continue;
            f.element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],true) +=
                element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],oo,os,oo,os,true,true)-
                element(offset(i),m_orbital_symmetries[i],oo,os,oo,os,offset(j),m_orbital_symmetries[j],true,true);
          }
    }
  for (auto o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
    {
      // xout<< "f alpha, beta occ: " <<*o << std::endl;
      unsigned int os=m_orbital_symmetries[*o];
      unsigned int oo=offset(*o);
      for (unsigned int i=1; i<=basisSize; i++)
        for (unsigned int j=1; j<=i; j++) {
            if (m_orbital_symmetries[i-1]!=m_orbital_symmetries[j-1]) continue;
            f.element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],true) +=
                element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],oo,os,oo,os,true,false);
          }
    }
  if (f.m_uhf) {
      for (auto o=refBetaOrbitals.begin(); o != refBetaOrbitals.end(); o++)
        {
          // xout<< "f beta, beta occ: " <<*o << std::endl;
          unsigned int os=m_orbital_symmetries[*o];
          unsigned int oo=offset(*o);
          for (unsigned int i=1; i<=basisSize; i++)
            for (unsigned int j=1; j<=i; j++) {
                if (m_orbital_symmetries[i-1]!=m_orbital_symmetries[j-1]) continue;
                f.element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],false) +=
                    element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],oo,os,oo,os,false,false)-
                    element(offset(i),m_orbital_symmetries[i],oo,os,oo,os,offset(j),m_orbital_symmetries[j],false,false);
              }
        }
      for (auto o=refAlphaOrbitals.begin(); o != refAlphaOrbitals.end(); o++)
        {
          // xout<< "f beta, alpha occ: " <<*o << std::endl;
          unsigned int os=m_orbital_symmetries[*o];
          unsigned int oo=offset(*o);
          for (unsigned int i=1; i<=basisSize; i++)
            for (unsigned int j=1; j<=i; j++) {
                if (m_orbital_symmetries[i-1]!=m_orbital_symmetries[j-1]) continue;
                f.element(offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],false) +=
                    element(oo,os,oo,os,offset(i),m_orbital_symmetries[i],offset(j),m_orbital_symmetries[j],true,false);
              }
        }
    }
  return f;
}


gci::Operator gci::Operator::sameSpinOperator(const Determinant &reference, const std::string description) const
{
  Operator result(m_dimensions,
             m_rank,
             true,
             m_hermiticity,
             m_exchange,
             m_symmetry,
             description);
  {
  Determinant ra = reference; ra.stringBeta.nullify();
  auto f = this->fockOperator(ra);
  for (size_t i=0; i<O1(true).size(); i++)
      (*result.O1(true).data())[i] = (*O1(true).data())[i] - (*f.O1(true).data())[i];
  }
  {
  Determinant ra = reference; ra.stringAlpha.nullify();
  auto f = this->fockOperator(ra);
  for (size_t i=0; i<O1(false).size(); i++)
      (*result.O1(false).data())[i] = (*O1(false).data())[i] - (*f.O1(false).data())[i];
  }

  *result.O2(true,true).data() = *O2(true,true).data();
  *result.O2(false,false).data() = *O2(false,false).data();
  for (auto& s : *result.O2(true,false).data()) s=0;
  result.m_dirac_out_of_date=true;
  return result;
}
