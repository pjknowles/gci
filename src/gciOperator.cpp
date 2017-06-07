#include "gci.h"
#include "gciOperator.h"

gci::Operator gci::Operator::construct(const FCIdump &dump)
{
  int verbosity=0;
  std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
  dim_t dim(8);
  for (auto& s : orbital_symmetries)
    dim.at(--s)++;

  Operator result(dims_t(4,dim),2,dump.parameter("IUHF")[0]>0,{1,1},{-1,-1},0,"Hamiltonian");
  result.m_orbital_symmetries = orbital_symmetries;
  result.m_fcidump = &dump;
//  for (auto i=0; i<orbital_symmetries.size(); i++)
//      xout << "i="<<i+1<<", symmetry="<<orbital_symmetries[i]<<", offset="<<result.offset(i+1)<<std::endl;;

  dump.rewind();
  double value;
  FCIdump::integralType type;
  int i,j,k,l;
  SMat& integrals_a = result.O1(true); integrals_a.assign(0);
  SMat& integrals_b = result.O1(false);integrals_b.assign(0);
  SMatMat& integrals_aa = result.O2(true,true);
  SMatMat& integrals_ab = result.O2(true,false);
  SMatMat& integrals_bb = result.O2(false,false);
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
      xout << "\nvalue: "<<value<<std::endl;
      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

      if (type == FCIdump::I2aa) {
          if (verbosity>2) xout << "aa("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          if (verbosity>2) xout << "aa("<< k << l <<"|"<< i << j <<") = " << value <<std::endl;
          (sij ? integrals_aa.smat(sij,si,oi,oj)->blockM(sk)(ok,ol) : integrals_aa.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) = value;
          (sij ? integrals_aa.smat(sij,sk,ok,ol)->blockM(si)(oi,oj) : integrals_aa.smat(sij,sk,ok,ol)->block(si)[oi*(oi+1)/2+oj]) = value;
          if (sij)
           xout << "aa("<< i << j <<"|"<< k << l <<") = " << value
                <<" "<< integrals_aa.smat(sij,si,oi,oj)->blockM(sk)(ok,ol)
                <<" "<< &integrals_aa.smat(sij,si,oi,oj)->blockM(sk)(ok,ol)
                 -&integrals_aa.smat(0,0,0,0)->block(0)[0]
                <<std::endl;
        } else if (type == FCIdump::I2ab) {
          if (verbosity>2) xout << "ab("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          sij ? integrals_ab.smat(sij,si,oi,oj)->blockM(sk)(ok,ol) : integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol] = value;
        } else if (type == FCIdump::I2bb) {
          if (verbosity>2) xout << "bb("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          if (verbosity>2) xout << "bb("<< k << l <<"|"<< i << j <<") = " << value <<std::endl;
          sij ? integrals_bb.smat(sij,si,oi,oj)->blockM(sk)(ok,ol) : integrals_bb.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol] = value;
          sij ? integrals_bb.smat(sij,sk,ok,ol)->blockM(si)(oi,oj) : integrals_bb.smat(sij,sk,ok,ol)->block(si)[oi*(oi+1)/2+oj] = value;
        } else if (type == FCIdump::I1a) {
          if (verbosity>1) xout << "ha("<< i <<","<< j <<") = " << value <<std::endl;
          integrals_a.block(si).at(oi*(oi+1)/2+oj)=value;
        } else if (type == FCIdump::I1b) {
          if (verbosity>1) xout << "hb("<< i <<","<< j <<") = " << value <<std::endl;
          integrals_b.block(si).at(oi*(oi+1)/2+oj)=value;
        } else if (type == FCIdump::I0)
        result.m_O0 = value;
    }
  if (verbosity>0) xout << result << std::endl;
  if (verbosity>-1) xout << "int1:\n" <<result.int1(1)<< std::endl;
  if (verbosity>-1) xout << "intJ:\n" <<result.intJ(1,1)<< std::endl;
  if (verbosity>-1) xout << "intK:\n" <<result.intK(1)<< std::endl;

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
             O2(spin>0,spin>0).smat(si^sj,si,oi,oj)->blockM(si)(oi,oj)
           :
             O2(spin>0,spin>0).smat(si^sj,sj,oj,oi)->blockM(sj)(oj,oi)
             );
    }
  }
  return result;
}
