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
  for (auto i=0; i<orbital_symmetries.size(); i++) {
      xout << "i="<<i+1<<", symmetry="<<orbital_symmetries[i]<<", offset="<<result.offset(i+1)<<std::endl;;
    }

  dump.rewind();
  double value;
  FCIdump::integralType type;
  int i,j,k,l;
  SMat& integrals_a = result.O1(true);
  SMat& integrals_b = result.O1(false);
  SMatMat& integrals_aa = result.O2(true,true);
  SMatMat& integrals_ab = result.O2(true,false);
  SMatMat& integrals_bb = result.O2(false,false);
  while ((type=dump.nextIntegral(i,j,k,l,value))!=FCIdump::endOfFile) {
      auto oi = result.offset(i);
      auto oj = result.offset(j);
      auto ok = result.offset(k);
      auto ol = result.offset(l);
      auto si = i<1 ? 0 : orbital_symmetries[i-1];
      auto sj = j<1 ? 0 : orbital_symmetries[j-1];
      auto sk = k<1 ? 0 : orbital_symmetries[k-1];
      auto sl = l<1 ? 0 : orbital_symmetries[l-1];
      xout << "ijkl "<<i<<j<<k<<l<<std::endl;
      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
      if (si>sj || (si==sj && oi < oj)) { std::swap(oi,oj); std::swap(si,sj);}
      if (sk>sl || (sk==sl && ok < ol)) { std::swap(ok,ol); std::swap(sk,sl);}
      unsigned int sij = si^sj;
      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

      if (type == FCIdump::I2aa) {
          if (verbosity>2) xout << "aa("<< i << j <<"|"<< k << l <<") = " << value <<std::endl;
          if (verbosity>2) xout << "aa("<< k << l <<"|"<< i << j <<") = " << value <<std::endl;
          if (sij) xout << "block "<<integrals_aa.smat(sij,si,oi,oj)->blockM(sk)<<std::endl;
          sij ? integrals_aa.smat(sij,si,oi,oj)->blockM(sk)(ok,ol) : integrals_aa.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol] = value;
          if (sij) xout << "block "<<integrals_aa.smat(sij,si,ok,ol)->blockM(sk)<<std::endl;
          sij ? integrals_aa.smat(sij,sk,ok,ol)->blockM(si)(oi,oj) : integrals_aa.smat(sij,sk,ok,ol)->block(si)[oi*(oi+1)/2+oj] = value;
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

  return result;
}
