#include "gciOrbitalSpace.h"

void OrbitalSpace::load(const Options &dump, int verbosity) {
  if (verbosity) xout <<"Load OrbitalSpace from Options object" <<std::endl;
  load(dump.parameter("ORBSYM",std::vector<int>{0}),
  dump.parameter("IUHF")!=0,
       verbosity);
}

void OrbitalSpace::load(const std::vector<int> syms, bool uhf, int verbosity) {
  spinUnrestricted=uhf;
  orbital_symmetries = std::vector<unsigned int>(syms.size(),0);
  for (auto s=syms.begin(); s!=syms.end(); s++) {
    orbital_symmetries[s-syms.begin()]=(*s)-1; // convert 1-8 to 0-7
    at((*s)-1)++;
  }
  calculateOffsets(); // set up the rest of the underlying SymmetrySpace object
  //    xout << this->str(2) << std::endl;

  pairSpace[-1] = SymmetrySpace("Antisymmetric pairs of orbitals");
  pairSpace[0] = SymmetrySpace("Pairs of orbitals");
  pairSpace[1] = SymmetrySpace("Symmetric pairs of orbitals");
  for (int isym=0; isym<8; isym++) {
    pairSpace[-1][isym] = total(isym,-1);
    pairSpace[0][isym] = total(isym,0);
    pairSpace[1][isym] = total(isym,1);
  }
  pairSpace[-1].calculateOffsets();
  pairSpace[0].calculateOffsets();
  pairSpace[1].calculateOffsets();
  //    xout << "symmetricPairSpace constructed: " << pairSpace[1];


}

size_t OrbitalSpace::pairIndex(unsigned int i, unsigned int j, int parity) const {
  unsigned int ii = orbitalIndex(i);
  unsigned int jj = orbitalIndex(j);
  unsigned int isym = orbital_symmetries[i-1];
  unsigned int jsym = orbital_symmetries[j-1];
  if (parity == 0 || isym > jsym ) {
//      xout < i<<" "<<j<<" "<< offset(isym^jsym,isym,parity) + jj*at(isym)+ii<<std::endl;
    return offset(isym^jsym,isym,parity) + jj*at(isym)+ii;
  } else {
    if (parity > 0 && isym == jsym)
      return offset(0,isym,parity) + ((ii>jj) ? (ii*(ii+1))/2+jj : (jj*(jj+1))/2+ii);
    else if (parity < 0 && isym == jsym)
      return offset(0,isym,parity) + ((ii>jj) ? (ii*(ii-1))/2+jj : (jj*(jj-1))/2+ii);
    else
      return offset(isym^jsym,jsym,parity) + ii*at(jsym) + jj;
  }
}

size_t OrbitalSpace::quadIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l, int parity, int parity2) const {
  unsigned int isym = orbital_symmetries[i-1];
  unsigned int jsym = orbital_symmetries[j-1];
  unsigned int ksym = orbital_symmetries[k-1];
  unsigned int lsym = orbital_symmetries[l-1];
  return pairSpace.find(parity)->second.offset(isym^jsym^ksym^lsym,isym^jsym,parity2)
      + pairIndex(i,j,parity)
      + pairIndex(k,l,parity) * total(isym^jsym,parity);
}

unsigned int OrbitalSpace::orbitalIndex(unsigned int i) const {
  unsigned int n=0;
  for (unsigned int j=1; j<i; j++)
    if (orbital_symmetries[j-1] == orbital_symmetries[i-1]) n++;
  return n;
}

std::string OrbitalSpace::str(int verbosity, unsigned int columns) const
{
  std::ostringstream o;
  o << SymmetrySpace::str(verbosity);
  if (verbosity>=0) {
    o << std::endl << "Basis size="<<total()<<" Spin unrestricted? "<<spinUnrestricted;
    for (int i=-1; i<2; i++)
      if (pairSpace.find(i) != pairSpace.end())
        o << std::endl << pairSpace.find(i)->second.str(verbosity);
  }
  if (verbosity>=1) {
    o << std::endl << "Orbital symmetries";
    for (std::vector<unsigned int>::const_iterator i=orbital_symmetries.begin(); i!=orbital_symmetries.end(); i++)
      o<<" "<<*i+1;
  }
  return o.str();
}
