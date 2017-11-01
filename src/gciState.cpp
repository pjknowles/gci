#include "gciState.h"
#include <fstream>
#include <iostream>
#include <sstream>

State::State(const Options& dump)
{
  load(dump);
}

State::State(OrbitalSpace *h, int n, int s, int m2)
  : orbitalSpace(new OrbitalSpace(*h)), nelec(n), ms2(m2), symmetry(s)
{
}

State::State(OrbitalSpace& h, int n, int s, int m2)
  : orbitalSpace(new OrbitalSpace(h)), nelec(n), ms2(m2), symmetry(s)
{
}


void State::load(const Options &dump)
{
  nelec = dump.parameter("NELEC");
  ms2 = dump.parameter("MS2");
  symmetry = dump.parameter("ISYM",1)-1; // D2h symmetries are 0-7 internally, 1-8 externally
  // xout <<"nelec="<<nelec<<std::endl;
  //    xout <<"ms2="<<ms2<<endl;
  orbitalSpace.reset(new OrbitalSpace(dump));
  //    xout << "State::load orbitalSpace=" << (orbitalSpace != nullptr) << std::endl;
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
}


std::string State::str(int verbosity, unsigned int columns) const
{
  //    xout << "State::printable orbitalSpace=" << (orbitalSpace != nullptr) << verbosity << std::endl;
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
  std::ostringstream s;
  if (verbosity >= 0) {
    s<< "nelec="<<nelec<<" ms2="<<ms2<<" symmetry="<<symmetry+1;
    if (orbitalSpace!=nullptr)
      s << std::endl << "OrbitalSpace: " << orbitalSpace->OrbitalSpace::str(verbosity > 2 ? verbosity-1 : 0);
  }
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
  return s.str();
}
