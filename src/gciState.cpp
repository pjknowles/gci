#include "gciState.h"
#include "FCIdump.h"
#include <fstream>
#include <iostream>
#include <sstream>

State::State(std::string filename)
{
  orbitalSpace=NULL;
  if (filename!="") load(filename);
}

State::State(FCIdump& dump)
{
  orbitalSpace=NULL;
  load(dump);
}

State::State(OrbitalSpace *h, int n, int s, int m2)
{
  orbitalSpace=h;
  nelec = n;
  symmetry = s;
  ms2 = m2;
}

State::State(State *s)
{
  orbitalSpace = s->orbitalSpace;
  nelec = s->nelec;
  symmetry = s->symmetry;
  ms2 = s->ms2;
}

State::~State()
{
}


void State::load(std::string filename) {
  FCIdump dump(filename);
  load(dump);
}

void State::load(const FCIdump &dump)
{
  nelec = dump.parameter("NELEC").at(0);
  ms2 = dump.parameter("MS2").at(0);
  symmetry = dump.parameter("ISYM",std::vector<int>(1,1)).at(0)-1; // D2h symmetries are 0-7 internally, 1-8 externally
  // xout <<"nelec="<<nelec<<std::endl;
  //    xout <<"ms2="<<ms2<<endl;
  orbitalSpace = new OrbitalSpace(dump);
  //    xout << "State::load orbitalSpace=" << (orbitalSpace != NULL) << std::endl;
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
}

std::string State::str(int verbosity, unsigned int columns) const
{
  //    xout << "State::printable orbitalSpace=" << (orbitalSpace != NULL) << verbosity << std::endl;
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
  std::ostringstream s;
  if (verbosity >= 0) {
    s<< "nelec="<<nelec<<" ms2="<<ms2<<" symmetry="<<symmetry+1;
    if (orbitalSpace!=NULL)
      s << std::endl << "OrbitalSpace: " << orbitalSpace->OrbitalSpace::str(verbosity > 2 ? verbosity-1 : 0);
  }
  //    xout << "basisSize=" << orbitalSpace->basisSize <<std::endl;
  return s.str();
}
