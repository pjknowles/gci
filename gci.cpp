#include "gci.h"
#include "gciFile.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciStringSet.h"
#include "gciExcitationSet.h"
#include "gciOrbitalSpace.h"
#include "gciRun.h"
#include "FCIdump.h"
#include <iostream>
#include <iomanip>
using namespace gci;



#ifdef MOLPRO
//#include "util/machines.h"
//#define GCI FORT_Extern(gci,GCI)
#ifdef __cplusplus
extern "C" {
#endif
  void gcirun(double* energies, int nenergies, char* fcidump) {
    Run run(fcidump);
    std::vector<double>e = run.run();
//    xout << "e after run:"; for (int i=0; i<e.size(); i++) xout <<" "<<e[i]; xout <<std::endl;
    for (int i=0; i < (nenergies > (int)e.size() ? (int)e.size() : nenergies); i++)
        energies[i]=e[i];
    return;
}
#ifdef __cplusplus
}
#endif
#endif


#ifndef MOLPRO
int main(int argc, char *argv[])
//int main()
{
  Run run("FCIDUMP");
  if (argc<2) {
    run.addParameter("METHOD","DAVIDSON");
    run.addParameter("PROFILER","0");
  }
  else
    for (int i=1; i<argc; i++) {
      std::string s(argv[i]);
      size_t equals = s.find("=");
      if (equals != std::string::npos)
        run.addParameter(s.substr(0,equals),s.substr(equals+1),true);
  }
  std::vector<double> e=run.run();
  xout << "e after run:"; for (size_t i=0; i<e.size(); i++) xout <<" "<<e[i]; xout <<std::endl;
  return 0;
}
#endif
