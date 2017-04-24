#include "gci.h"
#include "gciFile.h"
#include "gciOperator.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciStringSet.h"
#include "gciExcitationSet.h"
#include "gciOrbitalSpace.h"
#include "gciRun.h"
#include "FCIdump.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include "memory.h"
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
using namespace gci;



#ifdef MOLPRO
#ifdef __cplusplus
extern "C" {
#endif
  void gcirun(double* energies, int nenergies, char* fcidump) {
    Run run(fcidump);
    try {
    std::vector<double>e = run.run();
    for (int i=0; i < (nenergies > (int)e.size() ? (int)e.size() : nenergies); i++)
        energies[i]=e[i];
    }
    catch(char const* c) { xout << "caught error: " <<c<<std::endl;}
    return;
}
#ifdef __cplusplus
}
#endif
#endif

int gci::parallel_size=1;
int gci::parallel_rank=0;
bool gci::molpro_plugin=false;
MPI_Comm gci::molpro_plugin_intercomm=MPI_COMM_NULL;
sharedCounter* gci::_nextval_counter;
int64_t gci::__nextval_counter=0;
int64_t gci::__my_first_task=0;
int64_t gci::__task=0;
int64_t gci::__task_granularity=1;

#ifndef MOLPRO
#include <errno.h>
#include <string.h>
#include <sys/param.h>
#include "PluginGuest.h"
std::string get_working_path()
{
   char temp[MAXPATHLEN];
   return ( getcwd(temp, MAXPATHLEN) ? std::string( temp ) : std::string("") );
}
int main(int argc, char *argv[])
//int main()
{
  char fcidumpname[1024]="gci.fcidump";
  molpro_plugin=false;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&parallel_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&parallel_rank);
  gci::_nextval_counter= new sharedCounter();
  if (parallel_rank > 0) freopen("/dev/null", "w", stdout);
  PluginGuest plugin("MOLPRO");
  if (plugin.active()) {
      if (!plugin.send("GIVE OPERATOR HAMILTONIAN FCIDUMP GCI")) throw std::logic_error("Unexpected plugin failure");
      strcpy(fcidumpname,plugin.receive().c_str());
    }
  Run run(fcidumpname);
  if (argc<2 || 
     plugin.active()) {
    run.addParameter("METHOD","DAVIDSON");
    // run.addParameter("PROFILER","0");
  }
  else
    for (int i=1; i<argc; i++) {
      std::string s(argv[i]);
      size_t equals = s.find("=");
      if (equals != std::string::npos)
        run.addParameter(s.substr(0,equals),s.substr(equals+1),true);
  }
  memory_initialize(run.parameter("MEMORY",std::vector<int>{100000000})[0]);
  size_t memory_allocated=memory_remaining();
  std::vector<double> e=run.run();
  xout << "e after run:"; for (size_t i=0; i<e.size(); i++) xout <<" "<<e[i]; xout <<std::endl;

  if (plugin.active()) {
     // send the energy back
      if (plugin.send("TAKE PROPERTY ENERGY")) {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(16);
          for (auto ee=e.begin(); ee!=e.end(); ee++)
            ss << *ee << " ";
          plugin.send(ss.str());
     }
   }

  xout << "initial memory="<<memory_allocated<<", remaining memory="<<memory_remaining()<<std::endl;
  if (plugin.active())
    plugin.send("");
  delete _nextval_counter;
  MPI_Finalize();
  return 0;
}
#endif
