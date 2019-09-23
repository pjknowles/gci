#include "sharedCounter.h"
//TODO a cleaner separation of library-internal headers
#include "gci.h"
#include "gciRun.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory.h>
#include <unistd.h>
#include <cstring>
#include <fcntl.h>
#include <ga.h>
using namespace gci;


extern int gci::parallel_size;
extern int gci::parallel_rank;
extern bool gci::molpro_plugin;

#ifndef MOLPRO
#include <cerrno>
#include <cstring>
#include <sys/param.h>
#include "PluginGuest.h"
std::string get_working_path() {
  char temp[MAXPATHLEN];
  return (getcwd(temp, MAXPATHLEN) ? std::string(temp) : std::string(""));
}
int main(int argc, char *argv[])
//int main()
{
  char fcidumpname[1024] = "gci.fcidump";
  molpro_plugin = false;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &parallel_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &parallel_rank);
  GA_Initialize();
  if (parallel_rank > 0) freopen("/dev/null", "w", stdout);
  PluginGuest plugin("MOLPRO");
  if (plugin.active()) {
    if (!plugin.send("GIVE OPERATOR HAMILTONIAN FCIDUMP GCI")) throw std::logic_error("Unexpected plugin failure");
    strcpy(fcidumpname, plugin.receive().c_str());
  }
  if (argc > 1) strcpy(fcidumpname, argv[1]);

  size_t memory = 1000000000;
  for (int i = 2; i < argc; i++) {
    std::string s(argv[i]);
    size_t equals = s.find('=');
    if (equals != std::string::npos && s.substr(0, equals) == "MEMORY")
      memory = static_cast<size_t>(std::stol(s.substr(equals + 1)));
  }
  memory_initialize(memory);
  std::cout << "memory initialised to " << memory_remaining() << std::endl;
  size_t memory_allocated = memory_remaining();

  {
    Run run(fcidumpname);
    if (argc < 3 ||
        plugin.active()) {
      run.options.addParameter("METHOD", "DAVIDSON");
      // run.options.addParameter("PROFILER","0");
    } else
      for (int i = 2; i < argc; i++) {
        std::string s(argv[i]);
        size_t equals = s.find('=');
        if (equals != std::string::npos)
          run.options.addParameter(s.substr(0, equals), s.substr(equals + 1), true);
      }
    std::vector<double> e = run.run();
//  xout << "e after run:"; for (size_t i=0; i<e.size(); i++) xout <<" "<<e[i]; xout <<std::endl;
    if (plugin.active()) {
      // send the energy back
      if (plugin.send("TAKE PROPERTY ENERGY")) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(16);
        for (const double &ee : e)
          ss << ee << " ";
        plugin.send(ss.str());
      }
    }

    for (size_t state = 0; state < run.m_densityMatrices.size(); state++) {
      std::string densityname(fcidumpname);
      size_t pos = std::min(densityname.rfind('.'), densityname.size());
      densityname.replace(pos, densityname.size() - pos, ".density.");
      densityname += std::to_string(state + 1) + ".fcidump";
      FCIDump(run.m_densityMatrices[state],densityname);
      if (plugin.active()) { // send the density back
        if (plugin.send("TAKE DENSITY FCIDUMP"))
          plugin.send(densityname);
      }
    }
  }

  xout << "initial memory=" << memory_allocated << ", remaining memory=" << memory_remaining() << std::endl;
  if (plugin.active())
    plugin.send("");
  GA_Terminate();
  MPI_Finalize();
  return 0;
}
#endif
