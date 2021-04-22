// TODO a cleaner separation of library-internal headers
#include "molpro/gci/gci.h"
#include "molpro/gci/gciRun.h"
#include "molpro/gci/gciWavefunction.h"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <macdecls.h>
#include <molpro/memory.h>
#include <unistd.h>
#include <ga.h>

//#ifndef MOLPRO

#include <molpro/PluginGuest.h>

int main(int argc, char *argv[])
// int main()
{
  char fcidumpname[1024] = "gci.fcidump";
  molpro::gci::molpro_plugin = false;
  MPI_Init(&argc, &argv);
  GA_Initialize();
  molpro::gci::mpi_comm_compute = molpro::mpi::comm_global();
  MPI_Comm_size(molpro::gci::mpi_comm_compute, &molpro::gci::parallel_size);
  MPI_Comm_rank(molpro::gci::mpi_comm_compute, &molpro::gci::parallel_rank);
  if (molpro::gci::parallel_rank == 0)
    std::cout << "MPI_Comm_size = " << molpro::gci::parallel_size << std::endl;
  if (molpro::gci::parallel_rank > 0)
    freopen("/dev/null", "w", stdout);
  molpro::PluginGuest plugin("MOLPRO", molpro::gci::mpi_comm_compute);
  if (plugin.active()) {
    if (!plugin.send("GIVE OPERATOR HAMILTONIAN FCIDUMP GCI"))
      throw std::logic_error("Unexpected plugin failure");
    strcpy(fcidumpname, plugin.receive().c_str());
  }
  if (argc > 1)
    strcpy(fcidumpname, argv[1]);

  size_t memory = 4000000000; // TODO reduce this when memory leak has been found
  size_t ga_memory = 500000000;
  for (int i = 2; i < argc; i++) {
    std::string s(argv[i]);
    size_t equals = s.find('=');
    if (equals != std::string::npos && s.substr(0, equals) == "MEMORY")
      memory = static_cast<size_t>(std::stol(s.substr(equals + 1)));
    if (equals != std::string::npos && s.substr(0, equals) == "GAMEMORY")
      ga_memory = static_cast<size_t>(std::stol(s.substr(equals + 1)));
  }
  memory_initialize(memory);
  MA_init(C_CHAR, 10000000, ga_memory);
  if (molpro::gci::parallel_rank == 0)
    std::cout << "memory initialised to " << memory_remaining() << std::endl;
  size_t memory_allocated = memory_remaining();

  {
    molpro::gci::Run run(fcidumpname);
    if (argc < 3 || plugin.active()) {
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
    //  cout << "e after run:"; for (size_t i=0; i<e.size(); i++) cout <<" "<<e[i]; cout <<std::endl;
    if (plugin.active()) {
      // send the energy back
      if (plugin.send("TAKE PROPERTY ENERGY")) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(16);
        for (const double &ee : e) {
          ss << ee << " ";
        }
        plugin.send(ss.str());
      }
    }

    for (size_t state = 0; state < run.m_densityMatrices.size(); state++) {
      std::string densityname(fcidumpname);
      size_t pos = std::min(densityname.rfind('.'), densityname.size());
      densityname.replace(pos, densityname.size() - pos, ".density.");
      densityname += std::to_string(state + 1) + ".fcidump";
      molpro::gci::FCIDump(run.m_densityMatrices[state], densityname);
      if (plugin.active()) { // send the density back
        if (plugin.send("TAKE DENSITY FCIDUMP"))
          plugin.send(densityname);
      }
    }
  }

  cout << "initial memory=" << memory_allocated << ", remaining memory=" << memory_remaining() << std::endl;
  if (plugin.active())
    plugin.send("");
  GA_Terminate();
  MPI_Finalize();
  return 0;
}

//#endif
