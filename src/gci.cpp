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
#ifdef GCI_PARALLEL
MPI_Comm gci::molpro_plugin_intercomm=MPI_COMM_NULL;
#endif
#ifdef MOLPRO
#include "cic/ItfMpp.h"
itf::FMppInt gci::mpp=itf::FMppInt(itf::FMppInt::MPP_NeedSharedFs|itf::FMppInt::MPP_GlobalDeclaration);
#else
sharedCounter* gci::_nextval_counter;
int64_t gci::__nextval_counter=0;
int64_t gci::__my_first_task=0;
int64_t gci::__task=0;
int64_t gci::__task_granularity=1;
#endif

#ifndef MOLPRO
int main(int argc, char *argv[])
//int main()
{
  char fcidumpname[1024]="gci.fcidump";
  molpro_plugin=false;
#ifdef GCI_PARALLEL
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&parallel_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&parallel_rank);
  gci::_nextval_counter= new sharedCounter();
  if (parallel_rank > 0) freopen("/dev/null", "w", stdout);
  MPI_Comm_get_parent(&molpro_plugin_intercomm);
  if (parallel_rank==0 && molpro_plugin_intercomm != MPI_COMM_NULL) {
      int length;
      MPI_Status status;
      // expect plugin server to identify itself
      MPI_Recv(&length,1,MPI_INT,0,0,molpro_plugin_intercomm,&status);
      char* id = (char*) malloc(length);
      MPI_Recv(id,length,MPI_CHAR,0,1,molpro_plugin_intercomm,&status);
//      printf("Plugin server: %s\n",id);
      molpro_plugin = !strncmp(id,"MOLPRO",6);
      if (molpro_plugin) {
          char molpro_version[5];
          strncpy(molpro_version,&id[7],4);
          fflush(stdout);
          int n = open(&id[12],O_WRONLY);
          dup2(n,1);
          close(n);
          printf("Plugin for Molpro version %s\n",molpro_version);
        }
    }
  MPI_Bcast(&molpro_plugin,1,MPI_INT,0,MPI_COMM_WORLD);
  if (molpro_plugin && parallel_rank==0) { // communication should be handled by just the root process, exchanging messages with the root process on the server
      // ask for an FCIDUMP
      char cmd[]="GIVE OPERATOR HAMILTONIAN FCIDUMP";
      int length=sizeof(cmd);
      MPI_Send(&length,1,MPI_INT,0,0,molpro_plugin_intercomm);
      MPI_Send(cmd,length,MPI_CHAR,0,1,molpro_plugin_intercomm);
      MPI_Status status;
      MPI_Recv(&length,1,MPI_INT,0,0,molpro_plugin_intercomm,&status);
      if (length==0) throw std::logic_error("plugin request has failed");
      MPI_Recv(fcidumpname,length,MPI_CHAR,0,1,molpro_plugin_intercomm,&status);
    }
  int length=strlen(fcidumpname);
  MPI_Bcast(fcidumpname,length,MPI_CHAR,0,MPI_COMM_WORLD);
#endif
  Run run(fcidumpname);
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
  memory_initialize(run.parameter("MEMORY",std::vector<int>{100000000})[0]);
  size_t memory_allocated=memory_remaining();
  std::vector<double> e=run.run();
  xout << "e after run:"; for (size_t i=0; i<e.size(); i++) xout <<" "<<e[i]; xout <<std::endl;

  if (molpro_plugin && parallel_rank==0) {
     // send the energy back
     char cmd[]="TAKE ENERGY";
     int length=sizeof(cmd);
     MPI_Send(&length,1,MPI_INT,0,0,molpro_plugin_intercomm);
     MPI_Send(cmd,length,MPI_CHAR,0,1,molpro_plugin_intercomm);
     int nstate=e.size();
     MPI_Send(&nstate,1,MPI_INT,0,2,molpro_plugin_intercomm);
     MPI_Status status;
     MPI_Recv(&length,1,MPI_INT,0,0,molpro_plugin_intercomm,&status);
     if (length) { // 'yes' answer received
       MPI_Send(&e[0],nstate,MPI_DOUBLE,0,3,molpro_plugin_intercomm);
     }
   }

  xout << "initial memory="<<memory_allocated<<", remaining memory="<<memory_remaining()<<std::endl;
#ifdef GCI_PARALLEL
  if (molpro_plugin && parallel_rank==0) {
      int signal=0;
      MPI_Send(&signal,1,MPI_INT,0,0,molpro_plugin_intercomm);
    }
  delete _nextval_counter;
  MPI_Finalize();
#endif
  return 0;
}
#endif
