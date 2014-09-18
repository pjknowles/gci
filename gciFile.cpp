#include "gci.h"
#include "gciFile.h"
#include "gciMolpro.h"
#ifdef GCIMOLPROFILE
using namespace itf;
#else
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#endif
#ifdef GCIMOLPROFILE
int File::baseRecord=13000;
#ifdef __cplusplus
extern "C" {
#endif

void unusedRecord(int*, int*);

#ifdef __cplusplus
}
#endif
#endif

File::File()
{
#ifdef GCIMOLPROFILE

  int file=7;
  unusedRecord(&file,&baseRecord);
//  xout << "baseRecord="<<baseRecord<<std::endl;
  f = new FMolproStorageBlock(file,baseRecord++,0);
//  xout << "new file"<<std::endl;
#else
  char *tmpname = strdup("tmpfileXXXXXX");
  mkstemp(tmpname);
  f.open(tmpname,std::fstream::binary | std::fstream::in| std::fstream::out);
  remove(tmpname);
  free(tmpname);
#endif
}

File::~File()
{
#ifdef GCIMOLPROFILE
//  xout << "close file"<<std::endl;
  f->Delete();
  delete [] f;
  // need to delete Molpro records
#else
  f.close();
#endif
}

void File::read(std::vector<double> &buf, size_t address)
{
  if (parallel_rank==0) {
#ifdef GCIMOLPROFILE
//  xout << "read file"<<std::endl;
  f->Read(&buf[0], (FOffset)buf.size()*8,(FOffset)address*8);
#else
 f.seekg(address*8,std::ios_base::beg);
 f.read((char *) &buf[0],buf.size()*8);
#endif
  }
  //broadcast
#ifdef MOLPRO
  itf::GlobalBroadcast(&buf[0],buf.size(),0);
#elif GCI_PARALLEL
  int64_t type=1, root=0, size=buf.size();
  PPIDD_BCast(&buf[0],&size,&type,&root);
#endif
}


void File::write(std::vector<double> &buf, size_t address)
{
  if (parallel_rank==0) {
#ifdef GCIMOLPROFILE
//  xout << "write file"<<std::endl;
  f->Write(&buf[0], (FOffset)buf.size()*8,(FOffset)address*8);
#else
 f.seekp(address*8,std::ios_base::beg);
 f.write((char *) &buf[0],buf.size()*8);
#endif
}
}
