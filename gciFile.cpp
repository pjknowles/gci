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
//  char *tmpname = strdup("tmpfileXXXXXX");
  char tmpname[1];
  sprintf(tmpname,"tmp%6.6dXXXXXX",parallel_rank);
  mkstemp(tmpname);
//  xout << "tmpname="<<tmpname<<std::endl;
  f.open(tmpname,std::fstream::binary | std::fstream::in| std::fstream::out);
  remove(tmpname);
//  free(tmpname);
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
  read(&buf[0],buf.size(),address);
}

void File::read(double* buf, size_t length, size_t address)
{
#ifdef GCIMOLPROFILE
//  xout << "read file"<<std::endl;
  // Molpro I/O is collective
  f->Read(&buf[0], (FOffset)length*8,(FOffset)address*8);
#else
//  if (parallel_rank==0) {
 f.seekg(address*8,std::ios_base::beg);
 f.read((char *) &buf[0],length*8);
//  }
//  //broadcast
//#ifdef GCI_PARALLEL
//  int64_t type=1, root=0, size=buf.size();
//  PPIDD_BCast(&buf[0],&size,&type,&root);
//#endif
#endif
}


void File::write(std::vector<double> &buf, size_t address)
{
  write(&buf[0],buf.size(),address);
}

void File::write(double* buf, size_t length, size_t address)
{
#ifdef GCIMOLPROFILE
//  xout << "write file"<<std::endl;
  // Molpro I/O is collective
  f->Write(&buf[0], (FOffset)length*8,(FOffset)address*8);
#else
//  if (parallel_rank==0) {
    f.seekp(address*8,std::ios_base::beg);
    f.write((char *) &buf[0],length*8);
//  }
#endif
}
