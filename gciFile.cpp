#include "gciFile.h"
#ifdef GCIMOLPROFILE
#include "cic/ItfFortranInt.h"
#else
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#endif

File::File()
{
#ifdef GCIMOLPROFILE
  f = new FMolproStorageBlock(7,8000,0);
  xout << "new file"<<std::endl;
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
  xout << "close file"<<std::endl;
  f->Delete();
  delete [] f;
#else
  f.close();
#endif
}

void File::read(std::vector<double> &buf, size_t address)
{
#ifdef GCIMOLPROFILE
  xout << "read file"<<std::endl;
  f->Read(&buf[0], (FOffset)buf.size(),(FOffset)address);
#else
 f.seekg(address*8,std::ios_base::beg);
 f.read((char *) &buf[0],buf.size()*8);
#endif
}


void File::write(std::vector<double> &buf, size_t address)
{
#ifdef GCIMOLPROFILE
  xout << "write file"<<std::endl;
  f->Write(&buf[0], (FOffset)buf.size(),(FOffset)address);
#else
 f.seekp(address*8,std::ios_base::beg);
 f.write((char *) &buf[0],buf.size()*8);
#endif
}
