#include "common/molpro_config.h"
#include "gci.h"
#include "gciFile.h"
#include "gciMolpro.h"
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>

File::File()
{
  char tmpname[16];
  sprintf(tmpname,"tmp%6.6dXXXXXX",(int)parallel_rank);
  mkstemp(tmpname);
  f.open(tmpname,std::fstream::binary | std::fstream::in| std::fstream::out);
  remove(tmpname);
}

File::~File()
{
  f.close();
}

void File::read(std::vector<double> &buf, size_t address)
{
  read(&buf[0],buf.size(),address);
}

void File::read(double* buf, size_t length, size_t address)
{
 f.seekg(address*8,std::ios_base::beg);
 f.read((char *) &buf[0],length*8);
}


void File::write(std::vector<double> &buf, size_t address)
{
  write(&buf[0],buf.size(),address);
}

void File::write(double* buf, size_t length, size_t address)
{
    f.seekp(address*8,std::ios_base::beg);
    f.write((char *) &buf[0],length*8);
}
