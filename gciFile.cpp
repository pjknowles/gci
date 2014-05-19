#include "gciFile.h"
#include <iostream>

File::File()
{
  char *tmpname = strdup("tmpfileXXXXXX");
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
 f.seekg(address*8,std::ios_base::beg);
 f.read((char *) &buf[0],buf.size()*8);
}


void File::write(std::vector<double> &buf, size_t address)
{
 f.seekp(address*8,std::ios_base::beg);
 f.write((char *) &buf[0],buf.size()*8);
}
