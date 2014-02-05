#include "gciSymmetry.h"
#include <sstream>

SymmetryOffset::SymmetryOffset()
{
    resize(8,0);
}

std::string SymmetryOffset::printable()
{
   std::ostringstream s;
   for (iterator x=begin(); x!=end(); x++) {
       s << *x; if (x != end()-1) s << " ";
   }
   return s.str();
}
