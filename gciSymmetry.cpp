#include "gciSymmetry.h"
#include <sstream>

SymmetryOffset::SymmetryOffset()
{
    resize(8,0);
    Title="";
}

SymmetryOffset::SymmetryOffset(std::string title)
{
    resize(8,0);
    Title=title;
}

std::string SymmetryOffset::toString(std::string title)
{
   std::ostringstream s;
   s << ((title=="") ? Title : title);
   for (iterator x=begin(); x!=end(); x++) {
       s << " " << *x;
   }
   return s.str();
}
