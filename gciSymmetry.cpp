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

std::string SymmetryOffset::toString(int verbosity) const
{
    if (verbosity < 0) return std::string("");
    std::ostringstream s;
    s << Title;
    for (const_iterator x=begin(); x!=end(); x++) {
        s << " " << *x;
    }
    return s.str();
}
