#include "gciPrintable.h"
#include <sstream>
Printable::Printable()
{
}

std::ostream& gci::operator<<(std::ostream& os, Printable const& obj)
{
  return os << obj.str();
}
