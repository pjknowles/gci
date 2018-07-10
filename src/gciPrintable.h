#ifndef GCIPRINTABLE_H
#define GCIPRINTABLE_H

#include <string>
#include <ostream>
#include "gci.h"

namespace gci {
/**
 * @brief
 * Base class for printable representation of object
 *
 */

class Printable {
 public:
  Printable() {}
  /*!
     \brief
    printable form of the object.
     \return std::string
    */
  virtual std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const =0;
};

/*!
 \brief
 Overloaded <<
*/
inline std::ostream &operator<<(std::ostream &os, Printable const &obj) { return os << obj.str(); }
}
using namespace gci;

#endif // GCIPRINTABLE_H
