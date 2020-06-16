#ifndef GCIPRINTABLE_H
#define GCIPRINTABLE_H

#include "molpro/gci/gci.h"
#include <ostream>
#include <string>

namespace molpro {
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
  virtual std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const = 0;
};

/*!
 \brief
 Overloaded <<
*/
inline std::ostream &operator<<(std::ostream &os, Printable const &obj) { return os << obj.str(); }
} // namespace gci
} // namespace molpro

#endif // GCIPRINTABLE_H
