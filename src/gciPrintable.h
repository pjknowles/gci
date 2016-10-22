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

class Printable
{
public:
  Printable();
  /*!
     \brief
    printable form of the object.
     \return std::string
    */
  virtual std::string str(int verbosity=0) const=0;

  //    friend std::ostream& operator<<(std::ostream& os, Printable const& obj);

};

/*!
 \brief
 Overloaded <<
*/
std::ostream& operator<<(std::ostream& os, Printable const& obj);
}
using namespace gci;

#endif // GCIPRINTABLE_H
