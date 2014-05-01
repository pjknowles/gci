#ifndef GCISYMMETRY_H
#define GCISYMMETRY_H
#include <string>
#include <vector>
#include "gciPrintable.h"

namespace gci {

/*!
 * \brief General class to hold offsets of symmetry blocks in a matrix or more general tensor
 */
class SymmetryOffset :public std::vector<size_t>, public gci::Printable {
public:
    SymmetryOffset();
    /*!
     * \brief Construct a SymmetryOffset object
     * \param title String to associate with the object
     */
    SymmetryOffset(std::string title);
    std::string toString(int verbosity=0) const;
    /*!
     * \brief A string describing the object
     */
    std::string Title;
};
}

using namespace gci;

#endif // GCISYMMETRY_H
