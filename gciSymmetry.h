#ifndef GCISYMMETRY_H
#define GCISYMMETRY_H
#include <string>
#include <vector>
#include "gciPrintable.h"

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
    /*!
     * \brief Generate a printable representation of the object
     * \param verbosity How much to print
     * \return The string
     */
    std::string toString(int verbosity=0) const;
    /*!
     * \brief A string describing the object
     */
    std::string Title;
};

#endif // GCISYMMETRY_H
