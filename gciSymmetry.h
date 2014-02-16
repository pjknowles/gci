#ifndef GCISYMMETRY_H
#define GCISYMMETRY_H
#include <string>
#include <vector>

/*!
 * \brief General class to hold offsets of symmetry blocks in a matrix or more general tensor
 */
class SymmetryOffset :public std::vector<size_t> {
public:
    SymmetryOffset();
    /*!
     * \brief Construct a SymmetryOffset object
     * \param title String to associate with the object
     */
    SymmetryOffset(std::string title);
    /*!
     * \brief Generate a printable representation of the object
     * \param title Description, defaults to the object's title element
     * \return The string
     */
    std::string printable(std::string title="");
    /*!
     * \brief A string describing the object
     */
    std::string Title;
};

#endif // GCISYMMETRY_H
