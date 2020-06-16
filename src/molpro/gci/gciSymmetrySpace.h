#ifndef GCISYMMETRYSPACE_H
#define GCISYMMETRYSPACE_H
#include "molpro/gci/gciPrintable.h"
#include <string>
#include <vector>

namespace gci {

/*!
 * \brief General class to hold dimensions of symmetry blocks in a matrix or more general tensor.
 * The vector<> elements are the numbers of functions in each symmetry.
 */
class SymmetrySpace : public std::vector<size_t>, public gci::Printable {
public:
  /*!
   * \brief Construct a SymmetrySpace object
   * \param title String to associate with the object
   * \param maxrank Maximum rank of underlying tensor that will be addressed
   */
  explicit SymmetrySpace(std::string title = "", int maxrank = 2);
  /*!
   * \brief The maximum rank of tensor
   */
  int maxrank;
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const override;
  /*!
   * \brief A string describing the object
   */
  std::string Title;
  /*!
   * \brief offset of symmetry block for rank-1 tensor
   * \param sym1 symmetry of index in desired block
   * \return
   */
  size_t offset(int sym1) const;
  /*!
   * \brief offset of symmetry block for rank-2 tensor
   * \param sym1 symmetry of index pairs in desired block
   * \param sym2 symmetry of first index in desired block
   * \param parity 0 denotes general matrix, 1 symmetric, -1 antisymmetric
   * \return
   */
  size_t offset(int sym1, int sym2, int parity = 0) const;
  /*!
   * \brief calculateOffsets compute the offsets from the dimension
   */
  void calculateOffsets();
  /*!
   * \brief The total size of a rank-1 tensor
   * \return The size
   */
  size_t total() const;
  /*!
   * \brief The total size of a rank-2 tensor
   * \param sym1 The symmetry of the tensor
   * \param parity 0 denotes general matrix, 1 symmetric, -1 antisymmetric
   * \return The size
   */
  size_t total(int sym1, int parity = 0) const;

private:
  std::vector<size_t> offsets;
  std::vector<size_t> buffer;
};
} // namespace gci

#endif // GCISYMMETRYSPACE_H
