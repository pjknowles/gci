#ifndef SMATMAT_H
#define SMATMAT_H
#include "SMat.h"

namespace molpro {

/*!
 * \brief Class that contains a matrix of matrices.
 *
 * The data layout follows the loops
 * \code
 * size_t offset=0;
   for (unsigned int symij=0; symij<8; symij++) {
     symab = symij ^ symmetry();
     // block(symij) is a rectangular matrix, dimension m_templates[symab]->size(),block_size(symij)
     // can also be obtained as an Eigen::Map via blockM(symij)
     for (unsigned int symi=0; symi<8; symi++) {
       symj=symij^symi;
       if (parity()==0 || symj < symi) {
         for (size_t i=0; i<m_dimensions[0][symi]; i++)
           for (size_t j=0; j<((symi!=symj||parity()==0)m_dimensions[0][symj]:i+1); j++)
             // (*m_buffer)[offset] is the first data number in an SMat of symmetry symab indexed a,b holding abij
             // A pointer to the SMat can be obtained from smat(symij,symi,i,j)
             offset += m_templates[symab]->size();
       }
     }
   }
        \endcode
        \tparam The type of the contained matrices
 */
template <class T> class SMatMat_ {
public:
  using value_type = typename T::value_type;
  using scalar_type = typename T::scalar_type;
  using M = typename T::M;
  using Mconst = typename T::Mconst;
  using V = typename T::V;
  /*!
   * @brief Construct from a template contained matrix, and container dimensions
   * @param smat A prototype contained matrix. A copy is taken as a construction template, one for each symmetry
   * \param dimensions Numbers of rows and columns in each symmetry block.
   * dimensions.size() indicates the desired rank (1 or 2).
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric, default
   * 0. \param symmetry Overall symmetry of supermatrix (0-7), or -1 in the case of a rank-1 matrix with entries in all
   * symmetries \param description A string describing the object.
   */
  SMatMat_(const T& smat, dims_t dimensions, value_type* buffer = nullptr, parity_t parity = parityNone,
           int symmetry = 0, std::string description = "");
  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a bytestream produced by bytestream()
   * \param buffer Location of actual data. If provided, the new object will attach to this, otherwise
   * a new buffer will be constructed.
   */
  explicit SMatMat_(const char* dump, value_type* buffer = nullptr);
  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param bs The bytestream produced by bytestream()
   * \param buffer Location of actual data. If provided, the new object will attach to this, otherwise
   * a new buffer will be constructed.
   */
  explicit SMatMat_(const molpro::bytestream& bs, value_type* buffer = nullptr)
      : SMatMat_::SMatMat_((const char*)&(bs.data()[0]), buffer) {}
  SMatMat_(const SMatMat_& source);
  SMatMat_& operator=(const SMatMat_& other);
  SMatMat_& operator+=(const SMatMat_& other);
  SMatMat_& operator-=(const SMatMat_& other);
  SMatMat_& operator*=(typename T::value_type other);
  SMatMat_ operator+(const SMatMat_& other) {
    SMatMat_ copy(*this);
    return copy += other;
  }
  SMatMat_ operator-(const SMatMat_& other) {
    SMatMat_ copy(*this);
    return copy -= other;
  }
  SMatMat_ operator*(typename T::value_type other) {
    SMatMat_ copy(*this);
    return copy *= other;
  }
  /*!
   * \brief Form the trace-product of this operator with another, ie \f$\sum_{ijkl} A^{ij}_{kl} B^{ij}_{kl}\f$
   * \param other
   * \return
   */
  scalar_type operator&(const SMatMat_& other) const;

private:
  SMatMat_();

public:
  std::shared_ptr<T> smat(unsigned int ijsym, unsigned int isym, int i, int j = 0) const;

  /*!
   * \brief Return the buffer size of the object
   * \return
   */
  size_t size() const;

  /*!
   * \brief The symmetry of the object
   * \return
   */
  unsigned int symmetry() const;

  /*!
   * \brief The parity of the object.
   * \return
   */
  parity_t parity() const;

  /*!
   * \brief Return a dimension of a symmetry block
   * \param block_symmetry The symmetry of the row index
   * \param axis 0: number of rows ; 1: number of columns
   * \return
   */
  size_t dimension(unsigned int block_symmetry = 0, unsigned int axis = 0);

  /*!
   * \brief Return the rank of the object
   * \return
   */
  unsigned int rank() const;

  /*!
   * \brief Serialise the object to a stream of bytes
   * \param data If true, write out the data buffer as well as the meta-information
   * \return the serialised representation of the object
   */
  class molpro::bytestream bytestream(bool data = true);

  /*!
   * \brief Generate a printable summary of the object
   * \param title: Any desired title
   * \param level: Amount of information to report
   * \return
   */
  std::string str(std::string title = "", int level = 0) const;

  /*!
   * \brief scal: Scale a supermatrix by a constant
   * \param a: the constant factor
   * \param scaleDiagonal If false, then do not touch the diagonal elements of the underlying matrix components. The
   * matrices have to be square for this to make sense.
   */
  void scal(value_type a, bool scaleDiagonal = true);

  // private:
  /*!
   * \brief Get the size of a symmetry block specified by pair symmetry and row symmetry
   * \param ijsym The product symmetry of the row and column indices of the desired block
   * \param isym The symmetry of the row (first) index of the desired block
   * \return
   */
  size_t block_size(unsigned int ijsym, unsigned int isym) const;

  /*!
   * \brief Get the size of a symmetry block specified by pair symmetry
   * \param ijsym The product symmetry of the row and column indices of the desired block
   * \return The summed size of all blocks with the given index-pair symmetry
   */
  size_t block_size(unsigned int ijsym) const;

  /*!
   * \brief Get the offset of a symmetry block in the buffer
   * \param ijsym The product symmetry of the row and column indices of the desired block
   * \param isym The symmetry of the row (first) index of the desired block
   * \return
   */
  size_t block_offset(unsigned int ijsym, unsigned int isym) const;

  /*!
   * \brief Get the offset of a symmetry block specified by pair symmetry in the buffer
   * \param ijsym The product symmetry of the row and column indices of the desired block
   * \return
   */
  size_t block_offset(unsigned int ijsym) const;

public:
  /*!
   * \brief Get a container pointing to all the data in a symmetry block
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return A vector object containing the data
   *
   */
  molpro::array<value_type> block(unsigned int block_symmetry) const;

#ifdef EIGEN_CORE_H
  /*!
   * \brief Get an Eigen Matrix mapping to all the data in a symmetry block
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return A Matrix object containing the data
   */
  M blockM(unsigned int block_symmetry) const;
#endif

  bool compatible(const SMatMat_& other) const {
    return rank() == other.rank() && parity() == other.parity() && m_dimensions == other.m_dimensions &&
           m_symmetry == other.m_symmetry;
  }
  void checkCompatible(const SMatMat_& other) const {
    if (!compatible(other))
      throw std::logic_error("Incompatible shapes in assigment of SMatMat_ object");
  }

  molpro::array<value_type>* data() const { return m_buffer; }

  unsigned int max_symmetry_ = 8;
  dims_t m_dimensions; //!< numbers of rows and columns in each symmetry block
  parity_t m_parity;
  unsigned int m_symmetry;

public:
  /*!
   * \brief A string describing the object.
   */
  std::string m_description;

private:
private:
  bool m_managed_buffer;
  molpro::array<value_type>* m_buffer;
  std::shared_ptr<molpro::array<value_type>> m_bufferp;
  std::vector<std::shared_ptr<T>> m_templates;
  std::vector<value_type> m_zero_buffer;
};

/*!
 * \brief operator << Stream the default printable representation of an \ref SMatMat.
 * \param os
 * \param obj
 * \return
 */
template <class T> inline std::ostream& operator<<(std::ostream& os, SMatMat_<T> const& obj) { return os << obj.str(); }
template <class T> inline unsigned int SMatMat_<T>::rank() const { return this->m_dimensions.size(); }
template <class T> inline parity_t SMatMat_<T>::parity() const { return this->m_parity; }

using SMatMat = typename molpro::SMatMat_<molpro::SMat_<double>>;
} // namespace molpro

#endif // SMATMAT_H
