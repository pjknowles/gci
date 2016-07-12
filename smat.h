#ifndef SMAT_H
#define SMAT_H
#include <string>
#include <vector>

/*!
 * \brief The smat class holds the meta-information and a pointer to data for a
 * symmetry-packed matrix object.  The \c smat is characterised by the symmetry, parity, and rank of the matrix, and its row and column dimensions in each symmetry.
 */
class smat
{
public:
  /*!
   * \brief smat
   * \param dimensions Numbers of rows and columns in each symmetry block
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric, default 0
   * \param symmetry Overall symmetry of matrix (0-7), or -1 in the case of a rank-1 matrix with entries in all symmetries, eg eigenvalues of a matrix
   * \param rank 1=vector 2=matrix
   * \param create 0 or 1: whether to allocate buffer memory (1 not compatible with buffer); -1 means create if and only if buffer is not passed
   * \param buffer The actual data
   * \param description A string describing the object
   */
  smat(std::vector<std::vector<size_t> > dimensions , int parity=0, int symmetry=0, unsigned int rank=2, int create=-1, double* buffer=NULL, std::string description="");

  /*!
   * \brief Construct \c smat templated by another \c smat. The actual data are not copied.
   * \param source the \c smat to copy.
   * \param symmetry Overall symmetry of matrix (0-7), or -1 in the case of a rank-1 matrix with entries in all symmetries, eg eigenvalues of a matrix. If not given, \c source.symmetry is used
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric. If not given, \c source.parity is used
   * \param rank 1=vector 2=matrix. If not given, \c source.rank is used
   * \param create 0 or 1: whether to allocate buffer memory (1 not compatible with buffer); -1 means create if and only if buffer is not passed
   * \param buffer The actual data
   * \param description A string describing the object. Defaults to source->description.
   */
  smat(smat const & source, int parity=9, int symmetry=9, unsigned int rank=0, int create=-1, double* buffer=NULL, std::string description="");

  /*!
   * \brief Construct an empty \c smat.
   */
  smat();

#ifdef MOLPRO
  /*!
   * \brief smat
   * \param space Designator of which orbital space ('X', 'Y', 'Z' or blank (default)).
   * Two characters can be given to indicate a rank-2 matrix with possibly different spaces for rows and columns, e.g. 'XY'
   * \param symmetry Overall symmetry of matrix (0-7), or -1 in the case of a rank-1 matrix with entries in all symmetries, eg eigenvalues of a matrix
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric, default 0
   * \param create 0 or 1: whether to allocate buffer memory (1 not compatible with buffer); -1 means create if and only if buffer is not passed
   * \param buffer The actual data
   * \param description A string describing the object
   */
  smat(std::string space, int parity=0, int symmetry=0, int create=-1, double* buffer=NULL, std::string description="");
#endif

  ~smat();

  /*!
   * \brief Return the buffer size of an \c smat object
   * \return
   */
  size_t size() const;

  /* /\*! */
  /*  * \brief Get a pointer to a symmetry block */
  /*  * \param block_symmetry The symmetry of the row (first) index of the desired block */
  /*  * \return */
  /*  *\/ */
  /* double* block(unsigned int block_symmetry) const; */

  //  /*!
  //   * \brief Get a 1-dimensional container of all the data in a symmetry block
  //   * \param block_symmetry The symmetry of the row (first) index of the desired block
  //   * \return A Standard Template Library container that points to the data. It should not be resized or
  //   * otherwise reallocated, but the data can be both read and written.
  //   */
  //  std::vector<double> * blockv(unsigned int block_symmetry) const;
  //
  //  /*!
  //   * \brief Get a 2-dimensional container of all the data in a symmetry block
  //   * \param block_symmetry The symmetry of the row (first) index of the desired block
  //   * \return A Standard Template Library container that points to the data. It should not be resized or
  //   * otherwise reallocated, but the data can be both read and written.
  //   */
  //  std::vector<std::vector<double> > * blockvv(unsigned int block_symmetry) const;

  /*!
   * \brief Get the offset of a symmetry block in the buffer
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return
   */
  size_t block_offset(unsigned int block_symmetry) const;

  /*!
   * \brief Get the dimensions of a symmetry block in the buffer
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return The first two numbers returned are the numbers of rows and columns (Fortran convention). The third number is zero if the rows are drawn from \c block_symmetry, or 1 if the columns are drawn from \c block_symmetry
   */
  std::vector<size_t> block_dimensions(unsigned int block_symmetry) const;

  /*!
   * \brief Calculate the trace of a matrix
   * \return
   */
  double trace() const;

  /*!
   * \brief Set the elements of the matrix to zero
   */
  void zero();

  /*!
   * \brief Generate a printable summary of the object
   * \param title: Any desired title
   * \param number: Any desired number to identify the print
   * \return
   */
  std::string str(std::string title="", int number=999999) const;

  /*!
   * @brief Copy an smat, if necessary converting parity.
   * If the copy is from square to triangular, the parity of the result is enforced, ie
   * half \c in plus or minus the transpose of \c in is constructed.
   * Even if no conversions are done, the routine does a 'deep' copy, in other words allocating
   * a new buffer if necessary. A shallow copy, where only the pointer to the buffer is copied, can be achieved
   * by using a simple assignment statement instead.
   * \param source The source of data.
   * \param parity Force parity of result
   * \param description If present, replace this.description.
   */
  void copy(smat const & source, int parity=999999, std::string description="");

  /*!
   * \brief scal: Scale a matrix by a constant
   * \param a: the constant factor
   */
  void scal(double a);

  /*!
   * \brief axpy Constant times a matrix plus a matrix, y=y+a*x
   * \param a Constant to multiply \c x
   * \param x Matrix to be added to \c y
   */
  void axpy(double a, smat& x);

  /*!
   * \brief ger Rank 1 update, a = a + alpha * x * y^T
   * \param alpha Constant to multiply vector outer product
   * \param x vector to be added to \c a
   * \param y vector to be added to \c a
   */
  void ger(double alpha, smat& x , smat& y );

  /*!
   * \brief gemm
   * \param a
   * \param b
   * \param transa
   * \param transb
   * \param alpha
   * \param beta
   */
  void gemm(smat const & a , smat const & b , char transa='N', char transb='N', double alpha=1.0, double beta=1.0);

  /*!
   * \brief gemv
   * \param a
   * \param x
   * \param transa
   * \param alpha
   * \param beta
   */
  void gemv(smat const & a , smat const & x , char transa='N', double alpha=1.0, double beta=1.0);

  /*!
   * \brief Find the eigenvalues and, optionally, eigenvectors
   * \param val The eigenvalues
   * \param vec The right eigenvectors
   * \param vali The imaginary part of the eigenvalues. Not used if matrix is symmetric. Not required if it turns out that all the eigenvalues are real, but if supplied it will always be filled.
 If not supplied, but imaginary eigenvalues are found, an error is thrown.
   * \param vecl The left eigenvectors
   * \param algorithm 'lapack' or 'jacobi'
   * \param sort 'ascending', 'descending' or 'overlap'
   * \return If non-zero, an error code produced by Lapack dysev or dgeev
   */
  int ev(smat & val, smat* vec=NULL, smat* vali=NULL, smat* vecl=NULL, const std::string algorithm="lapack", const std::string sort="ascending") const;

  /*!
   * \brief Produce a copy of the object which has no symmetry
   */
  smat desymmetrise() const;

private:
  void initialise(std::vector<std::vector<size_t> > dimensions , int symmetry=0, int parity=0, unsigned int rank=2, int create=-1, double* buffer=NULL, std::string description="");
  void multiply();
  void ensure_buffer();
  bool dimdiff(std::vector<size_t> const v1, std::vector<size_t> const v2) const;

  size_t buffer_size;
  std::string description;
  unsigned int rank;
  int symmetry;
  int parity;
 public:
  double* buffer;
  std::vector<std::vector<size_t> > dimensions; //!< numbers of rows and columns in each symmetry block
  std::vector<double*> block; //!< pointers to the symmetry blocks

};

#ifndef DOXYGEN_SCAN
extern "C" { void smat_get_orbital_space(char c, size_t nt[]); }
#endif

#endif // SMAT_H
