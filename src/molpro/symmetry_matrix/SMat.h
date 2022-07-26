#ifndef SMAT_H
#define SMAT_H
#include <Eigen/Dense>
#include <climits>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>

#include <molpro/bytestreamC.h>
#include <molpro/memory.h>

namespace molpro {

/*!
 * \brief Dimension for each symmetry in one axis of an SMat
 */
typedef std::vector<size_t> dim_t;
/*!
 * \brief dims_t Dimensions for an entire SMat
 */
typedef std::vector<dim_t> dims_t;

/*!
 * \brief Index transposition symmetry of the SMat
 */
enum parity_t : int {
  parityNone = 0,       ///< Unsymmetric
  parityEven = 1,       ///< Symmetric
  parityOdd = -1,       ///< Antisymmetric, storage as for Symmetric
  parityOddPacked = -2, ///< Antisymmetric, diagonal elements not stored. Incompatible with Molpro matrices.
  parityUnspecified = 9 ///< Unspecified.
};

/*!
 * \brief The SMat class holds the meta-information and a pointer to data for a
 * symmetry-packed matrix object.  The \ref SMat is characterised by the symmetry, parity, and rank of the matrix, and
 * its row and column dimensions in each symmetry.
 *
 * Internally, the data are stored in a \ref molpro::array<T> with exactly the same layout as Molpro's symmetry-packed
 * matrices. There are facilities to construct an \ref SMat attached to an existing data buffer, and to return a pointer
 * to the buffer, to allow interoperability with legacy code. The intended normal usage mode is either (and preferably)
 * via the high-level functions and operators provided with the class, or using the symmetry blocks provided by \ref
 * block().
 *
 * Ranks 0 (vector) and 1 (matrix) are supported. In the case of matrices, the class stores a notion of whether the
 * matrix is to be used transposed or not; this allows transposition of the matrix without moving data. In addition, in
 * the case of symmetric or antisymmetric matrices, redundant symmetry blocks are not stored, and if they are requested,
 * the transpose block will be returned. Therefore all code that accesses symmetry blocks via \ref block() must enquire
 * via \ref block_transposed() or \ref block_dimensions() whether the block is transposed or not.
 *
 * Example:
 * \code
 *  dim_t dim; dim.push_back(3); dim.push_back(2);
 *  dims_t dims; dims.push_back(dim); dims.push_back(dim); // dimensions of a 3*3 + 2*2 matrix (if symmetry=0)
 *  SMat_<double> m(dims,-1,1); // antisymmetric matrix of symmetry 1
 *  std::cout << m.block_transposed(0) <<","<<m.block_transposed(1)<<std::endl; // returns 1,0 because only blocks with
 * sym(row) > sym(col) are actually stored. std::cout << "Matrix symmetry="<<m.symmetry()<<"; matrix
 * parity="<<m.parity()<<std::endl; std::cout << "block.dimensions(0): "<<m.block_dimensions(0)[0] <<
 * m.block_dimensions(0)[1] << m.block_dimensions(0)[2] << std::endl; std::cout << "block.dimensions(1):
 * "<<m.block_dimensions(1)[0] << m.block_dimensions(1)[1] << m.block_dimensions(1)[2] << std::endl; for (size_t k=0;
 * k<m.block(1).size(); k++) m.block(1)[k]=k+100.0; // some data for (unsigned int sym=0; sym<8; sym++) { size_t
 * rowstride=m.block_transposed(sym)?m.dimension(sym^m.symmetry(),1):1; size_t
 * colstride=m.block_transposed(sym)?1:m.dimension(sym,0); double factor=m.block_transposed(sym)?m.parity():1; for
 * (size_t row=0; row<m.dimension(sym,0);row++) { for (size_t col=0; col<m.dimension(sym^m.symmetry(),1);col++)
 *              std::cout <<
 * "m("<<row<<"."<<sym<<","<<col<<"."<<(sym^m.symmetry())<<")="<<factor*(m.block(sym)[col*colstride+row*rowstride])<<";
 * "; std::cout << std::endl;
 *      }
 *  }
 * \endcode
 * @tparam T data type of contained matrix elements
 */
template <class T> class SMat_ {
public:
  //  typedef typename Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> > M;
  //  typedef typename Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
  //  > Mconst; typedef typename Eigen::Map<Eigen::VectorXd> V;
  typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned,
                              Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
      M;
  typedef typename Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Unaligned,
                              Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
      Mconst;
  typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> V;
  using value_type = T;
  using scalar_type = T; // TODO deduce properly for complex
  /*!
   * \brief Construct using explicit symmetry-block dimensions.
   * The data buffer is allocated and maintained internally.
   * \param dimensions Numbers of rows and columns in each symmetry block.
   * dimensions.size() indicates the desired rank (1 or 2).
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric, default 0.
   * \param symmetry Overall symmetry of matrix (0-7).
   * \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object.
   */
  explicit SMat_(dims_t dimensions, parity_t parity = parityNone, int symmetry = 0, bool diagonal = false,
                 std::string description = "");
  /*!
   * \brief Construct using explicit symmetry-block dimensions.
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. \param dimensions Numbers of rows and columns in each symmetry block. dimensions.size()
   * indicates the desired rank (1 or 2). \param parity Index permutation symmetry: 0=none, 1=symmetric,
   * -1=antisymmetric, default 0. \param symmetry Overall symmetry of matrix (0-7). \param diagonal Whether this is a
   * diagonal matrix \param description A string describing the object.
   */
  explicit SMat_(dims_t dimensions, molpro::array<value_type>& buffer, parity_t parity = parityNone, int symmetry = 0,
                 bool diagonal = false, std::string description = "");
  /*!
   * \brief Construct using explicit symmetry-block dimensions.
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. \param dimensions Numbers of rows and columns in each symmetry block. dimensions.size()
   * indicates the desired rank (1 or 2). \param parity Index permutation symmetry: 0=none, 1=symmetric,
   * -1=antisymmetric, default 0. \param symmetry Overall symmetry of matrix (0-7). \param diagonal Whether this is a
   * diagonal matrix \param description A string describing the object.
   */
  explicit SMat_(dims_t dimensions, value_type* buffer, parity_t parity = parityNone, int symmetry = 0,
                 bool diagonal = false, std::string description = "");

  /*!
   * \brief Construct using symmetry-block dimensions obtained from standard orbital space definitions.
   * The data buffer is allocated and maintained internally.
   * \param space Designator of which orbital space ('X', 'Y', 'Z' or blank (default)).
   * Two characters can be given to indicate a rank-2 matrix with possibly different spaces for rows and columns, e.g.
   * 'XY'. \param symmetry Overall symmetry of matrix (0-7). \param parity Index permutation symmetry: 0=none,
   * 1=symmetric, -1=antisymmetric, default 0. \param diagonal Whether this is a diagonal matrix \param description A
   * string describing the object.
   */
  explicit SMat_(std::string space, parity_t parity = parityNone, int symmetry = 0, bool diagonal = false,
                 std::string description = "");
  /*!
   * \brief Construct using symmetry-block dimensions obtained from standard orbital space definitions.
   * \param space Designator of which orbital space ('X', 'Y', 'Z' or blank (default)).
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. Two characters can be given to indicate a rank-2 matrix with possibly different spaces for
   * rows and columns, e.g. 'XY'. \param symmetry Overall symmetry of matrix (0-7). \param parity Index permutation
   * symmetry: 0=none, 1=symmetric, -1=antisymmetric, default 0. \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object.
   */
  explicit SMat_(std::string space, molpro::array<value_type>& buffer, parity_t parity = parityNone, int symmetry = 0,
                 bool diagonal = false, std::string description = "");
  /*!
   * \brief Construct using symmetry-block dimensions obtained from standard orbital space definitions.
   * \param space Designator of which orbital space ('X', 'Y', 'Z' or blank (default)).
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. Two characters can be given to indicate a rank-2 matrix with possibly different spaces for
   * rows and columns, e.g. 'XY'. \param symmetry Overall symmetry of matrix (0-7). \param parity Index permutation
   * symmetry: 0=none, 1=symmetric, -1=antisymmetric, default 0. \param description A string describing the object.
   * \param diagonal Whether this is a diagonal matrix
   */
  explicit SMat_(std::string space, value_type* buffer, parity_t parity = parityNone, int symmetry = 0,
                 bool diagonal = false, std::string description = "");
  /*!
   * \brief Construct \ref SMat templated by another \ref SMat. The actual data are not copied.
   * The new data buffer is allocated and maintained internally.
   * \param source the \ref SMat to copy.
   * \param symmetry Overall symmetry of matrix (0-7). If not given, \c source.symmetry is used
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric. If not given, \c source.parity is
   * used. \param rank 1=vector, 2=matrix, 3=matrix that is transpose of source. If not given, \c source.rank is used.
   * \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object. Defaults to source->m_description.
   */
  explicit SMat_(SMat_ const* source, parity_t parity = parityUnspecified, int symmetry = 9, unsigned int rank = 0,
                 bool diagonal = false, std::string description = "");
  /*!
   * \brief Construct \ref SMat templated by another \ref SMat. The actual data are not copied.
   * \param source the \ref SMat to copy.
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. \param symmetry Overall symmetry of matrix (0-7). If not given, \c source.symmetry is used
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric. If not given, \c source.parity is
   * used. \param rank 1=vector, 2=matrix, 3=matrix that is transpose of source. If not given, \c source.rank is used.
   * \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object. Defaults to source->m_description.
   */
  explicit SMat_(SMat_ const* source, molpro::array<value_type>& buffer, parity_t parity = parityUnspecified,
                 int symmetry = 9, unsigned int rank = 0, bool diagonal = false, std::string description = "");
  /*!
   * \brief Construct \ref SMat templated by another \ref SMat. The actual data are not copied.
   * \param source the \ref SMat to copy.
   * \param buffer A pointer to an existing buffer of the right size, which this object will attach to, instead of
   * creating its own buffer. \param symmetry Overall symmetry of matrix (0-7). If not given, \c source.symmetry is used
   * \param parity Index permutation symmetry: 0=none, 1=symmetric, -1=antisymmetric. If not given, \c source.parity is
   * used. \param rank 1=vector, 2=matrix, 3=matrix that is transpose of source. If not given, \c source.rank is used.
   * \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object. Defaults to source->description.
   */
  explicit SMat_(SMat_ const* source, value_type* buffer, parity_t parity, int symmetry = 9, unsigned int rank = 0,
                 bool diagonal = false, std::string description = "");

  /*!
   * \brief Copy constructor.
   * @param source from where to copy
   * @param option
   *  - 1: offline storage if possible
   *  - 2: distributed storage if possible
   *  If offline or distributed storage is used, most of the class functions that access data cannot be used. The
   * exceptions are axpy, dot, scal
   */
  SMat_(SMat_ const& source, int option = 0);

  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a bytestream produced by bytestream()
   * \param buffer Location of actual data. If provided, the new object will attach to this, otherwise
   * a new buffer will be constructed.
   */
  explicit SMat_(const char* dump, value_type* buffer = nullptr);
  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param bs The bytestream produced by bytestream()
   * \param buffer Location of actual data. If provided, the new object will attach to this, otherwise
   * a new buffer will be constructed.
   */
  explicit SMat_(const molpro::bytestream& bs, value_type* buffer = nullptr)
      : SMat_::SMat_((const char*)&(bs.data()[0]), buffer) {}

  ~SMat_();

  /*!
   * @brief Copy an SMat, if necessary converting parity.
   * If the copy is from square to triangular, the parity of the result is enforced, ie
   * half \c in plus or minus the transpose of \c in is constructed.
   * If the dimensions of source and destination do not match, only addresses valid in source and destination will be
   * copied. \param source The source of data.
   * @return A reference to the result
   */
  SMat_& operator=(SMat_ const& source) { return copy(source); }
  /*!
   * @brief Copy an SMat, if necessary converting parity, and with the possibility of a positive offset
   * in the source row and column indices.
   * If the copy is from square to triangular, the parity of the result is enforced, ie
   * half \c in plus or minus the transpose of \c in is constructed.
   * If the dimensions of source and destination do not match, only addresses valid in source and destination will be
   * copied. \param source The source of data.
   * @param sourceOffset The offset in source data
   * @return A reference to the result
   */
  SMat_& copy(SMat_ const& source, dims_t sourceOffset = {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}});

  /*!
   * \brief Copy part or all of another SMat object, corresponding to a reduced range of indices, into this one.
   * The copy range for either rows or columns is restricted so that it lies within both the source and destination.
   * \param source Another SMat of any dimensions
   * \param sourceOffset Row and column index offsets for source
   * \param offset Row and column index offsets for *this
   */
  void splice(SMat_ const& source, dims_t sourceOffset = {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
              dims_t offset = {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}});

  /*!
   * @brief Construct a new SMat as a copy of a submatrix of this one
   * @param dimensions dimensions of new matrix
   * @param offset offsets in this for the start of the result
   * \param description A string describing the object. Defaults to this->m_description.
   * @return  The new matrix
   */
  SMat_ slice(dims_t dimensions, dims_t offset = {{0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0}},
              std::string description = "");

  SMat_& operator*=(value_type a);
  const SMat_ operator*(value_type a) const;

  SMat_& operator+=(SMat_ const& other);
  SMat_& operator-=(SMat_ const& other);
  const SMat_ operator+(SMat_ const& other) const;
  const SMat_ operator-(SMat_ const& other) const;

  /*!
   * \brief Multiply two SMat objects with contraction such that the result is another SMat.
   * Depending on the ranks of the sources, the possible paradigms are
   * - this(i,j)*other(j,k)
   * - this(i,j)*other(j)
   * - this(j)*other(j,k)
   * - this(j)*other(k)
   * \param other The matrix to postmultiply with.
   * \return a new SMat containing the product.
   */
  const SMat_ operator*(SMat_ const& other) const;
  /*!
   * \brief Multiply two SMat objects with contraction such that the result is another SMat, which will be placed
   * in the current object, which must have the correct dimensions.
   * Depending on the ranks of the sources, the possible paradigms are
   * - this(i,k) = beta*this(i,k) + alpha*a(i,j)*b(j,k)
   * - this(i) = beta*this(i) + alpha*a(i,j)*b(j)
   * - this(k) = beta*this(k) + alpha*a(j)*b(j,k)
   * - this(j,k) = beta*this(j,k) + alpha*a(j)*b(k)
   * \param a The first matrix in the product.
   * \param b The second matrix in the product.
   * \param alpha Factor to multiply a*b.
   * \param beta The fraction of the current object to add to the product.
   */
  void multiply(SMat_ const& a, SMat_ const& b, value_type alpha = 1.0, value_type beta = 0.0);

  /*!
   * \brief Multiply two SMat objects with contraction such that the result is a scalar scalar_type.
   * Depending on the ranks of the sources, the possible paradigms are
   * - this(i,j)*other(j,i)
   * - this(j)*other(j)
   * \param other The matrix to postmultiply with.
   * \return the scalar contraction of the two matrices.
   */
  scalar_type operator&(SMat_ const& other) const;

  /*!
   * \brief Transpose the object. The implementation does not touch the data,
   * but internally marks the state as transposed, affecting the results of subsequent
   * calls to block_transposed() and block_dimensions().
   */
  void transpose();

  /*!
   * @brief Evaluate, including explicitly transposing if necessary
   */
  void eval();

  /*!
   * \brief Set to zero all elements of the object whose absolute value
   * is less than a specified threshold.
   * \brief cut The threshold. Defaults to machine epsilon.
   */
  void trim(scalar_type cut = 0);

  /*!
   * \brief Test whether two objects are equal
   * \param other
   * \return
   */
  bool operator==(SMat_ const& other) const;
  bool operator!=(SMat_ const& other) const;

  /*!
   * \brief Obtain a pointer to the raw data.
   *
   * The underlying data model is
   * - For parity()==0 (neither symmetric or antisymmetric matrices),
   *   + symmetry blocks appear in ascending order of, and are indexed by, the row index symmetry rs.
   *   + the symmetry of the column index is cs=rs^symmetry
   *   + within the symmetry block, the row index changes fastest, ie has stride 1, ie columns are stored contiguously.
   *   + If transposed(), the data layout is that of the transpose of the matrix.
   * - For parity()!=0 (symmetric or antisymmetric) and symmetry()==0,
   *   + symmetry blocks appear in ascending order of, and are indexed by, the row index symmetry rs.
   *   + within the symmetry block, addressing is triangular, i.e. k*(k+1)/2+l with k > l.
   * - For parity()!=0 (symmetric or antisymmetric) and symmetry()!=0,
   *   + only blocks with rs > cs=rs^symmetry are actually stored, with the row index changing fastest, and column index
   * having stride dimension(rs).
   *   + block_offset(rs) returns block_offset(cs) if rs < cs. In this case, the implication is then that the column
   * index changes fastest, i.e. row index stride is dimension(cs). This is signalled through block_transposed()
   * delivering true.
   *   + In the case of parity()<0, the data stored are related to the actual data by a factor of -1 if
   * block_transposed() != transposed(). \return A pointer to the data buffer.
   */
  molpro::array<value_type>* data() const;

  /*!
   * \brief Construction of an empty \ref SMat is allowed to support those clients that need it, but the result is not
   * meaningfully usable.
   */
  explicit SMat_() : SMat_(dims_t{{0}, {0}}) {}

  explicit SMat_(const std::map<size_t, T>& source) {
    throw std::logic_error("Attempt to call meaningless copy constructor");
  }

  std::map<size_t, T> select_max_dot(size_t n, const SMat_<T>& y) const {
    throw std::logic_error("Attempt to call unimplemented select_max_dot");
  }

  std::map<size_t, T> select_max_dot(size_t n, const std::map<size_t, T>& y) const {
    throw std::logic_error("Attempt to call unimplemented select_max_dot");
  }

  SMat_<T>& operator=(const std::map<size_t, T>& source) {
    throw std::logic_error("Attempt to call unimplemented operator=");
  }

  /*!
   * \brief Return the buffer size of an \ref SMat object
   * \return
   */
  size_t size() const;

  /*!
   * \brief The maximum symmetry of the object, i.e. the symmetry index above which all dimensions are zero
   * and are impossible to reach by direct products.
   * \return
   */
  int max_symmetry() const;

  /*!
   * \brief The parity of the object.
   * \return
   */
private:
  void check_max_symmetry();

public:
  /*!
   * \brief return the symmetry of the object
   * \return
   */
  int symmetry() const;

  /*!
   * \brief The parity of the object.
   * \return
   */
  parity_t parity() const;

  /*!
   * \brief Whether the \ref SMat is transposed.
   * \return
   */
  bool transposed() const;

  /*!
   * @brief Whether the \ref SMat is diagonal
   * @return
   */
  bool Diagonal() const;

  /*!
   * \brief Return a dimension of a symmetry block
   * \param block_symmetry The symmetry of the row index
   * \param axis 0: number of rows ; 1: number of columns
   * \return
   */
  size_t dimension(unsigned int block_symmetry = 0, unsigned int axis = 0) const;

  /*!
   * \brief Return the rank of the object
   * \return
   */
  unsigned int rank() const;

  /*!
   * \brief Get the size of a symmetry block
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return The summed size of the symmetry block
   */
  size_t block_size(unsigned int block_symmetry) const;

  /*!
   * \brief Get a container pointing to a symmetry block
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return A vector object containing the data
   *
   */
  molpro::array<value_type> block(unsigned int block_symmetry) const;

#ifdef EIGEN_CORE_H
  /*!
   * \brief Get an Eigen Matrix mapping to a symmetry block.
   * All aspects of transposition are handled,
   * so that what is returned can be addressed using the (i,j) overload,
   * i running over block_symmetry, j over m_symmetry^block_symmetry.
   * If the result would need a parity factor, or triangular unpacking, an exception is thrown.
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return A Matrix Map object mapping the data
   */
  Mconst blockMap(unsigned int block_symmetry) const;
  M blockMap(unsigned int block_symmetry);
  bool blockNotEmpty(unsigned int block_symmetry) const;
  bool blockMapPossible(unsigned int block_symmetry) const;

  /*!
   * \brief Get an Eigen Matrix that is a copy of a symmetry block.
   * All aspects of transposition and triangular unpacking are handled,
   * so that what is returned can be addressed using the (i,j) overload,
   * i running over block_symmetry, j over m_symmetry^block_symmetry
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return A Matrix object containing the data
   */
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> blockCopy(unsigned int block_symmetry) const;

  void blockImport(unsigned int block_symmetry, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& e);

  /*!
   * \brief Get an Eigen Vector mapping to a symmetry block
   * \param block_symmetry The symmetry of the index of the desired block
   * \return A Vector object containing the data
   */
  V blockV(unsigned int block_symmetry) const;
#endif

private:
  /*!
   * \brief Get the offset of a symmetry block in the buffer
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return
   */
  size_t block_offset(unsigned int block_symmetry) const;

public:
  /*!
   * \brief Get the dimensions of a symmetry block in the buffer
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return The first two numbers returned are the numbers of rows and columns (Fortran convention).
   * The third number is zero if the rows are drawn from \c block_symmetry, or 1 if the columns are drawn from \c
   * block_symmetry, i.e. block_transposed().
   */
  std::vector<size_t> block_dimensions(unsigned int block_symmetry) const;

  dims_t dimensions() const;

  /*!
   * \brief Whether or not a particular symmetry block is stored as the transpose
   * \param block_symmetry The symmetry of the row (first) index of the desired block
   * \return
   */
  bool block_transposed(unsigned int block_symmetry) const;

  /*!
   * \brief Calculate the trace of a matrix
   * \return
   */
  value_type trace() const;

  /*!
   * \brief Calculate the norm. For vectors, it is the L2-norm; for matrices, the square root of the sum of the matrix
   * elements. \return
   */
  scalar_type norm() const;

  /*!
   * \brief Set the elements of the matrix to a value
   */
  void assign(value_type value);

  /*!
   * \brief Set equal to the unit matrix
   */
  void setIdentity();

  /*!
   * \brief Generate a printable summary of the object
   * \param title: Any desired title
   * \param level: Amount of information to report
   * \param precision How many digits to show
   * \return
   */
  std::string str(std::string title, int level = 0, int precision = Eigen::StreamPrecision) const;
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const { return str("", verbosity); }

  /*!
   * \brief scal: Scale a matrix by a constant
   * \param a: the constant factor
   * \param scaleDiagonal If false, then do not touch the diagonal elements. The matrix has to be square for this to
   * make sense.
   */
  void scal(value_type a, bool scaleDiagonal = true);

  /*!
   * @brief Assign a single scalar value to all elements of the matrix
   * @param a: the value to fill
   */
  void fill(value_type a) { assign(a); }

  /*!
   * \brief axpy Constant times a matrix plus a matrix, y=y+a*x
   * \param a Constant to multiply \c x
   * \param x Matrix to be added to \c y
   */
  void axpy(value_type a, const SMat_& x);

  /*!
   * \brief Perform a similarity transformation
   * \param q the transformation matrix
   * \param orthogonal Whether or not to assume that q is an orthogonal matrix
   * \return  transpose(q)*this*q
   */
  SMat_ transform(const SMat_& q, bool orthogonal = true) const;

  SMat_ diagonal() const;

  /*!
   * \brief Find the eigenvalues and, optionally, eigenvectors
   * \param val The eigenvalues
   * \param vec The right eigenvectors
   * \param vali The imaginary part of the eigenvalues. Not used if matrix is symmetric. Not required if it turns out
 that all the eigenvalues are real, but if supplied it will always be filled. If not supplied, but imaginary eigenvalues
 are found, an error is thrown.
   * \param vecl The left eigenvectors
   * \param algorithm 'lapack' or 'jacobi'
   * \param sort 'ascending', 'descending' or 'overlap'
   * \return If non-zero, an error code produced by Lapack dysev or dgeev
   */
  int ev(SMat_& val, SMat_* vec = nullptr, SMat_* vali = nullptr, SMat_* vecl = nullptr,
         std::string algorithm = "lapack", std::string sort = "ascending") const;

  /*!
   * @brief Orthonormalize the columns of the matrix
   * @param metric The metric, defaulting to unit matrix
   * @param algorithm
   */
  void orthogonalize(const SMat_* metric = nullptr, std::string algorithm = "Gram-Schmidt");
  /*!
   * \brief Calculate the exponential of a matrix.
   */
  SMat_ exp() const;

  /*!
   * \brief Calculate the logarithm of a matrix.
   */
  SMat_ log() const;

  /*!
   * \brief Calculate the square root of a matrix.
   */
  SMat_ sqrt() const;

  /*!
   * \brief Calculate a power of a matrix.
   */
  SMat_ pow(value_type p) const;

  /*!
   * \brief Calculate the inverse of a matrix.
   * \param SVThresh if specified, calculate the pseudoinverse via the singular-value decomposition of the matrix,
   * ignoring contributions from all singular values less than this parameter. \return The inverse
   */
  SMat_ inverse(scalar_type SVThresh = 0.0) const;

  /*!
   * \brief Solve one or more linear equation systems, where the current matrix is the coefficient matrix
   * \param rhs
   * \param algorithm One of the following, as described in the documentation for Eigen.
   * - "partialPivLu"
   * - "FullPivLU"
   * - "HouseholderQR"
   * - "ColPivHouseholderQR"
   * - "FullPivHouseholderQR"
   * - "LLT"
   * - "LDLT"
   * \return
   */
  SMat_ solve(const SMat_& rhs, std::string algorithm = "ColPivHouseholderQR") const;

  /*!
   * \brief Container for a singular-value decomposition of an SMat
   */
  class SVD {
  private:
    const SMat_& m_matrix;
    std::unique_ptr<SMat_> m_S;
    std::unique_ptr<SMat_> m_U;
    std::unique_ptr<SMat_> m_V;

  public:
    const SMat_& singularValues() { return *m_S; }
    const SMat_& matrixU() { return *m_U; }
    const SMat_& matrixV() { return *m_V; }
    explicit SVD(const SMat_& matrix, std::string algorithm = "BDC",
                 unsigned int computationOptions = Eigen::ComputeThinU | Eigen::ComputeThinV);
  };

public:
  /*!
   * \brief Produce a copy of the object which has no symmetry
   */
  SMat_ desymmetrise() const;

  /*!
   * \brief Serialise the object to a stream of bytes
   * \param data If true, write out the data buffer as well as the meta-information
   * \return the serialised representation of the object
   */
  class molpro::bytestream bytestream(bool data = true);

private:
  bool compatible(const SMat_& other, bool transpose = false) const {
    return rank() == other.rank() && parity() == other.parity() && m_dimensions == other.m_dimensions &&
           transposed() == (transpose ? (!other.transposed()) : other.transposed());
  }

  unsigned int max_symmetry_ = 8;
  dims_t m_dimensions; //!< numbers of rows and columns in each symmetry block
  enum parity_t m_parity;
  int m_symmetry;

public:
  /*!
   * \brief A string describing the object.
   */
  std::string m_description;

private:
  bool m_managed_buffer;
  bool m_transposed;
  bool m_diagonal; ///< if the matrix is diagonal
  molpro::array<T>* m_buffer;
  std::shared_ptr<molpro::array<T>> m_bufferp;

public:
  // implementation of LinearAlgebra::vector
  void axpy(value_type a, const std::map<size_t, value_type>& other) {
    throw std::logic_error("Unimplemented function");
  }
  scalar_type dot(const SMat_& other) const {
    if (other.m_parity == parityNone && m_parity != parityNone) {
      SMat_<T> x(m_dimensions, parityNone, m_symmetry, m_diagonal);
      x = *this;
      x.transpose();
      return other & x;
    } else if (m_parity == parityNone && other.m_parity != parityNone) {
      SMat_<T> x(m_dimensions, parityNone, m_symmetry, m_diagonal);
      x = other;
      x.transpose();
      return *this & x;
    } else if (&other == this and
               m_parity == parityNone) { // if the same object, transposition needs to be done carefully!
      auto x = other;                    // could do it without a copy
      x.transpose();
      return *this & x;
    } else {
      auto& x = const_cast<SMat_&>(other);
      x.transpose();
      auto result = *this & x;
      x.transpose();
      if (&other == this and (m_parity == parityOdd or m_parity == parityOddPacked)) // in that case
        return -result;
      else
        return result;
    }
  }

  scalar_type dot(const std::map<size_t, value_type>& other) const { throw std::logic_error("Unimplemented function"); }

  std::tuple<std::vector<size_t>, std::vector<value_type>> select(const molpro::array<value_type>& measure,
                                                                  const size_t maximumNumber = 1000,
                                                                  const scalar_type threshold = 0) const {
    throw std::logic_error("Unimplemented function");
  }
};

template <class T> inline SMat_<T> operator*(typename SMat_<T>::value_type a, SMat_<T> const& b) { return b * a; }

/*!
 * \brief Construct a new \ref SMat that is the transpose of the given one. The implementation does not touch the data,
 * but internally marks the state as transposed. No new buffer is allocated - instead, the new
 * object has a pointer to the data of mat.
 * \param mat The matrix to transpose.
 * \return The transposed matrix.
 */
template <class T> SMat_<T> transpose(const SMat_<T>& mat);

template <class T> SMat_<T> eval(const SMat_<T>& mat) {
  SMat_<T> result(mat);
  result = mat;
  result.eval();
  return result;
}

/*!
 * \brief Construct a new \ref SMat that is a numerically-truncated copy of the given one.
 * \param mat The matrix to trim.
 * \param cut The threshold. Defaults to machine epsilon.
 * \return The trimmed matrix.
 */
template <class T> SMat_<T> trim(const SMat_<T>& mat, typename SMat_<T>::scalar_type cut = 0);

/*!
 * \brief operator << Stream the default printable representation of an \ref SMat.
 * \param os
 * \param obj
 * \return
 */
template <class T> inline std::ostream& operator<<(std::ostream& os, SMat_<T> const& obj) { return os << obj.str(); }

// implementation only below here

template <class T> inline int SMat_<T>::symmetry() const { return this->m_symmetry; }

template <class T> inline parity_t SMat_<T>::parity() const { return this->m_parity; }

template <class T> inline bool SMat_<T>::transposed() const { return this->m_transposed; }

template <class T> inline bool SMat_<T>::Diagonal() const { return this->m_diagonal; }

template <class T> inline unsigned int SMat_<T>::rank() const { return this->m_dimensions.size(); }

template <class T> inline size_t SMat_<T>::dimension(unsigned int block_symmetry, unsigned int axis) const {
  if (axis > m_dimensions.size())
    throw std::runtime_error("dimension axis too large");
  return m_dimensions[axis][block_symmetry];
}

template <class T> inline dims_t SMat_<T>::dimensions() const { return m_dimensions; }

template <class T> inline int SMat_<T>::max_symmetry() const { return max_symmetry_; }

// free functions
template <class T> SMat_<T> transpose(const molpro::SMat_<T>& mat) {
  SMat_<T> result(mat);
  result.transpose();
  return result;
}

template <class T> SMat_<T> trim(const molpro::SMat_<T>& mat, typename SMat_<T>::scalar_type cut) {
  SMat_<T> result(mat);
  result.trim(cut);
  return result;
}

#ifndef DOXYGEN_SCAN
/*!
 * \private
 */
// obtain dimensions of named orbital spaces
namespace SymmetryMatrix {
dims_t spaces(std::string space);
void register_get_orbital_space(void (*func)(char c, size_t nt[]));
} // namespace SymmetryMatrix
#endif
using SMat = typename molpro::SMat_<double>;

} // namespace molpro

#endif // SMAT_H
