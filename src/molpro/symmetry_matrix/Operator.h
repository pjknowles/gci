#ifndef OPERATOR_H
#define OPERATOR_H
#include <array>
#include <complex>
#include <memory>
#include "SMat.h"
#include "SMatMat.h"

namespace molpro {

/*!
 * \brief Container for a quantum-mechanical operator (including density matrix)
 * described by zero-, one- and two-particle matrix elements. The operator is of the
 * form
 * \f[
 * \hat O = O + O^i_a a^\dagger i + \frac14 O^{ij}_{ab} a^\dagger b^\dagger j i
 * \f]
 *
 * The underlying basis set consists of alpha- and beta-spin orbitals whose numbers
 * within each symmetry type are identical, and given by dimension parameters on
 * construction. It is acceptable for range and domain to be drawn from different spaces, and for
 * particle interchange antisymmetry to be either not assumed or observed, i.e.
 * i,j not belonging, or belonging, to the same one-particle space (similarly for a,b).
 *
 * The operator is constructed with a parity for exchanging bra and ket that is -1, 0 or 1.
 * For the two-particle part this is specified separately for the first (ai) and second (bj)
 * indices separately, through the constructor parameter hermiticity. The overall operator parity,
 * and the parity of the one-particle part is the product of these.
 * If the parity is
 * non-zero, then the operator is assumed hermitian (\f$O^{ij}_{ab}=O_{ij}^{ab*}\f$); if it is negative,
 * the real stored values are to be
 * interpreted as the imaginary parts of the actual matrix elements.
 * For non-zero parity, the ket (i,j) and bra (a,b) dimensions must be the same.
 *
 * If the first element of construction parameter exchange is non-zero,
 * then it will be assumed that the
 * two bra indices a, b are drawn from the same space, and that the operator is symmetric (exchange=1)
 * or antisymmetric (exchange[0]=-1) with respect to swapping them,
 * with similar assumptions for i,j based on exchange[1].
 *
 * The one-particle part of the operator is stored as a single \ref SMat whose parity must
 * match that of the operator.
 *
 * The two-particle part of the operator is represented by \ref SMatMat objects,
 * in Mulliken order ([aibj]=\f$O^{ij}_{ab}\f$) and/or Dirac order([abij]=\f$O^{ij}_{ab}\f$).
 * In Mulliken order, the particle label swap symmetry [aibj]=[bjai] that exists for
 * non-zero exchange is not exploited,
 * i.e., the \ref SMatMat data always consists of a
 * rectangular matrix for each symmetry block.
 *
 * On construction, data structures for the matrix elements are allocated but not filled with values.
 * @tparam T data type of contained matrix elements
 */
template <class T> class Operator_ {
public:
  using SMat = typename molpro::SMat_<T>;
  using SMatMat = typename molpro::SMatMat_<SMat_<T>>;
  using value_type = T;
  using scalar_type = typename SMat::scalar_type;
  using M = typename SMat::M;
  using Mconst = typename SMat::Mconst;
  using V = typename SMat::V;
  /*!
   * \brief Construct an operator from its symmetry-block dimensions.
   * \param dimensions Numbers of one-particle functions in each symmetry block
   * for each of the two electron spins (alpha, beta), and each of the 2*rank two-particle indices (bra, ket, bra, ket).
   * In the Mulliken representation of two-electron operator, the storage is implemented in SMatMat, where
   * dimensions[][0] and dimensions[][1] specify the underlying SMat objects, and dimensions[][2] and dimensions[][3]
   * specify the dimensioning in the SMatMat container. In the Dirac representation, dimensions[][0] and dimensions[][2]
   * (bra) are for the inner SMat objects, and dimensions[][1] and dimensions[][3] (ket) for the container. \param rank
   * The rank (0, 1 or 2) of the operator. \param uhf Whether the underlying 1-particle spaces are different for alpha
   * and beta spin. \param hermiticity Conjugation symmetry, [ai] -> [ia], [bj] -> [jb]: 0=none, 1=symmetric,
   * -1=antisymmetric, default 0. \param exchange For each of bra, ket, whether there is an interchange symmetry between
   * two indices in the two-particle part of the operator: [abij] -> [baij]: 0=none, 1=symmetric, -1=antisymmetric,
   * default -1. \param symmetry Overall symmetry of operator (0-7). \param covariant Whether the operator transforms
   * like spatial basis functions, rather than like a density matrix \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object.
   */
  explicit Operator_(std::array<dims_t, 2> dimensions, int rank = 2, bool uhf = false,
                     std::vector<int> hermiticity = {1, 1}, std::vector<int> exchange = {-1, -1},
                     unsigned int symmetry = 0, bool covariant = true, bool diagonal = false,
                     std::string description = "");
  /*!
   * \brief Construct an operator from its symmetry-block dimensions.
   * \param dimensions Numbers of one-particle functions in each symmetry block
   * for each of the 2*rank two-particle indices (bra, ket, bra, ket).
   * \param rank The rank (0, 1 or 2) of the operator.
   * \param uhf Whether the underlying 1-particle spaces are different for alpha and beta spin.
   * \param hermiticity Conjugation symmetry, [ai] -> [ia], [bj] -> [jb]: 0=none, 1=symmetric, -1=antisymmetric, default
   * 0. \param exchange For each of bra, ket, whether there is an interchange symmetry between two indices in the
   * two-particle part of the operator: [abij] -> [baij]: 0=none, 1=symmetric, -1=antisymmetric, default -1. \param
   * symmetry Overall symmetry of operator (0-7). \param covariant Whether the operator transforms like spatial basis
   * functions, rather than like a density matrix \param diagonal Whether this is a diagonal matrix \param description A
   * string describing the object.
   */
  explicit Operator_(dims_t dimensions, int rank = 2, bool uhf = false, std::vector<int> hermiticity = {1, 1},
                     std::vector<int> exchange = {-1, -1}, unsigned int symmetry = 0, bool covariant = true,
                     bool diagonal = false, std::string description = "")
      : Operator_(std::array<dims_t, 2>{dimensions, dimensions}, std::move(rank), std::move(uhf),
                  std::move(hermiticity), std::move(exchange), std::move(symmetry), std::move(covariant),
                  std::move(diagonal), std::move(description)) {}
  /*!
   * \brief Construct a hermitian fermionic operator
   * \param dimension Numbers of orbitals in each symmetry block.
   * \param rank The rank (0, 1 or 2) of the operator.
   * \param uhf Whether the underlying 1-particle spaces are different for alpha and beta spin.
   * \param symmetry Overall symmetry of operator (0-7).
   * \param covariant Whether the operator transforms like spatial basis functions, rather than like a density matrix
   * \param diagonal Whether this is a diagonal matrix
   * \param description A string describing the object.
   */
  explicit Operator_(dim_t dimension, int rank = 2, bool uhf = false, unsigned int symmetry = 0, bool covariant = true,
                     bool diagonal = false, std::string description = "")
      : Operator_(dims_t{dimension, dimension, dimension, dimension}, rank, uhf, std::vector<int>{1, 1},
                  std::vector<int>{-1, -1}, symmetry, covariant, diagonal, description) {}
  /*!
   * \brief Copy constructor. A complete (deep) copy is made.
   * \param source Object to be copied.
   */
  Operator_(const Operator_<T>& source);

  virtual ~Operator_() = default;

  /*!
   * \brief Assigment operator
   * \param source Object to be copied.
   * \return A reference to this.
   */
  Operator_<T>& operator=(const Operator_<T>& source);

  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param dump The raw buffer of a bytestream produced by bytestream()
   */
  static class Operator_ construct(const char* dump);
  /*!
   * \brief Construct an object from what is produced by bytestream(). If the bytestream
   * contains data, it will be loaded, otherwise the contents of the object are undefined,
   * and only the dimensions and parameters are loaded.
   * \param bs The bytestream produced by bytestream()
   */
  static class Operator_ construct(class molpro::bytestream& bs);

  /*!
   * @brief Construct a new operator as a copy of a submatrix of this one
   * @param dimensions dimensions of new matrix. First index 0..1 alpha or beta; second index 0..3 bra1 ket1 bra2 ket2;
   * third index 0..7 symmetry.
   * @param offset offsets in *this for the start of the result. If omitted, zero; if alpha only, beta is the same.
   * \param description A string describing the object. Defaults to this->m_description.
   * @return  The new matrix
   */
  class Operator_ slice(std::array<dims_t, 2> dimensions, std::array<dims_t, 2> offset = {},
                        std::string description = "") const;

  /*!
   * @brief Construct a new operator as a copy of a submatrix of this one
   * @param dimensions dimensions of new matrix
   * @param offset offsets in this for the start of the result
   * \param description A string describing the object. Defaults to this->m_description.
   * @return  The new matrix
   */
  class Operator_ slice(dims_t dimensions,
                        dims_t offset = {{0, 0, 0, 0, 0, 0, 0, 0},
                                         {0, 0, 0, 0, 0, 0, 0, 0},
                                         {0, 0, 0, 0, 0, 0, 0, 0},
                                         {0, 0, 0, 0, 0, 0, 0, 0}},
                        std::string description = "") const {
    return slice({dimensions, dimensions}, {offset, offset}, description);
  }

  /*!
   * \brief Obtain a reference to 1-particle matrix elements.
   * \param spinUp alpha or beta spin.
   * \return
   */
  const SMat& O1(bool spinUp = true) const;
  SMat& O1(bool spinUp = true) { return const_cast<SMat&>(static_cast<const Operator_*>(this)->O1(spinUp)); }
  /*!
   * \brief Obtain a read-only reference to 2-particle matrix elements.
   * \param mulliken Whether Mulliken or Dirac index ordering is wanted.
   * \param spinUp1 alpha or beta spin for first particle.
   * \param spinUp2 alpha or beta spin for second particle.
   * \return
   */
  const SMatMat& O2(bool spinUp1 = true, bool spinUp2 = true, bool mulliken = true, bool ensure = true) const;
  /*!
   * \brief Obtain a reference to 2-particle matrix elements.
   * Only Mulliken, not Dirac, index ordering can be obtained; use the const version of this function to get read-only
   * access to Dirac. \param spinUp1 alpha or beta spin for first particle. \param spinUp2 alpha or beta spin for second
   * particle. \return
   */
  SMatMat& O2(bool spinUp1 = true, bool spinUp2 = true) { // no write access to mulliken=false
    return const_cast<SMatMat&>(static_cast<const Operator_*>(this)->O2(spinUp1, spinUp2, true));
  }

  /*!
   * \brief Get a reference to a single matrix element
   * \param i
   * \param isym
   * \param j
   * \param jsym
   * \param k If less than zero, then a one-electron value
   * \param ksym
   * \param l
   * \param lsym
   * \param spinUp1
   * \param spinUp2
   * \param mulliken
   * \return
   */
  const value_type& element(int i, int isym, int j, int jsym, int k, int ksym, int l, int lsym, bool spinUp1 = true,
                            bool spinUp2 = true, bool mulliken = true) const;
  value_type& element(int i, int isym, int j, int jsym, int k, int ksym, int l, int lsym, bool spinUp1 = true,
                      bool spinUp2 = true);
  const value_type& element(int i, int isym, int j, int jsym, bool spinUp1 = true) const {
    return element(i, isym, j, jsym, -1, -1, -1, -1, spinUp1);
  }
  value_type& element(int i, int isym, int j, int jsym, bool spinUp1 = true) {
    return element(i, isym, j, jsym, -1, -1, -1, -1, spinUp1);
  }

public:
  /*!
   * \brief Generate a printable summary of the object
   * \param title: Any desired title
   * \param level: Amount of information to report
   * \return
   */
  std::string str(std::string title = "", int level = 1) const;

  /*!
   * \brief Return the dimension of a symmetry block of basis functions
   * \param symmetry The symmetry
   * \param axis 0: bra1 ; 1: ket1; 2: bra2; 3: ket2
   * \param spinUp alpha or beta spin.
   * \return
   */
  size_t dimension(unsigned int symmetry = 0, unsigned int axis = 0, bool spinUp = true) const {
    return m_dimensions[spinUp ? 0 : 1][axis][symmetry];
  }

  /*!
   * \brief Serialise the object to a stream of bytes
   * \return the serialised representation of the object
   */
  class molpro::bytestream bytestream() const;

  /*!
   * \brief Construct the Fock operator for this operator with a given density matrix
   * \param density The density matrix.
   * \param oneElectron Whether to use the one-electron part of this operator, or zero
   * \param description A string describing the object.
   * \return
   */
  class Operator_ fock(const class Operator_& density, bool oneElectron = true, std::string description = "") const;

  void zero();
  class Operator_& operator+=(const class Operator_& other);
  class Operator_& operator-=(const class Operator_& other);
  class Operator_& operator*=(value_type other);
  scalar_type operator&(const class Operator_& other) const;
  const class Operator_ operator+(const class Operator_& other) {
    Operator_ copy(*this);
    copy += other;
    return copy;
  }
  const class Operator_ operator-(const class Operator_& other) {
    Operator_ copy(*this);
    copy -= other;
    return copy;
  }
  const class Operator_ operator*(value_type other) {
    Operator_ copy(*this);
    copy *= other;
    return copy;
  }

  bool compatible(const class Operator_& other) const {
    return (m_rank < 2 || other.m_rank < 2) ? true : O2().compatible(other.O2()) && m_covariant == other.m_covariant;
  } // FIXME make stronger
  void checkCompatible(const class Operator_& other) const {
    if (!compatible(other))
      throw std::logic_error("Incompatible shapes in assigment of Operator_ object");
  }
  void mulliken_from_dirac();

  /*!
   * \brief Construct M(pq,rs)=this.mulliken(q,p,r,s)+delta(p,r)*this(q,s).
   * The composite indices pq, rs correspond to the memory layout of this(p,q)
   * Implementation defined only when p,q,r,s are in the same space.
   * Implementation defined only for a spin-summed operator ie m_uhf=false
   * \return
   */
  SMat metric();

  /*!
   * @brief generates dirac integrals from Mulliken, if they are out of date
   */
  void ensure_dirac() const;

  /*!
   * @brief Marks Dirac integrals as out of date
   */
  void set_dirty() const;
protected:
  size_t offset(int isym, int jsym, int ksym, int lsym, int i, int j, int k, int l, bool mulliken = true,
                bool spinUp1 = true, bool spinUp2 = true);

private:
  const std::array<dims_t, 2> m_dimensions;

public:
  const int m_rank; //! the rank of the operator
  const bool m_uhf; //!< Whether the underlying 1-particle spaces are different for alpha and beta spin.
  const std::vector<int> m_hermiticity; //!< Conjugation symmetry, [ai] -> [ia], [bj] -> [jb]: 0=none, 1=symmetric,
                                        //!< -1=antisymmetric, default 0.
  const std::vector<int> m_exchange;
  //!< For each of bra, ket, whether there is an interchange symmetry between two
  //! indices in the two-particle part of the operator: [abij] -> [baij]: 0=none, 1=symmetric, -1=antisymmetric, default
  //! -1.
  const unsigned int m_symmetry; //!< The spatial symmetry of the operator.
  std::string m_description;     //!< A string describing the object.
  value_type m_O0;               //!< The constant part of the operator.
  const bool
      m_covariant; //!< Whether the operator transforms like spatial basis functions, rather than like a density matrix
  const bool m_diagonal; //!< Whether the representation of the operator is diagonal (1-electron only)
protected:
  std::shared_ptr<SMat> m_O1a, m_O1b;
  std::shared_ptr<SMatMat> m_O2aa_mulliken, m_O2ab_mulliken, m_O2ba_mulliken, m_O2bb_mulliken;
  mutable std::shared_ptr<SMatMat> m_O2aa_dirac, m_O2ab_dirac, m_O2ba_dirac, m_O2bb_dirac;
  mutable bool m_dirac_out_of_date;
};
template <class T> inline std::ostream& operator<<(std::ostream& os, Operator_<T> const& obj) {
  return os << obj.str();
}

using Operator = typename molpro::Operator_<double>;

} // namespace molpro

#endif // OPERATOR_H
