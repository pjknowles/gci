#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include <vector>
#include "gci.h"
#include "gciOperator.h"
#include "gciState.h"
#include "gciStringSet.h"
#include "gciPrintable.h"
#include "gciFile.h"
#include "gciDeterminant.h"
#include "memory.h"
#include "SMat.h"
#include "LinearAlgebra.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State, public LinearAlgebra::vector<double>
{
public:
  /**
 * @brief
 *
 * @param filename is the file containing the FCIDUMP. If present, load is called.
 */
  Wavefunction(std::string filename="");
  /**
     * @brief
     *
     * @param dump points to an FCIdump object. If present, load is called.
     */
  Wavefunction(FCIdump* dump);

  /*!
     * \brief Construct a Wavefunction object linked to an OrbitalSpace
     * \param h The hamiltonian
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
  Wavefunction(OrbitalSpace *h, int nelec, int symmetry, int ms2);

  /*!
     * \brief Construct a Wavefunction object from a State prototype
     * \param state the State prototype
     */
  Wavefunction(const State& state);

  /*!
     * \brief Wavefunction copy constructor
     * \param other The object to be copied
     */
  Wavefunction( const Wavefunction& other);

  std::vector<StringSet> alphaStrings; ///< The alpha-spin strings defining the CI basis
  std::vector<StringSet> betaStrings; ///< The beta-spin strings defining the CI basis

  void allocate_buffer(); ///< allocate buffer to full size
  size_t size() const; ///< the size of the space

  void diagonalOperator(const Operator& oper); ///< set this object to the diagonal elements of the hamiltonian

  /*!
     * \brief find the index of the smallest component
     * \param n the n'th smallest will be found
     * \return offset in buffer
     */
  size_t minloc(size_t n=1);

  /*!
     * \brief Get a component of the wavefunction
     * \param offset which component to get
     * \return  the value of the component
     */
  double at(size_t offset);
  /*!
     * \brief get the determinant corresponding to a particular component of the wavefunction
     * \param offset which component to get
     * \return  the determinant
     */
  Determinant determinantAt(size_t offset);

  /*!
     * \brief Fill this object with the action of an operator on another wavefunction
     * \param h the operator
     * \param w the wavefunction
     * \param parallel_stringset whether to use parallel algorithm in StringSet construction
     */
  void operatorOnWavefunction(const Operator& h, const Wavefunction &w, bool parallel_stringset=false);

  /*!
   * \brief construct 1- and 2-particle density matrices
   * \param den1 The 1-particle density matrix
   * \param den2 The 2-particle density matrix
   */
  void density(memory::vector<double>& den1, memory::vector<double>& den2);
  void density(memory::vector<double>& den1);
  /*!
   * \brief construct 1- and 2-particle transition density matrices
   * \param den1 The 1-particle density matrix
   * \param den2 The 2-particle density matrix
   * \param bra Wavefunction defining the bra state (the ket is this)
   */
  void density(memory::vector<double>& den1, memory::vector<double>& den2, const Wavefunction& bra);
  void density(memory::vector<double>& den1, const Wavefunction& bra);
  /*!
   * \brief Calculate natural orbitals
   * \return
   */
  SMat naturalOrbitals();
private:
  void density(memory::vector<double>& den1, memory::vector<double>& den2, bool d1, bool d2, const Wavefunction& bra);
public:

  /*!
   * \brief blockOffset gives the address of the start of a symmetry block of the Wavefunction object
   * \param syma the symmetry of alpha strings
   * \return  the offset
   */

  size_t blockOffset(const unsigned int syma) const;

  std::string str(int verbosity=0) const;

  /*!
   * \brief axpy Add a multiple of another Wavefunction object to this one
   * \param a the factor defining the multiple
   * \param x the other wavefunction
   */
  void axpy(double a, const LinearAlgebra::vector<double> *x);
  void scal(double a);
    // Every child of ParameterVector needs exactly this
    Wavefunction* clone() const { return new Wavefunction(*this); }

  /*!
   * \brief push the object's buffer to a file
   * \param f the file
   * \param index where on the file, in units of the size of the object
   */
  void put(File& f, int index=0) override;

  /*!
   * \brief pull the local part of the object's buffer from a file
   * \param f the file
   * \param index where on the file, in units of the size of the object
   */
  void get(File& f, int index=0);

  /*!
   * \brief pull the object's buffer from a file
   * \param f the file
   * \param index where on the file, in units of the size of the object
   */
  void getAll(File& f, int index=0);

  /*!
   * \brief Construct a cumulative histogram of absolute values
   * \param edges the values defining bin edges
   * \param parallel whether to calculate the histogram in parallel
   * \param start first element of buffer
   * \param stop last element of buffer plus one
   * \return the numbers of coefficients whose absolute value is greater than the corresponding edge
   */
  std::vector<std::size_t> histogram(const std::vector<double>edges, bool parallel=true, std::size_t start=0, std::size_t stop=(size_t)(-1));

  /*!
   * \brief distributed whether the mode of calculation is to be distributed across processors
   */
  bool distributed;

  /*!
   * \brief gather give each process a full copy of buffer
   */
  void gather();

  //    Wavefunction& operator=(const double &value);
  void set(size_t offset, const double val);///< set one element to a scalar
  void set(const double val);///< set all elements to a scalar
  //    Wavefunction& operator=(const Wavefunction &other); ///< copy
//  namespace IterativeSolver {
  Wavefunction& operator*=(const double &value); ///< multiply by a scalar
  Wavefunction& operator+=(const Wavefunction &other); ///< add another wavefunction
  Wavefunction& operator-=(const Wavefunction &other); ///< subtract another wavefunction
  Wavefunction& operator-=(const double); ///< subtract a scalar from every element
  Wavefunction& operator-(); ///< unary minus
  Wavefunction& operator/=(const Wavefunction &other); ///< element-by-element division by another wavefunction
  double update(const Wavefunction &diagonalH, double & eTruncated, const double dEmax=(double)0); ///< form a perturbation-theory update, and return the predicted energy change. eTruncated is the energy change lost by truncation

  double norm(const double k=2); ///< k-norm
  /*!
   * \brief addAbsPower Evaluate this[i] += factor * abs(c[I])^k * c[I]
   * \param c
   * \param k
   * \param factor
   * \return a pointer to this
   */
  Wavefunction& addAbsPower(const Wavefunction &c, const double k=0, const double factor=1);

//  double* data() { return &buffer[0];}
//  const double* cdata() const { return &buffer[0];}

  friend class TransitionDensity;
  friend double operator*(const Wavefunction &w1, const Wavefunction &w2);///< inner product of two wavefunctions
  double dot(const LinearAlgebra::vector<double> *other) const;
  /*!
   * \brief this[i] = a[i]*b[i]
   * \param a
   * \param b
   */
  void times(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b);
  /*!
   * \brief this[i] = a[i]/(b[i]+shift)
   * \param a
   * \param b
   * \param shift
   * \param append whether to do += or =
   * \param negative whether =- or =+
   */
  void divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b, double shift=0, bool append=false, bool negative=false);
  void zero();
private:
  void buildStrings(); ///< build alphaStrings and betaStrings
  size_t dimension; ///< the size of the space
  std::vector<double> buffer; ///< buffer to hold coefficients describing the object
  std::vector<double>::iterator begin(); ///< beginning of this processor's data
  std::vector<double>::iterator end(); ///< end of this processor's data
  bool compatible(const Wavefunction &other) const; ///< whether this wavefunction is on the same space as another
  std::vector<size_t> _blockOffset;

};
double operator*(const Wavefunction &w1, const Wavefunction &w2);///< inner product of two wavefunctions
Wavefunction& operator+(const Wavefunction &w1, const Wavefunction &w2); ///< add two wavefunctions
Wavefunction& operator-(const Wavefunction &w1, const Wavefunction &w2); ///< subtract two wavefunctions
Wavefunction& operator/(const Wavefunction &w1, const Wavefunction &w2); ///< element-by-element division of two wavefunctions
Wavefunction& operator*(const Wavefunction &w1, const double &value);///< multiply by a scalar
Wavefunction& operator*(const double &value, const Wavefunction &w1);///< multiply by a scalar
}
using namespace gci;

#endif // GCIWAVEFUNCTION_H
