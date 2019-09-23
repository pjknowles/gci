#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H

#include <vector>
#include <map>
#include <tuple>
#include <string>
#include "gci.h"
#include "gciState.h"
#include "gciStringSet.h"
#include "gciPrintable.h"
#include "gciFile.h"
#include "gciDeterminant.h"
#include "memory.h"
#include <SMat.h>
#include "gciOrbitals.h"

namespace gci {

/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State {
    friend class TransitionDensity;

    friend class MixedWavefunction;

public:
    using value_type = double;
    std::map<std::string, double> m_properties;
    std::map<size_t, double> buffer_sparse; ///< alternative storage to buffer, useful when very sparse
    bool m_sparse; ///< whether the coefficients are stored in buffer_sparse instead of buffer
    MPI_Comm m_communicator;
protected:
    memory::vector<double> buffer; ///< buffer to hold coefficients describing the object
    size_t dimension; ///< the size of the space
    std::vector<size_t> _blockOffset;
    int m_tilesize = -1, m_alphatilesize = -1, m_betatilesize = -1;
    static constexpr double m_activeStringTolerance = 1e-15;
public:
    /*!
       * \brief Construct a Wavefunction object linked to an OrbitalSpace
       * \param h The hamiltonian
       * \param nelec Number of electrons
       * \param symmetry Spatial symmetry
       * \param ms2 Sz quantum number times 2
       */
    Wavefunction(OrbitalSpace h, int nelec, int symmetry, int ms2, MPI_Comm communicator = MPI_COMM_COMPUTE);

    Wavefunction(OrbitalSpace *h, int nelec, int symmetry, int ms2, MPI_Comm communicator = MPI_COMM_COMPUTE);

    /*!
       * \brief Construct a Wavefunction object from a State prototype
       * \param state the State prototype
       */
    explicit Wavefunction(const State &state, MPI_Comm communicator = MPI_COMM_COMPUTE);

    Wavefunction(const Wavefunction &source, int option, MPI_Comm communicator = MPI_COMM_COMPUTE);

    std::vector<StringSet> alphaStrings; ///< The alpha-spin strings defining the CI basis
    std::vector<StringSet> betaStrings; ///< The beta-spin strings defining the CI basis

    void allocate_buffer(); ///< allocate buffer to full size
    size_t size() const {return m_sparse ? buffer_sparse.size() : dimension;} ///< the size of the space

    void diagonalOperator(
            const SymmetryMatrix::Operator &op); ///< set this object to the diagonal elements of the hamiltonian

    /*!
       * \brief find the index of the smallest component
       * \param n the n'th smallest will be found
       * \return offset in buffer
       */
    size_t minloc(size_t n = 1) const;

    /*!
       * \brief find the index of n smallest components
       * \param n number of smallest values to be found
       * \return offsets in buffer
       */
    std::vector<size_t> minlocN(size_t n = 1) const;

    /*!
       * \brief Get a component of the wavefunction
       * \param offset which component to get
       * \return  the value of the component
       */
    double at(size_t offset) const;

    /*!
       * \brief get the determinant corresponding to a particular component of the wavefunction
       * \param offset which component to get
       * \return  the determinant
       */
    Determinant determinantAt(size_t offset) const;

    /*!
     * \brief calculate address of one of the constituent strings in a Determinant
     * @param offset Determinant's address
     * @param axis
     * - 0: beta spin
     * - 1: alpha spin
     * @return string address
     */
    size_t stringAddress(size_t offset, unsigned int axis) const;

    /*!
     * \brief calculate symmetry of one of the constituent strings in a Determinant
     * @param offset Determinant's address
     * @param axis
     * - 0: beta spin
     * - 1: alpha spin
     * @return string symmetry
     */
    size_t stringSymmetry(size_t offset, unsigned int axis) const;

    /*!
       * \brief Add to this object the action of an operator on another wavefunction
       * \param h the operator
       * \param w the wavefunction
       * \param parallel_stringset whether to use parallel algorithm in StringSet construction
       */
    void
    operatorOnWavefunction(const SymmetryMatrix::Operator &h, const Wavefunction &w, bool parallel_stringset = false);

    /*!
     * \brief Construct a density matrix with this wavefunction
     * \param rank Can be 1 (1-particle density only) or 2 (1- and 2-particle densities)
     * \param uhf Whether to construct the spin-orbital or spin-summed density matrix
     * \param hermitian Whether to construct the average of the density and its hermitian conjugate
     * \param bra If given, the bra state to form a transition density with this as ket; otherwise this is bra and ket
     * \param description
     * \param parallel_stringset whether to use parallel algorithm in StringSet construction
     * \return
     */
    SymmetryMatrix::Operator density(int rank = 2,
                                     bool uhf = false,
                                     bool hermitian = true,
                                     const Wavefunction *bra = nullptr,
                                     std::string description = "",
                                     bool parallel_stringset = false) const;

    /*!
     * \brief Calculate natural orbitals
     * \return
     */
    Orbitals naturalOrbitals();

    /*!
     * \brief blockOffset gives the address of the start of a symmetry block of the Wavefunction object
     * \param syma the symmetry of alpha strings
     * \return  the offset
     */

    size_t blockOffset(unsigned int syma) const;

    std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const override;

    std::string values() const;

    /*!
     * \brief axpy Add a multiple of another Wavefunction object to this one
     * \param a the factor defining the multiple
     * \param x the other wavefunction
     */
    void axpy(double a, const Wavefunction &x);

    void axpy(double a, const std::shared_ptr<Wavefunction> &x) {
        axpy(a, *x);
    }

    void axpy(double a, const std::map<size_t, double> &x);

    std::tuple<std::vector<size_t>, std::vector<double> > select(const Wavefunction &measure,
                                                                 const size_t maximumNumber = 1000,
                                                                 const double threshold = 0) const;

    void scal(double a);

    /*!
     * \brief push the object's buffer to a file
     * \param f the file
     * \param index where on the file, in units of the size of the object
     */
    void putw(File &f, int index = 0);

    /*!
     * \brief pull the local part of the object's buffer from a file
     * \param f the file
     * \param index where on the file, in units of the size of the object
     */
    void getw(File &f, int index = 0);

    /*!
     * \brief pull the object's buffer from a file
     * \param f the file
     * \param index where on the file, in units of the size of the object
     */
    void getAll(File &f, int index = 0);

    /*!
     * \brief Construct a cumulative histogram of absolute values
     * \param edges the values defining bin edges
     * \param parallel whether to calculate the histogram in parallel
     * \param start first element of buffer
     * \param stop last element of buffer plus one
     * \return the numbers of coefficients whose absolute value is greater than the corresponding edge
     */
    std::vector<std::size_t> histogram(const std::vector<double> &edges,
                                       bool parallel = true,
                                       std::size_t start = 0,
                                       std::size_t stop = (size_t) (-1));

    /*!
     * \brief gather give each process a full copy of buffer
     */
    void replicate();

    void set(size_t offset, double val);///< set one element to a scalar
    void set(const double val);///< set all elements to a scalar
    Wavefunction &operator*=(const double &value); ///< multiply by a scalar
    Wavefunction &operator+=(const Wavefunction &other); ///< add another wavefunction
    Wavefunction &operator-=(const Wavefunction &other); ///< subtract another wavefunction
    Wavefunction &operator+=(double); ///< add a scalar to every element
    Wavefunction &operator-=(double); ///< subtract a scalar from every element
    Wavefunction &operator-(); ///< unary minus
    Wavefunction &operator/=(const Wavefunction &other); ///< element-by-element division by another wavefunction
    double update(const Wavefunction &diagonalH,
                  double &eTruncated,
                  const double dEmax = (double) 0); ///< form a perturbation-theory update, and return the predicted energy change. eTruncated is the energy change lost by truncation

    double norm(double k = 2); ///< k-norm
    /*!
     * \brief addAbsPower Evaluate this[i] += factor * abs(c[I])^k * c[I]
     * \param c
     * \param k
     * \param factor
     * \return a pointer to this
     */
    Wavefunction &addAbsPower(const Wavefunction &c, double k = 0, double factor = 1);


    friend double operator*(const Wavefunction &w1, const Wavefunction &w2);///< inner product of two wavefunctions
    double dot(const Wavefunction &other) const;

    double dot(const std::shared_ptr<Wavefunction> other) const {
        return dot(*other);
    }

    double dot(const std::unique_ptr<Wavefunction> other) const {
        return dot(*other);
    }

    double dot(const std::map<size_t, double> &other) const;

    /*!
     * \brief this[i] = a[i]*b[i]
     * \param a
     * \param b
     */
    void times(const Wavefunction *a, const Wavefunction *b);
    /*!
     * \brief this[i] = a[i]/(b[i]+shift)
     * \param a
     * \param b
     * \param shift
     * \param append whether to do += or =
     * \param negative whether =- or =+
     */
    void divide(const Wavefunction *a,
                const Wavefunction *b,
                double shift = 0,
                bool append = false,
                bool negative = false);
    void zero();

    void settilesize(int t = -1, int a = -1, int b = -1) {
        m_tilesize = t;
        m_alphatilesize = a;
        m_betatilesize = b;
    }

    std::vector<StringSet> activeStrings(bool spinUp = true) const;

protected:
    void buildStrings(); ///< build alphaStrings and betaStrings
    memory::vector<double>::iterator begin(); ///< beginning of this processor's data
    memory::vector<double>::iterator end(); ///< end of this processor's data
    memory::vector<double>::const_iterator cbegin() const; ///< beginning of this processor's data
    memory::vector<double>::const_iterator cend() const; ///< end of this processor's data
public:
    bool compatible(const Wavefunction &other) const; ///< whether this wavefunction is on the same space as another

};

double operator*(const Wavefunction &w1, const Wavefunction &w2);///< inner product of two wavefunctions
Wavefunction operator+(const Wavefunction &w1, const Wavefunction &w2); ///< add two wavefunctions
Wavefunction operator-(const Wavefunction &w1, const Wavefunction &w2); ///< subtract two wavefunctions
Wavefunction operator/(const Wavefunction &w1,
                       const Wavefunction &w2); ///< element-by-element division of two wavefunctions
Wavefunction operator*(const Wavefunction &w1, const double &value);///< multiply by a scalar
Wavefunction operator*(const double &value, const Wavefunction &w1);///< multiply by a scalar
}
#endif // GCIWAVEFUNCTION_H
