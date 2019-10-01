#ifndef GCI_TENSOR_H
#define GCI_TENSOR_H

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <mpi.h>

//typedef struct _gci_mpi_comm_dummy *MPI_Comm;

namespace gci {

/*!
 * @brief Wrapper over globa array exposing functionality required by IterativeSolver and Array
 */
class Array {
public:
    using value_type = double;
    MPI_Comm m_communicator; //!< Outer communicator
protected:
    size_t m_dimension; //!< Overall dimension of the direct product Fock space
    int m_ga_handle; //!< Global Array handle, needed by GA libary
    int m_ga_pgroup; //!< Global Array processor group handle
    int m_ga_chunk; //!< GA chunck size
    bool m_ga_allocated; //!< Flags that GA has been allocated
public:

    Array() : m_communicator(nullptr), m_dimension(0), m_ga_handle(0), m_ga_pgroup(0), m_ga_chunk(1),
              m_ga_allocated(false) { }

    explicit Array(MPI_Comm comm) : m_communicator(comm), m_dimension(0), m_ga_handle(0), m_ga_pgroup(0), m_ga_chunk(1),
                                    m_ga_allocated(false) { }

    Array(size_t dimension, MPI_Comm commun);

    Array(const Array &source);

    ~Array();

    bool empty() const; //!< check if GA has been allocated

    void allocate_buffer(); //!< allocates GA buffer
    void copy_buffer(const Array &source); //!< duplicates GA buffer. Requires communicators to be the same

    double at(size_t offset) const; ///< get element at the offset
    std::vector<double> gather(std::vector<int> indices); ///< gets elements with indices from GA
    void scatter(std::vector<int> indices, std::vector<double> vals);///< puts vals of elements with indices into GA

    /*!
       * \brief find the index of n smallest components
       * \param n number of smallest values to be found
       * \return offsets in buffer
       */
    std::vector<size_t> minlocN(size_t n = 1) const;


    std::string str() const; ///< Writes elements of wavefunction vector into a string

    std::vector<double> vec() const; ///< Returns all coefficients in a single vector


    size_t size() const {return m_dimension;} ///< total number of elements

    bool compatible(const Array &other) const; ///< Checks that arrays of the same dimensionality

    /*!
     * \brief axpy Add a multiple of another Wavefunction object to this one
     * \param a the factor defining the multiple
     * \param x the other wavefunction
     */
    void axpy(double a, const Array &x);

    void axpy(double a, const Array *other) {
        axpy(a, *other);
    }

    void axpy(double a, const std::map<size_t, double> &x);


    /*!
     * @copydoc IterativeSolver::vector::scal
     */
    void scal(double a); ///< Scale by a constant
    void add(const Array &other); ///< Add another array to this
    void add(double a); ///< Add a constant
    void sub(const Array &other); ///< Subtract another array from this
    void sub(double a); ///< Subtract a constant
    void recip(); ///< Take element-wise reciprocal of this

    /*!
     * @copydoc IterativeSolver::vector::dot
     */
    double dot(const Array &other) const; ///< Take scalar product with another array
    double dot(const Array *other) const {
        return dot(*other);
    }

    double dot(const std::map<size_t, double> &other) const;///< Scalary product with a sparse array

    void zero(); ///< Set all elements to zero


    void set(double val);///< set all elements to a scalar
    void set(size_t ind, double val);///< set one element to a scalar
    Array &operator=(const Array &source) noexcept;
    Array &operator*=(double value); //!< multiply by a scalar
    Array &operator+=(const Array &other); //!< add another wavefunction
    Array &operator-=(const Array &other); //!< subtract another wavefunction
    Array &operator+=(double); //!< add a scalar to every element
    Array &operator-=(double); //!< subtract a scalar from every element
    Array &operator-(); //!< unary minus
    Array &operator/=(const Array &other); //!< element-by-element division

    void times(const Array *a, const Array *b); ///< this[i] = a[i]*b[i]

    /*!
     * \brief this[i] = a[i]/(b[i]+shift)
     * \param a
     * \param b
     * \param shift
     * \param append Whether to do += or =
     * \param negative Whether - or +
     */
    void divide(const Array *a,
                const Array *b,
                double shift = 0,
                bool append = false,
                bool negative = false);


};  // class Array

double operator*(const Array &w1, const Array &w2);///< inner product of two wavefunctions
Array operator+(const Array &w1, const Array &w2); ///< add two wavefunctions
Array operator-(const Array &w1, const Array &w2); ///< subtract two wavefunctions
Array operator/(const Array &w1, const Array &w2); ///< element-by-element division
Array operator*(const Array &w1, const double &value);///< multiply by a scalar
Array operator*(const double &value, const Array &w1);///< multiply by a scalar


} // namespace gci
#endif //GCI_TENSOR_H
