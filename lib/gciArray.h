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
    int m_comm_rank; //!< rank in process group
    int m_comm_size; //!< size of process group
protected:
    size_t m_dimension; //!< Overall dimension of the direct product Fock space
    int m_ga_handle; //!< Global Array handle, needed by GA libary
    int m_ga_pgroup; //!< Global Array processor group handle
    int m_ga_chunk; //!< GA chunck size
    bool m_ga_allocated; //!< Flags that GA has been allocated
public:

    Array() : m_communicator(), m_comm_rank(0), m_comm_size(0), m_dimension(0), m_ga_handle(0),
              m_ga_pgroup(0), m_ga_chunk(1), m_ga_allocated(false) { }

    explicit Array(MPI_Comm comm);

    Array(size_t dimension, MPI_Comm commun);

    Array(const Array &source);

    ~Array();

    bool empty() const; //!< check if GA has been allocated

    void allocate_buffer(); //!< allocates GA buffer. Blocking, collective operation.

    /*!
     * @brief Duplicates GA buffer
     * Requires communicators to be the same. Blocking, collective operation
     */
    void copy_buffer(const Array &source);

    void sync();//!< synchronises all processes in this group and ensures GA operations completed
    double at(size_t offset) const; //!< get element at the offset. Blocking with one-sided communication.
    void zero(); ///< Set all elements to zero. Collective communication
    void set(double val);///< set all elements to a scalar. Collective communication
    // Use put instead.
    void set(size_t ind, double val);///< set one element to a scalar

    /*!
     * @brief gets buffer[lo:hi] from global array (hi inclusive, i.e. not pythonic)
     * Blocking with one-sided communication.
     */
    std::vector<double> get(int lo, int hi);

    /*!
     * @brief  puts data array into GA's buffer[lo:hi] (hi inclusive, i.e. not pythonic)
     * Blocking with one-sided communication
     */
    void put(int lo, int hi, double *data);

    /*!
     * @brief gets elements with discontinuos indices from GA
     * Blocking with one-sided communication
     * @return res[i] = ga_buffer[indices[i]]
     */
    std::vector<double> gather(std::vector<int> &indices) const;

    /*!
     * @brief ga_buffer[indices[i]] = vals[i]
     * Puts vals of elements with discontinuos indices into GA.
     * Blocking with one-sided communication.
     */
    void scatter(std::vector<int> &indices, std::vector<double> &vals);

    /*!
     * @brief ga_buffer[indices[i]] += alpha * vals[i]
     * Accumulates vals of elements with discontinuos indices into GA.
     * Atomic, blocking, with on-sided communication
     */
    void scatter_acc(std::vector<int> &indices, std::vector<double> &vals, double alpha);

    /*!
       * \brief find the index of n smallest components
       * Collective operation, must be called by all processes in the group.
       * \param n number of smallest values to be found
       * \return offsets in buffer
       */
    std::vector<size_t> minlocN(size_t n = 1) const;

    /*!
     * @brief Copies GA buffer into a local vector
     * Blocking, one-sided communication.
     */
    std::vector<double> vec() const;


    size_t size() const {return m_dimension;} ///< total number of elements

    virtual bool compatible(const Array &other) const; ///< Checks that arrays of the same dimensionality

    /*!
     * \brief axpy Add a multiple of another Wavefunction object to this one
     * Blocking, collective communication
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
    void scal(double a); ///< Scale by a constant. Collective communication
    void add(const Array &other); ///< Add another array to this. Collective communication
    void add(double a); ///< Add a constant. Collective communication
    void sub(const Array &other); ///< Subtract another array from this. Collective communication
    void sub(double a); ///< Subtract a constant. Collective communication
    void recip(); ///< Take element-wise reciprocal of this. Collective communication

    /*!
     * @brief Scalar product with another array.
     * Collective communication. Both arrays should be part of the same processor group (same communicator).
     */
    double dot(const Array &other) const;

    double dot(const Array *other) const {
        return dot(*other);
    }

    double dot(const std::map<size_t, double> &other) const;

    Array &operator=(const Array &source) noexcept;
    Array &operator*=(double value); //!< multiply by a scalar. Collective communication
    Array &operator+=(const Array &other); //!< add another array. Collective communication
    Array &operator-=(const Array &other); //!< subtract another array. Collective communication
    Array &operator+=(double); //!< add a scalar to every element. Collective communication
    Array &operator-=(double); //!< subtract a scalar from every element. Collective communication
    Array &operator-(); //!< unary minus. Collective communication
    Array &operator/=(const Array &other); //!< element-by-element division. Collective communication

    void times(const Array *a, const Array *b); ///< this[i] = a[i]*b[i]. Collective communication

    /*!
     * \brief this[i] = a[i]/(b[i]+shift)
     * Collective communication.
     * \param a
     * \param b
     * \param shift
     * \param append Whether to do += or =
     * \param negative Whether - or +
     */
    void divide(const Array *a, const Array *b, double shift = 0, bool append = false, bool negative = false);
};  // class Array

double operator*(const Array &w1, const Array &w2);///< inner product of two wavefunctions. Collective communication
Array operator+(const Array &w1, const Array &w2); ///< add two wavefunctions. Collective communication
Array operator-(const Array &w1, const Array &w2); ///< subtract two wavefunctions. Collective communication
Array operator/(const Array &w1, const Array &w2); ///< element-by-element division. Collective communication
Array operator*(const Array &w1, const double &value);///< multiply by a scalar. Collective communication
Array operator*(const double &value, const Array &w1);///< multiply by a scalar. Collective communication


} // namespace gci
#endif //GCI_TENSOR_H
