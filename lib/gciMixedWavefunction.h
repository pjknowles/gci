#ifndef GCI_GCIMIXEDWAVEFUNCTION_H
#define GCI_GCIMIXEDWAVEFUNCTION_H

#include <vector>

#include "gciWavefunction.h"
#include "gciHProductSet.h"
#include "gciMixedOperator.h"
#include "gciMixedOperatorSecondQuant.h"

namespace gci {


/*!
 * @brief Mixed bosonic and fermionic wavefunction represented as configuration expansion in the direct product
 * space of Slater Determenants and Hartree Product of modals. It implements storage of the wavefunction
 * coefficients and application of the Hamiltonitan ( r = H c).
 *
 * FCIdump parameters:
 *      NMODE
 *          -- number of vibrational modes
 *          default = 0
 *      NMODAL
 *          -- number of modal basis functions per mode (for ground state NMODAL=1)
 *          default = 1
 *      VIB_EXC_LVL
 *          -- maximum number of modes simultaneously excited in the wavefunction
 *          default = 1
 *
 * Vibrational basis is specified by mode coupling level, CIS, CISD etc, and maximum excitation level.
 *
 * Nomenclature:
 *      modal - one particle vibrational basis function
 *
 * The full wavefunction is stored in the global array buffer. It is ordered by the vibrational basis, with each
 * vibrational basis having a corresponding full electronic CI vector.
 * Psi = {psi_{p,q}_0, psi_{p,q}_1A, psi_{p,q}_2A,..., psi_{p,q}_1B, ..., psi_{p,q}_1A_1B, ...}
 * where psi_{p,q} is the electronic CI wavefunction and indices IA imply direct product with HO basis function I of
 * mode A.
 *
 * Prallelism
 * ----------
 * The full wavefunction is stored in a global array. Relevant chunks are copied into the Wavefunction buffer
 * for computation. The buffer is copied (or accumulated etc) back into GA after which it is released.
 *
 *
 */
class MixedWavefunction {
public:
    using value_type = double;
    MPI_Comm m_head_communicator; //!< Outer communicator
    MPI_Comm m_child_communicator; //!< Communicator for children Wavefunction objects
protected:

    VibSpace m_vibSpace; //!< Parameters defining the vibrational space of current wavefunction
    HProductSet m_vibBasis; //!< Vibrational basis for the full space of current wavefunction
    size_t m_elDim; //!< Dimension of the electronic (slater determinant) space
    size_t m_dimension; //!< Overall dimension of the direct product Fock space

    /*!
     * @brief Full wavefunction.
     *
     * Psi = {psi_{p,q}_0, psi_{p,q}_1A, psi_{p,q}_2A,..., psi_{p,q}_1B, ..., psi_{p,q}_1A_1B, ...}
     * where psi_{p,q} is the electronic CI wavefunction and indices IA imply direct product with HO basis function I of
     * mode A.
     */
//    std::vector<Wavefunction> m_wfn;
    /*!
     * @brief Prototype electronic wavefunction
     *
     * It's buffer is populated with relevant section from GA before computation
     */
    Wavefunction m_prototype;
    int m_ga_handle; //!< Global Array handle, needed by GA libary
    int m_ga_chunk; //!< GA chunck size
    bool m_ga_allocated; //!< Flags that GA has been allocated
public:

    /*!
     * @brief Constructs the mixed wavefunction.
     */
    explicit MixedWavefunction(const Options &options, const State &prototype,
                               MPI_Comm head_commun = MPI_COMM_COMPUTE);

    MixedWavefunction(const MixedWavefunction &source, int option = 0, MPI_Comm head_commun = MPI_COMM_COMPUTE);

    ~MixedWavefunction();

    //! Flags if wavefunction vector has been allocated. Does not guarantee that Wavefunction buffers are allocated.
    bool empty() const; //!< flags if GA has been created

    void allocate_buffer(); //!< allocates GA buffer
    void copy_buffer(const MixedWavefunction &source); //!< duplicates GA buffer

    /*!
       * \brief find the index of n smallest components
       * \param n number of smallest values to be found
       * \return offsets in buffer
       */
    std::vector<size_t> minlocN(size_t n = 1) const;

    /*!
       * \brief Get a component of the wavefunction
       * \param offset Which component to get
       * \return  The value of the component
       */
    double at(size_t offset) const;

    /*!
     * @brief Returns a wavefunction corresponding to vibrational product under offset
     * @param iVib Index to the vibrational product
     * @return Reference to a wavefunction under offset
     */
    Wavefunction wavefunctionAt(size_t iVib) const;

    //! @brief Writes elements of wavefunction vector into a string
    std::string str() const;

protected:
    /*!
     * @brief Updates boundaries in GA for a block corresponding to electronic wavefunction under vibrational
     * index ``indKet``
     * @param iVib vibrational index of the electronic Wavefunction
     * @param lo index of the start of the block
     * @param hi index of the end of the block (inclusive)
     */
    static void ga_wfn_block_bound(int iVib, int *lo, int *hi, int dimension);
    static void ga_copy_to_local(int ga_handle, int iVib, Wavefunction &wfn, int dimension);
    static void ga_accumulate(int ga_handle, int iVib, Wavefunction &wfn, int dimension, int scaling_constant = 1.0);
public:

    /*!
     * \brief Add to this object the action of an operator on another wavefunction
     * \param ham Fully second quantized mixed Hamiltonian operator
     * \param w Other mixed Wavefunction
     * \param parallel_stringset whether to use parallel algorithm in StringSet construction
     */
    void operatorOnWavefunction(const MixedOperatorSecondQuant &ham, const MixedWavefunction &w,
                                bool parallel_stringset = false);

    /*!
     * @brief Set this object to the diagonal elements of the hamiltonian
     * \param ham Fully second quantized mixed Hamiltonian operator
     */
    void diagonalOperator(const MixedOperatorSecondQuant &ham, bool parallel_stringset = false);

    /*!
     * @return Returns all coefficients in a single vector
     */
    std::vector<double> vec() const;

    //! A copy of the vibrational space
    VibSpace vibSpace() const {return m_vibSpace;}

    size_t size() const {return m_dimension;}

    size_t elDim() const {return m_elDim;}

    size_t vibDim() const {return m_vibBasis.vibDim();}

public:
    /*!
     * @brief Checks that the two wavefunctions are of the same electronic State and of the same dimension.
     */
    bool compatible(const MixedWavefunction &other) const;

    /*! @copydoc IterativeSolver::vector::axpy
     */
    void axpy(double a, const MixedWavefunction &x);

    /*!
     * @copydoc IterativeSolver::vector::axpy(scalar,const vector<scalar>&)
     */
    void axpy(double a, const std::shared_ptr<MixedWavefunction> &other) {
        axpy(a, *other);
    }

    /*!
     * @copydoc IterativeSolver::vector::axpy(scalar,const std::map<size_t,scalar>& )
     */
    void axpy(double a, const std::map<size_t, double> &x) {
        throw std::logic_error("Cannot assume sparse vector is a map");
    }

    /*!
     * @copydoc IterativeSolver::vector::select
     * @todo Implement
     */
    std::tuple<std::vector<size_t>, std::vector<double>>
    select(const std::vector<double> &measure, const size_t maximumNumber = 1000,
           const double threshold = 0) const {return std::tuple<std::vector<size_t>, std::vector<double>>{{0}, {0}};};

    /*!
     * @copydoc IterativeSolver::vector::scal
     */
    void scal(double a);
    void add(const MixedWavefunction &other);
    void add(double a);
    void sub(const MixedWavefunction &other);
    void sub(double a);
    void recip();

    /*!
     * @copydoc IterativeSolver::vector::dot
     */
    double dot(const MixedWavefunction &other) const;

    /*!
     * @overload
     */
    double dot(const std::shared_ptr<MixedWavefunction> &other) const {
        return dot(*other);
    }

    /*!
     * @overload
     */
    double dot(const std::unique_ptr<MixedWavefunction> &other) const {
        return dot(*other);
    }

    /*!
     * @overload
     */
    double dot(const std::map<size_t, double> &other) const {
        throw std::logic_error("Cannot assume sparse vector is a map");
    }

    /*!
     * @copydoc IterativeSolver::vector::zero
     */
    void zero();

    /*!
     * @copydoc IterativeSolver::vector::clone
     * @todo change to managed pointer
     */
    MixedWavefunction *clone(int option = 0) const {return new MixedWavefunction(*this);}

    void set(double val);///< set all elements to a scalar
    void set(size_t ind, double val);///< set one element to a scalar
    MixedWavefunction &operator*=(const double &value); //!< multiply by a scalar
    MixedWavefunction &operator+=(const MixedWavefunction &other); //!< add another wavefunction
    MixedWavefunction &operator-=(const MixedWavefunction &other); //!< subtract another wavefunction
    MixedWavefunction &operator+=(double); //!< add a scalar to every element
    MixedWavefunction &operator-=(double); //!< subtract a scalar from every element
    MixedWavefunction &operator-(); //!< unary minus
    MixedWavefunction &operator/=(const MixedWavefunction &other); //!< element-by-element division

    /*!
     * @brief form a perturbation-theory update, and return the predicted energy change.
     * @param diagonalH Diagonal of the Hamiltonian matrix
     * @param eTruncated Energy change lost by truncation
     * @param dEmax
     * @return
     */
    double update(const Wavefunction &diagonalH,
                  double &eTruncated,
                  const double dEmax = 0.0) = delete;
    /*!
     * \brief this[i] = a[i]*b[i]
     * \param a
     * \param b
     */
    void times(const MixedWavefunction *a, const MixedWavefunction *b);

    /*!
     * \brief this[i] = a[i]/(b[i]+shift)
     * \param a
     * \param b
     * \param shift
     * \param append Whether to do += or =
     * \param negative Whether - or +
     */
    void divide(const MixedWavefunction *a,
                const MixedWavefunction *b,
                double shift = 0,
                bool append = false,
                bool negative = false);

//    std::map<std::string, double> m_properties;

    //! @todo Implement
    void settilesize(int t = -1, int a = -1, int b = -1) { };
};  // class MixedWavefunction

double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2);///< inner product of two wavefunctions
MixedWavefunction operator+(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< add two wavefunctions
MixedWavefunction
operator-(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< subtract two wavefunctions
MixedWavefunction
operator/(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< element-by-element division
MixedWavefunction operator*(const MixedWavefunction &w1, const double &value);///< multiply by a scalar
MixedWavefunction operator*(const double &value, const MixedWavefunction &w1);///< multiply by a scalar


}  // namespace gci
#endif //GCI_GCIMIXEDWAVEFUNCTION_H
