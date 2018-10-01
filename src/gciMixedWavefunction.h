#ifndef GCI_GCIMIXEDWAVEFUNCTION_H
#define GCI_GCIMIXEDWAVEFUNCTION_H

#include <vector>

#include <LinearAlgebra.h>
#include "gciWavefunction.h"

namespace gci {
/*!
 * @brief Mixed bosonic and fermionic wavefunction represented as configuration expansion in the direct product
 * space of Slater Determenants and Hartree Product of modals. It implements storage of the wavefunction
 * coefficients and application of the Hamiltonitan ( r = H c).
 *
 * Nomenclature:
 *      modal - one particle vibrational basis function
 *
 * Currently hardcoded for working with up to second order truncated Molecular Hamiltonian (BBO Hamiltonian).
 *
 * H = H_el + H_s + H_vib + H1_int + H2_int
 *
 * H_el -- BO electronic Hamiltonian
 *
 * H_s -- diagonal BO part of the Hamiltonian (purely electronic)
 *
 * H_vib -- vibrational Harmonic oscillator Hamiltonian
 *
 * H1_int = (d/dX_A H_el) dX_A + (m_dg n d_mn) d/dX_A
 * couples electronic degrees of freedom with a single vibrational mode.
 * Due to the nature of HO basis, only adjacent vibrational states are coupled
 * Requires use of vibrational basis with no mode coupling.
 *
 * H2_int = (d/dX_A d/dX_B H_el) dX_A dX_B
 * couples electronic degrees of freedom with two vibrational modes
 * Requires use of vibrational basis with two mode coupling.
 *
 * Wavefunction with only 1 mode coupling:
 * Psi = sum_A C_{p,q,I_A} |p,q> |I_A> |0_B> |0_C> ...
 *
 * In practice the coefficients C_{p,q,I_A} are stored as a vector<C_{p,q}>(n = nModes * nModals).
 * This class is in fact mostly a wrapper of gci::Wavefunction with a few hardcoded expressions for the vibrational
 * parts of the Hamiltonian.
 */
class MixedWavefunction : public LinearAlgebra::vector<double> {
    // Inhereting the vector class brings the functions that need to be implemented for linera algebra solver to work.
public:
    /*!
     * @brief Constructs the wavefunction using a State prototype for the electronic component.
     *
     * Should take some parameters defining the electronic state and excitation level of the vibrational basis,
     * together with a BBO Hamiltonian.
     *
     * BBO Hamiltonian determines whether two-mode coupling is required or not in the vibrational basis.
     */
    MixedWavefunction(const State &state, int nMode, int nModal, int modeCoupling = 1);
    MixedWavefunction(const MixedWavefunction &other, int option) {*this = other;}
    ~MixedWavefunction() override = default;

    size_t size() {return m_dimension;}
    void allocate_buffer(); //!< allocate buffer to full size

    /*!
       * \brief find the index of the smallest component
       * \param n the n'th smallest will be found
       * \return offset in buffer
       */
    size_t minloc(size_t n = 1) const;

    /*!
       * \brief Get a component of the wavefunction
       * \param offset which component to get
       * \return  the value of the component
       */
    double at(size_t offset) const;

    /*!
       * \brief Add to this object the action of an operator on another wavefunction
       * \param h the operator
       * \param w the wavefunction
       * \param parallel_stringset whether to use parallel algorithm in StringSet construction
       */
    void operatorOnWavefunction(const MixedOperator &h, const MixedWavefunction &w, bool parallel_stringset = false);

    void diagonalOperator(const Operator &op); //!< set this object to the diagonal elements of the hamiltonian

private:
    int m_nMode; //!< number of vibrational modes
    int m_nModal; //!< number of modals per mode
    int m_modeCoupling; //!< level of mode-mode coupling = 1, or 2
    size_t m_elecWfnSize; //!< size of the electronic Wavefunction vector. When we include symmetry each wfn might have different size
    size_t m_dimension; //!< Overall dimension of the direct product Fock space

    /*!
     * @brief Overall wavefunction.
     *
     * Psi = {psi_{p,q}_0, psi_{p,q}_1A, psi_{p,q}_2A,..., psi_{p,q}_1B, ..., psi_{p,q}_1A_1B, ...}
     * where psi_{p,q} is the electronic CI wavefunction and indices IA imply direct product with HO basis function I of
     * mode A.
     */
    std::vector<Wavefunction> m_wfn;
public:
    /*!
     * @brief Checks that the two wavefunctions are of the same electronic State and of the same dimension.
     */
    bool compatible(const MixedWavefunction &other) const;
    /*! @copydoc gci::LinearAlgebra::axpy
     */
    void axpy(double a, const LinearAlgebra::vector<double> &other) override;

    /*!
     * @brief
     * @param a
     * @param x
     */
    void axpy(double a, const std::shared_ptr<LinearAlgebra::vector<double> > &x) {
        axpy(a, *x);
    }

    /*!
     * @brief
     * @param a
     * @param other
     */
    void axpy(double a, const std::map<size_t, double> &other) override {
        throw std::logic_error("Cannot assume sparse vector is a map");
    }

    /*!
     * @copydoc LinearAlgebra::vector::select
     */
    std::tuple<std::vector<size_t>, std::vector<double>>
    select(const vector<double> &measure, const size_t maximumNumber = 1000,
           const double threshold = 0) const override {return {{0}, {0}};};

    /*!
     * @brief
     * @param a
     */
    void scal(double a) override;
    /*!
     * @brief
     * @param other
     * @return
     */
    double dot(const LinearAlgebra::vector<double> &other) const override;

    double dot(const std::shared_ptr<LinearAlgebra::vector<double> > other) const {
        return dot(*other);
    }

    double dot(const std::unique_ptr<LinearAlgebra::vector<double> > other) const {
        return dot(*other);
    }

    double dot(const std::map<size_t, double> &other) const override {
        throw std::logic_error("Cannot assume sparse vector is a map");
    }

    void zero() override;

    /*!
     * @brief
     * @param option
     * @return
     */
    MixedWavefunction *clone(int option = 0) const override {return new MixedWavefunction(*this, option);}

    MixedWavefunction &operator*=(const double &value); ///< multiply by a scalar
    MixedWavefunction &operator+=(const MixedWavefunction &other); ///< add another wavefunction
    MixedWavefunction &operator-=(const MixedWavefunction &other); ///< subtract another wavefunction
    MixedWavefunction &operator-=(double); ///< subtract a scalar from every element
    MixedWavefunction &operator-(); ///< unary minus
    MixedWavefunction &
    operator/=(const MixedWavefunction &other); ///< element-by-element division by another wavefunction

    double update(const Wavefunction &diagonalH,
                  double &eTruncated,
                  const double dEmax = (double) 0) { }; ///< form a perturbation-theory update, and return the predicted energy change. eTruncated is the energy change lost by truncation
    /*!
     * \brief addAbsPower Evaluate this[i] += factor * abs(c[I])^k * c[I]
     * \param c
     * \param k
     * \param factor
     * \return a pointer to this
     */
    MixedWavefunction &addAbsPower(const MixedWavefunction &c, double k = 0, double factor = 1);
    friend double
    operator*(const MixedWavefunction &w1, const MixedWavefunction &w2);///< inner product of two wavefunctions

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
    void divide(const LinearAlgebra::vector<double> *a,
                const LinearAlgebra::vector<double> *b,
                double shift = 0,
                bool append = false,
                bool negative = false);

    std::map<std::string, double> m_properties;

    void settilesize(int t = -1, int a = -1, int b = -1) { };
};

double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2);///< inner product of two wavefunctions
MixedWavefunction operator+(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< add two wavefunctions
MixedWavefunction operator-(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< subtract two wavefunctions
MixedWavefunction operator/(const MixedWavefunction &w1,
                            const MixedWavefunction &w2); ///< element-by-element division of two wavefunctions
MixedWavefunction operator*(const MixedWavefunction &w1, const double &value);///< multiply by a scalar
MixedWavefunction operator*(const double &value, const MixedWavefunction &w1);///< multiply by a scalar


}  // namespace gci
#endif //GCI_GCIMIXEDWAVEFUNCTION_H
