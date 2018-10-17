#ifndef GCI_GCIMIXEDWAVEFUNCTION_H
#define GCI_GCIMIXEDWAVEFUNCTION_H

#include <vector>

#include <LinearAlgebra.h>
#include "gciWavefunction.h"
#include "gciMixedOperator.h"
#include "gciHProductSet.h"

namespace gci {


/*!
 * @brief Mixed bosonic and fermionic wavefunction represented as configuration expansion in the direct product
 * space of Slater Determenants and Hartree Product of modals. It implements storage of the wavefunction
 * coefficients and application of the Hamiltonitan ( r = H c).
 *
 * Vibrational basis is specified by mode coupling level, CIS, CISD etc, and maximum excitation level.
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
 *
 * H2_int = (d/dX_A d/dX_B H_el) dX_A dX_B
 * couples electronic degrees of freedom with two vibrational modes
 *
 * This introduces a lot of sparsity in Hamiltonian, but the wavefunction remains coupled.
 *
 * Wavefunction with only 1 mode coupling:
 * Psi = sum_A C_{p,q,I_A} |p,q> |I_A> |0_B> |0_C> ...
 *
 * Order of vibrational product basis:
 *  - index(CIS) < index(CISD)
 *  - size(CIS, mode_A) = nModal_A
 *  - Let Phi_mode_A be hartree products with mode A excited
 *    - CIS = {Phi_mode_1, Phi_mode_2, ...}
 *    - CISD = {Phi_mode_12, Phi_mode_13, ..., Phi_mode_23, ...}
 *  - and continued for CISDT, CISDTQ etc.
 *
 * In practice the coefficients C_{p,q,I_A} are stored as a vector<C_{p,q}>(n = nModes * nModals).
 * This class is in fact mostly a wrapper of gci::Wavefunction with a few hardcoded expressions for the vibrational
 * parts of the Hamiltonian.
 */
class MixedWavefunction : public LinearAlgebra::vector<double> {
    // Inhereting the vector class brings the functions that need to be implemented for linera algebra solver to work.
public:
    /*!
     * @brief Constructs the mixed wavefunction.
     * @param state Electronic state
     * @param vibSate Vibrational state
     */
    MixedWavefunction(const State &state, const VibSpace &vibSpace);
//    MixedWavefunction(const MixedWavefunction &other, int option) {*this = other;}
    ~MixedWavefunction() override = default;

    size_t size() const {return m_dimension;}
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

    /*!
     * @brief set this object to the diagonal elements of the hamiltonian
     * @param hamiltonian
     */
    void diagonalOperator(const MixedOperator &hamiltonian);

protected:

    VibSpace m_vibSpace; //!< Parameters defining the vibrational space of current wavefunction
    HProductSet m_vibBasis; //!< Vibrational basis for the full space of current wavefunction
    size_t m_elDim; //!< Dimension of the electronic (slater determinant) space
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
    /*! @copydoc LinearAlgebra::vector::axpy
     */
    void axpy(double a, const LinearAlgebra::vector<double> &other) override;
    /*!
     * @copydoc LinearAlgebra::vector::axpy(scalar,const vector<scalar>&)
     */
    void axpy(double a, const std::shared_ptr<LinearAlgebra::vector<double> > &other) {
        axpy(a, *other);
    }
    /*!
     * @copydoc LinearAlgebra::vector::axpy(scalar,const std::map<size_t,scalar>& )
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
     * @copydoc LinearAlgebra::vector::scal
     */
    void scal(double a) override;

    /*!
     * @copydoc LinearAlgebra::vector::dot
     */
    double dot(const LinearAlgebra::vector<double> &other) const override;

    /*!
     * @overload
     */
    double dot(const std::shared_ptr<LinearAlgebra::vector<double> > &other) const {
        return dot(*other);
    }

    /*!
     * @overload
     */
    double dot(const std::unique_ptr<LinearAlgebra::vector<double> > &other) const {
        return dot(*other);
    }

    /*!
     * @overload
     */
    double dot(const std::map<size_t, double> &other) const override {
        throw std::logic_error("Cannot assume sparse vector is a map");
    }

    /*!
     * @copydoc LinearAlgebra::vector::zero
     */
    void zero() override;

    /*!
     * @copydoc LinearAlgebra::vector::zero
     */
    MixedWavefunction *clone(int option = 0) const override {return new MixedWavefunction(*this, option);}

    void set(double val);///< set all elements to a scalar
    MixedWavefunction &operator*=(const double &value); //!< multiply by a scalar
    MixedWavefunction &operator+=(const MixedWavefunction &other); //!< add another wavefunction
    MixedWavefunction &operator-=(const MixedWavefunction &other); //!< subtract another wavefunction
    MixedWavefunction &operator+=(double); //!< add a scalar to every element
    MixedWavefunction &operator-=(double); //!< subtract a scalar from every element
    MixedWavefunction &operator-(); //!< unary minus
    MixedWavefunction &
    operator/=(const MixedWavefunction &other); //!< element-by-element division by another wavefunction

    /*!
     * @brief form a perturbation-theory update, and return the predicted energy change.
     * @param diagonalH Diagonal of the Hamiltonian matrix
     * @param eTruncated Energy change lost by truncation
     * @param dEmax
     * @return
     */
    double update(const Wavefunction &diagonalH,
                  double &eTruncated,
                  const double dEmax = 0.0) { };

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
};  // class MixedWavefunction

double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2);///< inner product of two wavefunctions
MixedWavefunction operator+(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< add two wavefunctions
MixedWavefunction operator-(const MixedWavefunction &w1, const MixedWavefunction &w2); ///< subtract two wavefunctions
MixedWavefunction operator/(const MixedWavefunction &w1,
                            const MixedWavefunction &w2); ///< element-by-element division of two wavefunctions
MixedWavefunction operator*(const MixedWavefunction &w1, const double &value);///< multiply by a scalar
MixedWavefunction operator*(const double &value, const MixedWavefunction &w1);///< multiply by a scalar


}  // namespace gci
#endif //GCI_GCIMIXEDWAVEFUNCTION_H
