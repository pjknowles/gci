#ifndef GCI_GCIMIXEDWAVEFUNCTION_H
#define GCI_GCIMIXEDWAVEFUNCTION_H

#include <vector>

#include <LinearAlgebra.h>
#include "gciWavefunction.h"

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

    ~MixedWavefunction() = default;

private:
    int nMode_; //!< number of vibrational modes
    int nModal_; //!< number of modals per vibrational mode
    int modeCoupling_; //!< level of mode-mode coupling = 1, or 2

    /*!
     * @brief Overall wavefunction.
     *
     * Psi = {psi_{p,q}_0, psi_{p,q}_1A, psi_{p,q}_2A,..., psi_{p,q}_1B, ..., psi_{p,q}_1A_1B, ...}
     * where psi_{p,q} is the electronic CI wavefunction and indices IA imply direct product with HO basis function I of
     * mode A.
     */
    std::vector<Wavefunction> wfn_;
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
     * @brief
     * @param a
     */
    void scal(double a) override;

    /*!
     * @brief
     * @param option
     * @return
     */
    MixedWavefunction *clone(int option = 0) const override {return new MixedWavefunction(*this, option);}

    double operator*(const MixedWavefunction &w1, const MixedWavefunction &w2);///< inner product of two wavefunctions
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
    void times(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b);
    void divide(const LinearAlgebra::vector<double> *a,
                const LinearAlgebra::vector<double> *b,
                double shift = 0,
                bool append = false,
                bool negative = false);
    void zero() override;
};


#endif //GCI_GCIMIXEDWAVEFUNCTION_H
