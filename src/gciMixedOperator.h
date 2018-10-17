#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>

#include "gciOperator.h"
#include "gciHProduct.h"

namespace gci {
/*!
 * @brief Mixed fermionic-bosonic Hamiltonian operator. Specialised to second-order expansion of the molecular
 * Hamiltonian in the BBO project.
 *
 * H = H'el + Hvib + Hel(A) Q_A + Hs(A) dQ_A
 *
 * H`el = Hel + Hs(A,A)
 *  -- first term is the electronic Hamiltonian at expansion point
 *  -- second term is the diagonal born oppenheimer correction <dQ_A Ph_I| dQ_A Ph_J>
 *
 * Hvib -- vibrational Harmonic oscillator Hamiltonian
 *
 * Hel(A) = dHel/dQ_A
 *  -- derivative of the electronic Hamiltonian
 *
 *
 */
class MixedOperator {
public:
    MixedOperator(const FCIdump &fcidump);
    ~MixedOperator() = default;

    /*!
     * @brief enumerators for vibrational operators
     *
     * HO -- Harmonic oscillator operator
     * The rest are self-explanatory
     */
    enum class VibOp {
        HO, Q, dQ, Qsq
    };

    /*!
     * @brief Returns the expectation value for one of the implemented operators
     */
    double expectVal(const HProduct &bra, const HProduct &ket, VibOp);


    int nMode; //!< Number of vibrational modes
    std::vector<double> freq; //!< Harmonic frequencies
    double zpe; //!< Vibrational zero point energy at HO level
    gci::Operator Hel; //!< Electronic Hamiltonian
    std::vector<gci::Operator> Hel_A; //!< First order expansion of electronic Hamiltonian
    std::vector<gci::Operator> Hs_A; //!< First order kinetic energy coupling term
    std::vector<gci::Operator> Hs_AA; //!< Second order kinetic energy coupling term
protected:
    /*!
     * @brief Expectation value of vibrational Hamiltonian <HO| Hvib | HO>
     * @param hamiltonian Coupled electronic-vibrational Hamiltonian
     * @param bra Vibrational mode
     * @param modal Harmonic Oscillator modal index
     */
    double O_Hvib(const HProduct &bra, const HProduct &ket) const;

    /*!
     * @brief Expectation value of Q_A
     * @param hamiltonian Coupled electronic-vibrational Hamiltonian
     * @param bra Vibrational Hartree product for the bra
     * @param ket Vibrational Hartree product for the ket
     * @return
     */
    double O_Q(const HProduct &bra, const HProduct &ket) const;

    /*!
     * @brief Expectation value of d/dQ_A
     * @param hamiltonian Coupled electronic-vibrational Hamiltonian
     * @param bra Vibrational Hartree product for the bra
     * @param ket Vibrational Hartree product for the ket
     * @return
     */
    double O_dQ(const HProduct &bra, const HProduct &ket) const;

    /*!
     * @brief Body for evaluation of operators Q and dQ
     * @tparam Func Function evaluating the operator, given the mode, highest modal and the difference (bra - ket)
     * @param bra Vibrational Hartree product for the bra
     * @param ket Vibrational Hartree product for the ket
     * @return
     */
    template<class Func>
    double QtypeOperator(const HProduct &bra, const HProduct &ket, const Func &func) const {
        if (bra.modeCouplingLvl() != ket.modeCouplingLvl())
            throw std::logic_error("Bra and Ket are of different mode-coupling level. Always 0.");
        double O_Q = 0.0;
        auto nExc = bra.modeCouplingLvl();
        for (int i = 0; i < nExc; ++i) {
            if (bra[i][0] != ket[i][0])
                throw std::logic_error("Bra and Ket excite different modes. Always 0.");
            auto mode = bra[i][0];
            auto iModal = bra[i][1];
            auto jModal = ket[i][1];
            auto diff = iModal - jModal;
            if (diff == 0) continue;
            if (std::abs(diff) != 1)
                throw std::logic_error("Bra and Ket are separated by more than 1 excitation. Always 0.");
            double n = diff > 0 ? iModal : jModal;
            O_Q += func(mode, n, diff);
        }
        return O_Q;

    }
};

}  // namespace gci

#endif //GCI_GCIMIXEDOPERATOR_H
