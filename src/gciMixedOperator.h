#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>

#include "gciOperator.h"
#include "gciHProduct.h"

namespace gci {

/*!
 * @brief enumerators for vibrational operators
 *
 * HO -- Harmonic oscillator operator
 * The rest are self-explanatory
 */
enum class VibOpType {
    HO, Q, dQ, Qsq
};

/*!
 * @brief Specifies one of ther hardcoded vibrational operators
 */
class VibOp {
public:
    VibOp(const VibOp &other) : type(other.type), mode(other.mode) {
        bool pass = false;
        switch (type) {
            case VibOpType::HO: if (mode.empty()) pass = true;
            case VibOpType::Q: if (mode.size() == 1) pass = true;
            case VibOpType::dQ: if (mode.size() == 1) pass = true;
            case VibOpType::Qsq: if (mode.size() <= 2) pass = true;
//            default: throw std::logic_error("Operator not implemented");
        }
        if (!pass) throw std::logic_error("Inconsistent dimensionality of the operator");
    }

    const VibOpType type; //!< Operator type
    const std::vector<int> mode; //!< indices of the operator (i.e for Q_A, mode = A)
};

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
     * @brief Returns the expectation value for one of the implemented operators
     * @result <bra| vibOp | ket>
     */
    double expectVal(const HProduct &bra, const HProduct &ket, const VibOp &vibOp);


    int nMode; //!< Number of vibrational modes
    std::vector<double> freq; //!< Harmonic frequencies
    double zpe; //!< Vibrational zero point energy at HO level
    gci::Operator Hel; //!< Electronic Hamiltonian
    std::vector<gci::Operator> Hel_A; //!< First order expansion of electronic Hamiltonian
    std::vector<gci::Operator> Hel_AB; //!< First order expansion of electronic Hamiltonian
    std::vector<gci::Operator> Hs_A; //!< First order kinetic energy coupling term
    std::vector<gci::Operator> Hs_AA; //!< Second order kinetic energy coupling term
protected:
    /*!
     * @brief Expectation value of vibrational Hamiltonian <HO| Hvib | HO>
     */
    double O_Hvib(const HProduct &bra, const HProduct &ket) const;

    /*!
     * @brief Expectation value of Q_A
     * @return
     */
    double O_Q(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const;

    /*!
     * @brief Expectation value of d/dQ_A
     * @return
     */
    double O_dQ(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const;

    /*!
     * @brief Expectation value of Q_A*Q_B
     */
    double O_Qsq(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const;

    /*!
     * @brief Body for evaluation of operators Q and dQ
     * @tparam Func Function evaluating the operator, given the mode, highest modal and the difference (bra - ket)
     * @param bra Vibrational Hartree product for the bra
     * @param ket Vibrational Hartree product for the ket
     * @return
     */
    double QtypeOperator(const HProduct &bra, const HProduct &ket, const std::function<double(double, int, int)> &func,
                         int targetMode) const;
};

}  // namespace gci

#endif //GCI_GCIMIXEDOPERATOR_H
