#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>
#include <Operator.h>
#include <FCIdump.h>

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
 * @brief Specifies one of the hardcoded vibrational operators
 */
class VibOp {
public:
    VibOp(const VibOpType &type, const std::vector<int> &mode = {}) : type(type), mode(mode) {
        switch (type) {
            case VibOpType::HO:
                if (!mode.empty())
                    throw std::logic_error("HO is a full operator, no modes need to be specified");
                break;
            case VibOpType::Q: if (mode.size() != 1) throw std::logic_error("Q operator is one-dimensional");
                break;
            case VibOpType::dQ: if (mode.size() != 1) throw std::logic_error("dQ operator is one-dimensional");
                break;
            case VibOpType::Qsq: if (mode.size() != 2) throw std::logic_error("Qsq operator is two-dimensional");
                break;
            default: throw std::logic_error("Operator not implemented");
        }
    }

    VibOpType type; //!< Operator type
    std::vector<int> mode; //!< indices of the operator (i.e for Q_A, mode = A)
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
    double expectVal(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const;


    int nMode; //!< Number of vibrational modes
    std::vector<double> freq; //!< Harmonic frequencies
    double zpe; //!< Vibrational zero point energy at HO level
    SymmetryMatrix::Operator Hel; //!< Electronic Hamiltonian
    std::vector<SymmetryMatrix::Operator> Hel_A; //!< First order expansion of electronic Hamiltonian
    std::vector<SymmetryMatrix::Operator> Hel_AB; //!< First order expansion of electronic Hamiltonian
    std::vector<SymmetryMatrix::Operator> Hs_A; //!< First order kinetic energy coupling term
    std::vector<SymmetryMatrix::Operator> Hs_AA; //!< Second order kinetic energy coupling term
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
