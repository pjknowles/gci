#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>
#include <Operator.h>
#include <FCIdump.h>

#include "gciHProduct.h"
#include "gciRun.h"

namespace gci {

/*!
 * @brief Enumerators for vibrational operators
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
    std::vector<int> mode; //!< Indices of the operator (e.g. for Q_A, mode = {A})
};

/*!
 * @brief One term in the mixed operator
 */
struct MixedOpTerm {
    MixedOpTerm(VibOp vibOp, const FCIdump &fcidump) : vibOp(std::move(vibOp)), Hel(constructOperator(fcidump)) { };
    MixedOpTerm(VibOp vibOp, const SymmetryMatrix::Operator &op) : vibOp(std::move(vibOp)), Hel(op) { };
    VibOp vibOp;
    SymmetryMatrix::Operator Hel;
};

/*!
 * @brief Mixed fermionic-bosonic Hamiltonian operator. Specialised to second-order expansion of the molecular
 * Hamiltonian in the BBO project.
 *
 * H = H'el + Hvib + Hel(A) Q_A + Ht(A) dQ_A
 *
 * H`el = Hel + Ht(A,A)
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
    /*!
     * @brief Constructs the full Hamiltonian from FCIdump files.
     *
     * In the following:
     *      - "fcidump" is the name of fcidump file passed to gci
     *      - modes are counted from 0, i.e. "{fcidump_t_0, fcidump_t_1, ...}"
     *
     * Naming converntion:
     *     "fcidump" -- Pure Electronic (Hel) and Vibrational (Hvib) Hamiltonian
     *     "fcidump"_d1_`A` -- first derivative of Hel wrt mode `A`
     *     "fcidump"_d2_`A`_`B` -- second derivative of Hel wrt modes `A` and `B`
     *                          -- @warning A >= B
     *     "fcidump"_t1_`A` -- first order kinetic energy coupling term with mode `A`
     *     "fcidump"_t2_`A` -- second order kinetic energy coupling term with mode `A`
     *
     * FCIdump parameters:
     *      int NMODE -- number of vibrational modes
     *      int INC_D1 -- include first derivative of Hel term, will search for "fcidump"_d1_`A` files
     *      int INC_D2 -- include second derivative of Hel term, will search for "fcidump"_d1_`AA` files
     *      int INC_T1 -- include first order kinetic energy coupling term, will search for "fcidump"_t1_`A` files
     *      int INC_T2 -- include second order kinetic energy coupling term, will search for "fcidump"_t2_`A` files
     *      (INC* statements are bool, 0 = false, 1 = true)
     * @warning If INC_* statement is found, but any of the fcidump files do not exist, than the corresponding
     *          term is assumed to be zero (fcidump's that ARE found are still included).
     *
     * @param fcidump Name of FCIdump file passed to gci
     */
    explicit MixedOperator(const FCIdump &fcidump);
    MixedOperator() : nMode(0), freq({}), zpe(0), Hel({0, 0}, 0, false, 0, true, "dummy"), Hmix({}), m_inc_d1(false),
                      m_inc_d2(false), m_inc_T1(false), m_inc_T2(false) { }
    ~MixedOperator() = default;

    /*!
     * @brief Returns the expectation value for one of the implemented operators
     * @result <bra| vibOp | ket>
     */
    double expectVal(const HProduct &bra, const HProduct &ket, const VibOp &vibOp) const;

    auto begin() const {return Hmix.begin();}
    auto end() const {return Hmix.end();}
    auto cbegin() const {return Hmix.cbegin();}
    auto cend() const {return Hmix.cend();}


    int nMode; //!< Number of vibrational modes
    std::vector<double> freq; //!< Harmonic frequencies
    double zpe; //!< Vibrational zero point energy at HO level
    SymmetryMatrix::Operator Hel; //!< Electronic Hamiltonian
    std::map<VibOpType, std::vector<MixedOpTerm>> Hmix; //!< All mixed vibrational-electronic terms, mapped by VibOpType
    bool inc_d1() {return m_inc_d1;}
    bool inc_d2() {return m_inc_d2;}
    bool inc_T1() {return m_inc_T1;}
    bool inc_T2() {return m_inc_T2;}
protected:
    bool m_inc_d1;
    bool m_inc_d2;
    bool m_inc_T1;
    bool m_inc_T2;

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
