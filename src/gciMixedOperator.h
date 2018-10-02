#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>

#include "gciOperator.h"

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
    ~MixedOperator()=default;

    int nMode; //!< Number of vibrational modes
    std::vector<double> freq; //!< Harmonic frequencies
    double zpe; //!< Vibrational zero point energy at HO level
    gci::Operator Hel; //!< Electronic Hamiltonian
    std::vector<gci::Operator> Hel_A; //!< First order expansion of electronic Hamiltonian
    std::vector<gci::Operator> Hs_A; //!< First order kinetic energy coupling term
    std::vector<gci::Operator> Hs_AA; //!< Second order kinetic energy coupling term
};


#endif //GCI_GCIMIXEDOPERATOR_H
