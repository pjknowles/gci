#ifndef GCI_GCIMIXEDOPERATOR_H
#define GCI_GCIMIXEDOPERATOR_H

#include <vector>

#include "gciOperator.h"

/*!
 * @brief Mixed fermionic-bosonic Hamiltonian operator. Specialised to second-order expansion of the molecular
 * Hamiltonian in the BBO project.
 *
 * H = H'_el + H_vib + H_el(A) Q_A + H_s(A) dQ_A
 *
 * H`_el = H_el + H_s(A,A)
 *  -- first term is the electronic Hamiltonian at expansion point
 *  -- second term is the diagonal born oppenheimer correction <dQ_A Ph_I| dQ_A Ph_J>
 *
 * H_vib -- vibrational Harmonic oscillator Hamiltonian
 *
 * H_el(A) = dH_el/dQ_A
 *  -- derivative of the electronic Hamiltonian
 *
 *
 */
class MixedOperator {
    MixedOperator(const FCIdump &fcidump);

    int nMode; //!< Number of vibrational modes
    Operator H_el; //!< Electronic Hamiltonian
    Operator H_vib; //!< Vibrational Hamiltonian
    std::vector<Operator> Hel_A; //!< First order expansion of electronic Hamiltonian
    std::vector<Operator> Hs_A; //!< First order kinetic energy coupling term
    std::vector<Operator> Hs_AA; //!< Second order kinetic energy coupling term
};


#endif //GCI_GCIMIXEDOPERATOR_H
