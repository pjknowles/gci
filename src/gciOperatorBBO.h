#ifndef GCI_GCIOPERATORBBO_H
#define GCI_GCIOPERATORBBO_H

#include "gciOperator.h"
#include "gciOptions.h"

#include <ostream>
#include <vector>
#include <string>

namespace gci {

/*!
 * @brief Full molecular Hamiltonian operator incorporating electronic, vibrational and interaction terms.
 */
class OperatorBBO {
public:
    std::string m_description; //! Description of the operator
    Operator m_Hel; //! Electronic Hamiltonian
    std::vector<Operator> m_Hvib; //! Vibrational Hamiltonian
    std::vector<Operator> m_HintEl; //! Electronic component of the interaction Hamiltonian
    std::vector<Operator> m_HintVib; //! Vibrational component of the interaction Hamiltonian
    std::vector<int> m_symMode; //! symmetry of each mode
    int m_nMode; //! number of vibrational modes
    int m_nModal; //! number of HO basis functions per mode
    std::vector<int> m_vibOcc; //! Occupancy of vibrational modals (which modal is occupied, ground state = 0)
    std::vector<double> m_freq; //! vibrational frequencies in a.u.
    std::string m_fcidump; //! Root name of the fcidump files defining this Hamiltonian

    /*!
     * @brief Initialises Hel, Hvib and Hint from fcidump files. The vibrational Hamiltonian is assumed to be harmonic.
     * @param symMode Symmetry of each mode
     * @param nMode Total number of modes
     * @param nModal Number of HO basis functions per mode
     * @param fcidump Root name of the fcidump files
     */
//    explicit OperatorBBO(std::vector<unsigned int> &symMode, std::vector<int> &vibOcc, int nMode, int nModal,
//                         std::string &fcidump, std::string &description);
    explicit OperatorBBO(Options &options, std::string description = "BBO Molecular Hamiltonian");
    ~OperatorBBO()=default;
};

std::ostream &operator<<(std::ostream &os, const OperatorBBO &obj);

}// namespace gci
#endif //GCI_GCIOPERATORBBO_H
