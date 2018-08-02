#include "gciOperatorBBO.h"
#include "utility"

#include "FCIdump.h"

namespace gci {
//OperatorBBO::OperatorBBO(std::vector<int> &symMode, std::vector<int> &vibOcc, int nMode, int nModal,
//                         std::string &fcidump, std::string &description) :
//        m_Hel(Operator::construct(FCIdump(fcidump))), m_vibOcc(vibOcc), m_symMode(symMode), m_nMode(nMode),
//        m_nModal(nModal), m_description(description) {
OperatorBBO::OperatorBBO(Options &options, std::string description) :
        m_description(std::move(description)),
        m_nMode(options.parameter("NMODE", m_nMode)),
        m_nModal(options.parameter("NMODAL", m_nModal)),
        m_vibOcc(options.parameter("VIBOCC", m_vibOcc)),
        m_symMode(options.parameter("SYMMODE", m_symMode)),
        m_freq(options.parameter("FREQ", m_freq)),
        m_fcidump(options.parameter("FCIDUMP", std::vector<std::string>(1, "")).at(0)),
        m_Hel(Operator::construct(FCIdump(m_fcidump))) {
    m_Hel.m_description = "Electronic Hamiltonian at the reference geometry";
    // initiate interaction and vibrational Hamiltonian
    for (int iMode = 0; iMode < m_nMode; ++iMode) {
        dim_t dimension{(unsigned int) m_nModal};
        std::vector<int> orbital_symmetries(m_nModal, m_symMode[iMode]);
        std::string descr;
        descr = "Vibrational Hamiltonian: mode = " + std::to_string(iMode);
        Operator hVib(dimension, orbital_symmetries, 1, false, (unsigned int) m_symMode[iMode], true, true, descr);
        descr = "Vibrational component of Interaction Hamiltonian: mode = " + std::to_string(iMode);
        Operator hIntVib(dimension, orbital_symmetries, 1, false, (unsigned int) m_symMode[iMode], true, true, descr);
        hVib.O1(true).assign(0.0);
        hIntVib.O1(true).assign(0.0);
        for (int jModal = 0; jModal < m_nModal; ++jModal) {
            hVib.element(jModal, m_symMode[iMode], jModal, m_symMode[iMode]) = m_freq[iMode] * (jModal + 0.5);
            hIntVib.element(jModal, m_symMode[iMode], jModal + 1, m_symMode[iMode]) = std::sqrt((double) jModal + 1);
        }
        std::string fcidumpN = m_fcidump + std::to_string(iMode);
        Operator hIntEl(Operator::construct(FCIdump(fcidumpN)));
        m_Hvib.push_back(hVib);
        m_HintEl.push_back(hIntEl);
        m_HintVib.push_back(hIntVib);
    }
}


std::ostream &operator<<(std::ostream &os, const OperatorBBO &obj) {
    os << "OperatorBBO: " << obj.m_description << std::endl;
    os << "nMode = " << obj.m_nMode << std::endl;
    os << "nModal = " << obj.m_nModal << std::endl;
    os << "vibOcc = ";
    for (auto iter = obj.m_vibOcc.begin(); iter < obj.m_vibOcc.end(); ++iter) os << *iter << " ";
    os << std::endl << "freq = ";
    for (auto iter = obj.m_freq.begin(); iter < obj.m_freq.end(); ++iter) os << *iter << " ";
    os << std::endl;
    os << obj.m_Hel << std::endl;
    for (int iMode = 0; iMode < obj.m_nMode; ++iMode) {
        os << obj.m_Hvib[iMode] << std::endl;
        os << obj.m_HintEl[iMode] << std::endl;
        os << obj.m_HintVib[iMode] << std::endl;
    }
    return os;
}

} // namespace gci
