#ifndef GCI_GCIVIBSPACE_H
#define GCI_GCIVIBSPACE_H

namespace gci {
/*!
 * @brief Parameters defining the vibrational space
 */
struct VibSpace {
    int nMode; //!< number of vibrational modes
    int nModal; //!< totalnumber of modals per mode (i.e. including the ground state)
    int modeCoupling; //!< level of mode-mode coupling

    VibSpace(int mode, int modal, int modeCouplingLvl);

    /*!
     * @brief Returns true if the two spaces are the same
     */
    bool operator==(const VibSpace &other) const;
};

}  // namespace gci
#endif //GCI_GCIVIBSPACE_H
