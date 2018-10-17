#ifndef GCI_GCIVIBSPACE_H
#define GCI_GCIVIBSPACE_H

/*!
 * @brief Parameters defining the vibrational space
 */
struct VibSpace {
    int nMode; //!< number of vibrational modes
    int nModal; //!< totalnumber of modals per mode (i.e. including the ground state)
    int modeCoupling; //!< level of mode-mode coupling

    bool operator==(const VibSpace &other) const {
        return (nMode == other.nMode) && (nModal == other.nModal) && (modeCoupling == other.modeCoupling);
    }
};

#endif //GCI_GCIVIBSPACE_H
