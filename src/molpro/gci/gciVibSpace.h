#ifndef GCI_GCIVIBSPACE_H
#define GCI_GCIVIBSPACE_H

namespace molpro {
namespace gci {
/*!
 * @brief Parameters defining the vibrational space
 */
struct VibSpace {
  int nMode;  //!< Number of vibrational modes
  int nModal; //!< Total number of modals per mode (i.e. including the ground state)
  int excLvl; //!< Maximum number of modes simultaneously excited
              //    int excQuanta; //!< Total quanta of siumultaneous excitations among all modes

  VibSpace(int mode, int modal, int excLevel);

  /*!
   * @brief Returns true if the two spaces are the same
   */
  bool operator==(const VibSpace &other) const;
};

} // namespace gci
} // namespace molpro
#endif // GCI_GCIVIBSPACE_H
