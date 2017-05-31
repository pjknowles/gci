#ifndef GCIORBITALSPACE_H
#define GCIORBITALSPACE_H
#include "gciSymmetrySpace.h"
#include "FCIdump.h"
#include <map>

namespace gci {
/*!
 * \brief The OrbitalSpace class: a container for an orbital space as specified via symmetries (1..8) of orbitals
 */
class OrbitalSpace : public SymmetrySpace
{
public:
  /*!
     \brief construct OrbitalSpace object
     \param filename : if present, call load
     \param verbosity : how much to print
    */
  OrbitalSpace(std::string filename="", int verbosity=0);
  /*!
     \brief construct OrbitalSpace object
     \param dump : if present, call load
     \param verbosity : how much to print
    */
  OrbitalSpace(FCIdump* dump, int verbosity=0);
  virtual ~OrbitalSpace(){}

  void load(std::string filename="FCIDUMP", int verbosity=0); /**< \brief load from FCIDUMP */
  void load(FCIdump* dump, int verbosity=0); /**< \brief load from FCIDUMP */

  bool spinUnrestricted; /**< \brief whether alpha and beta spin orbitals are different */
  std::vector<unsigned int> orbital_symmetries;///< \brief spatial symmetry of orbitals (0-7)
  //    SymmetrySpace antisymmetricPairSpace; ///< \brief pairs of orbitals, ij equiv ji, antisymmetric
  std::map<int,SymmetrySpace> pairSpace; ///< \brief pairs of orbitals, ij not equiv ji
  //    SymmetrySpace symmetricPairSpace; ///< \brief pairs of orbitals, ij equiv ji
  //    size_t std::map<int,SymmetrySpace>::operator [](size_t index) const;

  /*!
     * \brief Construct a printable representation of the object
     * \param verbosity how much information to include
     * \return printable representation of the object
     */
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;

  /*!
     * \brief calculate canonical index of an orbital within its symmetry.
     * \param i Absolute number (starting with 1) of orbital.
     * \return Number of orbitals of the same symmetry canonically before i.
     */
  unsigned int orbitalIndex(unsigned int i) const;

  /*!
     * \brief pairIndex: addressing of pairs of orbitals given as original orbital numbers
     * \param i orbital index in original numbering scheme
     * \param j orbital index in original numbering scheme
     * \param parity 0 i not equiv j, otherwise i equiv j
     * \return
     */
  size_t pairIndex(unsigned int i, unsigned int j, int parity=0) const;

  /*!
     * \brief quadIndex: addressing of quads of orbitals given as original orbital numbers
     * \param i orbital index in original numbering scheme
     * \param j orbital index in original numbering scheme
     * \param k orbital index in original numbering scheme
     * \param l orbital index in original numbering scheme
     * \param parity 0 i not equiv j, otherwise i equiv j
     * \param parity2 0 ij not equiv kl, otherwise ij equiv kl
     * \return
     */
  size_t quadIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l, int parity=0, int parity2=0) const;

};
}

using namespace gci;

#endif // GCIORBITALSPACE_H
