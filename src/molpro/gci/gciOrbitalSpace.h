#ifndef GCIORBITALSPACE_H
#define GCIORBITALSPACE_H
#include "gciOptions.h"
#include "gciSymmetrySpace.h"
#include <map>
#include <molpro/symmetry_matrix/Operator.h>

namespace molpro {
namespace gci {
/*!
 * \brief The OrbitalSpace class: a container for an orbital space as specified via symmetries (1..8) of orbitals
 */
class OrbitalSpace : public SymmetrySpace {
public:
  /*!
     \brief construct OrbitalSpace object
     \param dump : if present, call load
     \param verbosity : how much to print
    */
  OrbitalSpace(const Options &dump, int verbosity = 0) { load(dump); }

  OrbitalSpace(const std::vector<int> syms, bool uhf = false, int verbosity = 0) { load(syms, uhf, verbosity); }

  OrbitalSpace(const molpro::Operator &o, int verbosity = 0) {
    std::vector<int> syms;
    for (unsigned int sym = 0; sym < 8; sym++)
      for (size_t i = 0; i < o.dimension(sym, 0, 0); i++)
        syms.push_back(sym + 1);
    load(syms, o.m_uhf, verbosity);
  }

  virtual ~OrbitalSpace() {}

  void load(const Options &options, int verbosity = 0);

  void load(const std::vector<int> syms, bool uhf = false, int verbosity = 0);

  bool spinUnrestricted;                        /**< \brief whether alpha and beta spin orbitals are different */
  std::vector<unsigned int> orbital_symmetries; ///< \brief spatial symmetry of orbitals (0-7)
  //    SymmetrySpace antisymmetricPairSpace; ///< \brief pairs of orbitals, ij equiv ji, antisymmetric
  std::map<int, SymmetrySpace> pairSpace; ///< \brief pairs of orbitals, ij not equiv ji
  //    SymmetrySpace symmetricPairSpace; ///< \brief pairs of orbitals, ij equiv ji
  //    size_t std::map<int,SymmetrySpace>::operator [](size_t index) const;

  /*!
   * \brief Construct a printable representation of the object
   * \param verbosity how much information to include
   * \param columns Page width
   * \return printable representation of the object
   */
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const;

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
  size_t pairIndex(unsigned int i, unsigned int j, int parity = 0) const;

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
  size_t quadIndex(unsigned int i, unsigned int j, unsigned int k, unsigned int l, int parity = 0,
                   int parity2 = 0) const;
};
} // namespace gci
} // namespace molpro

#endif // GCIORBITALSPACE_H
