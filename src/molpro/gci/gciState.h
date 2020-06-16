#ifndef GCIState_H
#define GCIState_H
#include <memory>
#include <string>

#include "molpro/gci/gci.h"
#include "molpro/gci/gciOptions.h"
#include "molpro/gci/gciOrbitalSpace.h"
#include "molpro/gci/gciPrintable.h"

namespace gci {
/**
 * @brief
 * Class represents a quantum-mechanical state with or without reference to a particular representation
 *
 */
class State : public Printable {
public:
  /**
   * @brief
   *
   */
  State() : nelec(0), ms2(0), symmetry(0){};
  /**
   * @brief
   *
   * @param dump points to an Options object. If present, load is called.
   */
  explicit State(const Options &dump);
  /*!
   * \brief Construct a State object linked to a OrbitalSpace
   * \param os The orbital space
   * \param nelec Number of electrons
   * \param symmetry Spatial symmetry
   * \param ms2 Sz quantum number times 2
   */
  explicit State(OrbitalSpace *os, int nelec = 0, int symmetry = 1, int ms2 = 0);
  explicit State(OrbitalSpace &os, int nelec = 0, int symmetry = 1, int ms2 = 0);
  explicit State(const molpro::Operator &o, int n, int s, int m2)
      : orbitalSpace(new OrbitalSpace(o)), nelec(n), ms2(m2), symmetry(s) {}

  virtual ~State() = default;
  /*!
     \brief
    load number of electrons, spin
     \param dump is an Options object.
    */
  void load(const Options &dump);
  /*!
      \brief
       Pointer to orbital basis set, if any
      */
  std::shared_ptr<OrbitalSpace> orbitalSpace;
  /*! \brief Number of electrons */
  unsigned int nelec;
  /*! \brief Twice the spin quantum number, ie multiplicity minus one */
  int ms2;
  /*! \brief Spatial symmetry of state */
  unsigned int symmetry;
  std::string str(int verbosity = 0, unsigned int columns = UINT_MAX) const override;
};
} // namespace gci

#endif // GCIState_H
