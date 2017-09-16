#ifndef GCIState_H
#define GCIState_H
#include "gci.h"
#include "gciOrbitalSpace.h"
#include "gciPrintable.h"
#include "gciOptions.h"
#include <string>
#include <memory>

namespace gci {
/**
 * @brief
 * Class represents a quantum-mechanical state with or without reference to a particular representation
 *
 */
class State : public Printable
{
public:

  /**
 * @brief
 *
 */
  State() {};
  /**
     * @brief
     *
     * @param dump points to an Options object. If present, load is called.
     */
  State(const Options &dump);
  /*!
     * \brief Construct a State object linked to a OrbitalSpace
     * \param os The orbital space
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
  State(OrbitalSpace *os, int nelec=0, int symmetry=1, int ms2=0);
  State(OrbitalSpace &os, int nelec=0, int symmetry=1, int ms2=0);

  State(const State& source)
    : orbitalSpace(new OrbitalSpace(*source.orbitalSpace))
    , nelec(source.nelec), ms2(source.ms2), symmetry(source.symmetry)
  {
  }

  State& operator=(const State& source)
  {
    orbitalSpace.reset(new OrbitalSpace(*source.orbitalSpace));
    nelec = source.nelec;
    ms2 = source.ms2;
    symmetry = source.symmetry;
    return *this;
  }

  ~State();
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
  std::unique_ptr<OrbitalSpace> orbitalSpace;
  /*! \brief Number of electrons */
  unsigned int nelec;
  /*! \brief Twice the spin quantum number, ie multiplicity minus one */
  int ms2;
  /*! \brief Spatial symmetry of state */
  unsigned int symmetry;
  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;


};
}

using namespace gci;

#endif // GCIState_H
