#ifndef GCIState_H
#define GCIState_H
#include "gci.h"
#include "FCIdump.h"
#include "gciOrbitalSpace.h"
#include "gciPrintable.h"
#include <string>

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
 * @param filename is the file containing the FCIDUMP. If present, load is called.
 */
    State(std::string filename="");
    /**
     * @brief
     *
     * @param dump points to an FCIdump object. If present, load is called.
     */
    State(FCIdump* dump);
    /*!
     * \brief Construct a State object linked to a OrbitalSpace
     * \param os The orbital space
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
    State(OrbitalSpace *os, int nelec=0, int symmetry=1, int ms2=0);
    /*!
     * \brief Construct a State object with data copied from another State
     * \param s State to copy
     */
    State(State* s);

    ~State();
    /*!
     \brief
    load number of electrons, spin from FCIDUMP file.
     \param filename is the file containing the FCIDUMP.
    */
    void load(std::string filename="FCIDUMP");
    /*!
     \brief
    load number of electrons, spin from FCIDUMP file.
     \param dump is an FCIdump object.
    */
    void load(FCIdump* dump);
    /*!
      \brief
       Pointer to orbital basis set, if any
      */
    OrbitalSpace* orbitalSpace;
     /*! \brief Number of electrons */
    unsigned int nelec;
    /*! \brief Twice the spin quantum number, ie multiplicity minus one */
    int ms2;
    /*! \brief Spatial symmetry of state */
    unsigned int symmetry;
    std::string str(int verbosity=0) const;


};
}

using namespace gci;

#endif // GCIState_H
