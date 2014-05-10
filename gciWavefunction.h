#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"
#include "gciState.h"
#include "gciStringSet.h"
#include "gciPrintable.h"
#include "gciDeterminant.h"

namespace gci {
/*!
 * \brief The Wavefunction class, which holds a configuration expansion in the tensor space defined by a hamiltonian
 */
class Wavefunction : public State
{
public:
    /**
 * @brief
 *
 * @param filename is the file containing the FCIDUMP. If present, load is called.
 */
    Wavefunction(std::string filename="");
    /**
     * @brief
     *
     * @param dump points to an FCIdump object. If present, load is called.
     */
    Wavefunction(FCIdump* dump);

    /*!
     * \brief Construct a Wavefunction object linked to an OrbitalSpace
     * \param h The hamiltonian
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
    Wavefunction(OrbitalSpace *h, int nelec, int symmetry, int ms2);

    /*!
     * \brief Wavefunction copy constructor
     * \param other The object to be copied
     */
    Wavefunction( const Wavefunction& other);

    StringSet alphaStrings[8]; ///< The alpha-spin strings defining the CI basis
    StringSet betaStrings[8]; ///< The beta-spin strings defining the CI basis

    void allocate_buffer(); ///< allocate buffer to full size
    size_t size(); ///< the size of the space

    void diagonalHamiltonian(Hamiltonian& hamiltonian); ///< set this object to the diagonal elements of the hamiltonian

    /*!
     * \brief find the index of the smallest component
     * \return offset in buffer
     */
    size_t minloc();
    /*!
     * \brief find the index of the largest component
     * \return offset in buffer
     */
    size_t maxloc();

    /*!
     * \brief Get a component of the wavefunction
     * \param offset which component to get
     * \return  the value of the component
     */
    double at(size_t offset);
    /*!
     * \brief get the determinant corresponding to a particular component of the wavefunction
     * \param offset which component to get
     * \return  the determinant
     */
    Determinant determinantAt(size_t offset);

    /*!
     * \brief Fill this object with the action of a hamiltonian operator on another wavefunction
     * \param h the hamiltonian
     * \param w the wavefunction
     */
    void hamiltonianOnWavefunction(Hamiltonian& h, Wavefunction& w);

    std::string str(int verbosity=0) const;

//    Wavefunction& operator=(const double &value);
    void set(size_t offset, const double val);///< set one element to a scalar
    void set(const double val);///< set all elements to a scalar
//    Wavefunction& operator=(const Wavefunction &other); ///< copy
    Wavefunction& operator*=(const double &value); ///< multiply by a scalar
    Wavefunction& operator+=(const Wavefunction &other); ///< add another wavefunction
    Wavefunction& operator-=(const Wavefunction &other); ///< subtract another wavefunction

    friend double operator*(const Wavefunction &w1, const Wavefunction &w2);///< inner product of two wavefunctions
private:
    void buildStrings(); ///< build alphaStrings and betaStrings
    size_t dimension; ///< the size of the space
    std::vector<double> buffer; ///< buffer to hold coefficients describing the object
    bool compatible(const Wavefunction &other) const; ///< whether this wavefunction is on the same space as another

};
    Wavefunction operator+(const Wavefunction &w1, const Wavefunction &w2); ///< add two wavefunctions
    Wavefunction operator-(const Wavefunction &w1, const Wavefunction &w2); ///< subtract two wavefunctions
    Wavefunction operator*(const Wavefunction &w1, const double &value);///< multiply by a scalar
    Wavefunction operator*(const double &value, const Wavefunction &w1);///< multiply by a scalar
}
using namespace gci;
#endif // GCIWAVEFUNCTION_H
