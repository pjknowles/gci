#ifndef GCIWAVEFUNCTION_H
#define GCIWAVEFUNCTION_H
#include "gci.h"
#include "gciHamiltonian.h"
#include "gciState.h"
#include "gciStringSet.h"
#include "gciPrintable.h"

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
     * \brief Construct a Wavefunction object linked to a Hamiltonian
     * \param h The hamiltonian
     * \param nelec Number of electrons
     * \param symmetry Spatial symmetry
     * \param ms2 Sz quantum number times 2
     */
    Wavefunction(Hamiltonian *h, int nelec, int symmetry, int ms2);

    /*!
     * \brief Wavefunction copy constructor
     * \param other The object to be copied
     */
    Wavefunction( const Wavefunction& other);

    StringSet alphaStrings[8]; ///< The alpha-spin strings defining the CI basis
    StringSet betaStrings[8]; ///< The beta-spin strings defining the CI basis

    void allocate_buffer(); ///< allocate buffer to full size
    size_t size(); ///< the size of the space

    void diagonalHamiltonian(); ///< set this object to the diagonal elements of the hamiltonian

    std::string toString(int verbosity=0) const;

//    Wavefunction& operator=(const double &value);
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
