#ifndef GCIEXCITATIONSET_H
#define GCIEXCITATIONSET_H
#include "gciString.h"
#include <vector>

namespace gci {
/*!
 * \brief Holds a set of excitations from a single String.
 *
 * Not stored in the class, but in the context in which it is used, are
 * the number of annihilations/creations in the excitation (which affects
 * the interpretation of orbitalAddress), and the StringSet into which the Excitation
 * points (via stringIndex).
 */
class Excitation
{
public:
    /*!
     * \brief Constructor
     * \param StringIndex points to the destination of the excitation in a StringSet
     * \param Phase plus or minus one, giving the parity of the line-up permutation
     * \param OrbitalAddress An address representing the orbital(s) involved.
     */
    Excitation(size_t StringIndex, int Phase, size_t OrbitalAddress);
    /*!
     * \brief stringIndex points to the destination of the excitation in a StringSet
     */
    size_t stringIndex;
    /*!
     * \brief phase plus or minus one, giving the parity of the line-up permutation
     */
    int phase;
    /*!
     * \brief orbitalAddress An address representing the orbital(s) involved.
     */
    size_t orbitalAddress;
};


class StringSet; // forward declaration
/*!
 * \brief Container for a number of Excitation objects all arising from the same base String
 */
class ExcitationSet : public std::vector<Excitation>
{
public:
    /*!
     * \brief Construct the ExcitationSet containing all excitations
     * from a given String
     * with a specified number of annihilations and creations.
     * \param from The String to be excited.
     * \param to StringSet against which results will be indexed.
     * \param annihilations How many annihilations.
     * \param creations How many creations.
     */
    ExcitationSet(String &from, StringSet &to, int annihilations, int creations);
    /*!
     * \brief The String to which this set relates
     */
    String From;
    /*!
     * \brief The StringSet containing the excited String objects
     */
    StringSet* To;
    /*!
     * \brief Generate printable representation of the object
     * \return Printable representation of the object
     */
    std::string str(int verbosity=0) const;

};

std::ostream& operator<<(std::ostream& os, ExcitationSet const& obj);

}
#include "gciStringSet.h"

#endif // GCIEXCITATIONSET_H
