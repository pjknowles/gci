#ifndef GCISTRINGSET_H
#define GCISTRINGSET_H
#include "gciString.h"
#include "gciExcitationSet.h"
#include "gciPrintable.h"
#include <vector>
#include <map>

namespace gci {

/*!
 * \brief The StringSet class holds a set of String objects, possibly the complete set for a given number of objects and boxes
 */
class StringSet : public std::vector<String>, public Printable
{
public:
    /*!
     * \brief StringSet default constructor produces an empty object
     */
    StringSet();
    /*!
     * \brief StringSet constructor
     * \param prototype A String object that will be copied into this as a source of the number of objects and boxes, and symmetry information
     * \param all Whether or not to construct the complete set of String objects
     * \param sym If all, specify symmetry of String objects
     */
    StringSet(String prototype, bool all=true, int sym=-1);
    /*!
     * \brief PartialWeightArray holds the partial weight array for addressing the full set of String objects
     */
    std::vector<std::vector<int> > PartialWeightArray;
    /*!
     * \brief Map from the summed partial weights to the canonical index of a String in this set
     */
    std::map<long,long> addressMap;
    /*!
     * \brief Populate the StringSet with the complete set of Strings
     * \param sym Restrict to those String objects with this symmetry if not negative
     */
    void complete(int sym=-1);
    /*!
     * \brief The symmetry of the StringSet, or -1 if no definite symmetry
     */
    int symmetry;
    /*!
     * \brief Calculate addressMap
     */
    void calculateAddressMap();
    /*!
     * \brief calculate the address of String in the Stringset
     * \param s Pointer to the String
     * \return index in this, or -1 if not found
     */
    long offset(String &s);
    /*!
     * \brief Generate all excitations from this StringSet to StringSet to.
     * \param to StringSet against which results will be indexed.
     * \param annihilations How many annihilations.
     * \param creations How many creations.
     * \return The vector of ExcitationSet objects.
     */
    std::vector<ExcitationSet> allExcitations(StringSet &to, int annihilations, int creations);

    /*!
     * \brief occupationNumbers
     * \return a one-dimensional array representing a two-dimensional matrix
     * dimensioned number of Strings in this Stringset (running fastest)
     * by number of orbitals (running slowest)
     */
    std::vector<double> occupationNumbers();

    std::string str(int verbosity=0) const;
private:
    String proto;
    static long binomial_coefficient(unsigned long n, unsigned long k) ;
};

}
#endif // GCISTRINGSET_H
