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
     * \param sym If all, specify symmetry of String objects; -1 denotes all symmetries
     */
  StringSet(String prototype, bool all=true, int sym=-1);
  /*!
     * \brief Construct a StringSet from another one by applying annihilations and creations
     * \param referenceSpace
     * \param annihilations number of annihilations to make
     * \param creations number of creations to make
     * \param sym specify symmetry of String objects; -1 denotes all symmetries
     * \param parallel whether to use parallel algorithm
     */
  StringSet(const StringSet &referenceSpace, int annihilations, int creations, int sym=-1, bool parallel=false);
  /*!
     * \brief Construct a StringSet from a vector of other ones by applying annihilations and creations
     * \param referenceSpaces vector of reference spaces
     * \param annihilations number of annihilations to make
     * \param creations number of creations to make
     * \param sym specify symmetry of String objects; -1 denotes all symmetries
     * \param parallel whether to use parallel algorithm
     */
  StringSet(const std::vector<StringSet> &referenceSpaces, int annihilations, int creations, int sym=-1,bool parallel=false);
  /*!
     * \brief PartialWeightArray holds the partial weight array for addressing the full set of String objects
     */
  std::vector<std::vector<int> > PartialWeightArray;
  /*!
     * \brief Map from the summed partial weights to the canonical index of a String in this set
     */
  std::map<size_t,size_t> addressMap;
  /*!
     * \brief Populate the StringSet with the complete set of String objects
     * \param sym Restrict to those String objects with this symmetry if not negative
     */
  void complete(int sym=-1);
  /*!
     * \brief The symmetry of the StringSet, or -1 if no definite symmetry
     */
  int symmetry;
  /*!
     * \brief Append to a StringSet from another one by applying annihilations and creations
     * \param referenceSpace
     * \param annihilations
     * \param creations
     * \param sym specify symmetry of String objects; -1 denotes all symmetries
     * \param parallel whether to run in parallel
     */
  void addByOperators(const StringSet &referenceSpace, int annihilations, int creations, int sym=-1, bool parallel=false);
  /*!
     * \brief Append to a StringSet from a collection of other ones by applying annihilations and creations
     * \param referenceSpaces
     * \param annihilations
     * \param creations
     * \param sym specify symmetry of String objects; -1 denotes all symmetries
     * \param parallel whether to use parallel algorithm
     */
  void addByOperators(const std::vector<StringSet> &referenceSpaces, int annihilations, int creations, int sym=-1,bool parallel=false);
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
     * dimensioned number of String objects in this StringSet (running fastest)
     * by number of orbitals (running slowest)
     */
  std::vector<double> occupationNumbers();
  /*!
     * \brief Add a new element to the end of the set
     * \param s The element to add
     */
  void insert(String &s);

  /*!
     * \brief Find the locations of the String objects of this in a given other StringSet
     * \param set the StringSet that hopefully contains this String
     * \return the offsets in set or StringNotFound if not in set
     */
  std::vector<size_t> index(const StringSet& set) const;

  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;
  String proto; ///< prototype String
private:
  void setupPartialWeightArray();
  static long binomial_coefficient(unsigned long n, unsigned long k) ;
  void makekey(String &s);
};

}
#endif // GCISTRINGSET_H
