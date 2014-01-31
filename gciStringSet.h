#ifndef GCISTRINGSET_H
#define GCISTRINGSET_H
#include "gciString.h"
#include <vector>

namespace gci {

/*!
 * \brief The StringSet class holds a set of String objects, possibly the complete set for a given number of objects and boxes
 */
class StringSet : public std::vector<String>
{
public:
    /*!
     * \brief StringSet default constructor produces an empty object
     */
    StringSet();
    /*!
     * \brief StringSet constructor
     * \param prototype A String object that will be copied into this as a source of the number of objects and boxes, and symmetry information
     * \param all Whether or not to construct the complete set of Strings
     */
    StringSet(String prototype, bool all=true);
    /*!
     * \brief PartialWeightArray holds the partial weight array for addressing the full set of String objects
     */
    std::vector<std::vector<int> > PartialWeightArray;
    /*!
     * \brief Populate the StringSet with the complete set of Strings
     */
    void complete();
private:
    String proto;
    long binomial_coefficient(unsigned long n, unsigned long k) ;
};

}
#endif // GCISTRINGSET_H
