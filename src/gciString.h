#ifndef GCISTRING_H
#define GCISTRING_H
#include "gci.h"

#include "gciState.h"
#include <vector>
#include <climits>
#include <cstdint>

namespace gci {
class StringSet;
/*!
 \brief
A string, which is an ordered set of orbitals
*/
class String : public State
{
  friend class ExcitationSet;
  friend class StringSet;
public:
  /*!
   * \brief String Construct a vacuum object
   * \param State Some State object from which to copy number of electrons etc for bound checking, and to define orbital symmetries
   * \param Spin 1=alpha, -1=beta
*/
  String(const State *State=nullptr, const int Spin=1);
  /*!
   * \brief String Construct from a serialised representation of data
   * \param bytestream representation produced by previous invocation of serialise()
   * \param State Some State object from which to copy number of electrons etc for bound checking, and to define orbital symmetries
   */
  String(const std::vector<char> bytestream, const State* State=nullptr);
  /*!
     \brief

     \param orbital Add an orbital to the string.
     \return int On exit, the phase change required to bring the determinant into canonical form is returned (plus or minus 1), or else zero if the orbital was already present in the determinant.
    */
  int create(unsigned int orbital)
  {
  //      xout << "String::create before="<<str()<<", orbital="<<orbital<<std::endl;
    //        xout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
    //        xout << "hamiltonian "<<(hamiltonian!=nullptr)<<std::endl;
//    if (orbitalSpace==nullptr)
//      throw std::logic_error("String::create missing orbitalSpace");
    //        xout << "basisSize "<<hamiltonian->total()<<std::endl;
//    if (orbitalSpace==nullptr || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total()) throw std::range_error("invalid orbital");
    assert(orbitalSpace!=nullptr);
    assert(orbital != 0);
    assert(orbital <= orbitalSpace->total());
    //    xout <<"make iterator "<<std::endl;
  //  std::vector<unsigned int>::iterator ilast=orbitals_.begin();
    //    xout <<"iterator OK"<<std::endl;

//    int phase=((m_orbitals.size()/2)*2 == m_orbitals.size()) ? 1 : -1;
    int phase=1-2*(m_orbitals.size()%2);
    //    xout <<"phase="<<phase<<std::endl;
    //    xout <<"spin="<<spin<<std::endl;
    for (std::vector<orbital_type>::const_iterator i = m_orbitals.begin(); i!=m_orbitals.end(); ++i) {
      if (*i==orbital) return 0; // exclusion principle
      if (*i > orbital){
        ms2+=m_spin;
        nelec++;
        symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
  //                  xout <<"create orbital="<<*i <<" with symmetry="<<orbitalSpace->orbital_symmetries[*i-1]<<", giving total symmetry"<<symmetry<<std::endl;
        m_orbitals.insert(i,orbital);
  //                  xout << "String::create inserts, after="<<str()<<", phase="<<phase<<std::endl;
        return phase;
      }
      phase=-phase;
    }
    ms2+=m_spin;
    nelec++;
    symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
    m_orbitals.insert(m_orbitals.end(),orbital);
  //      xout << "String::create final append, after="<<str()<<", phase="<<phase<<std::endl;
    return phase;
  }

  /*!
     \brief

     \param orbital Remove an orbital from the string.
     \return int the phase change required to bring the determinant into canonical form before annihilation is returned (plus or minus 1), or else zero if the orbital was not present in the determinant.
    */
  int destroy(unsigned int orbital)
  {
    if (orbitalSpace==nullptr || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total() ) throw std::range_error("invalid orbital");
    if (m_orbitals.size() <= 0) return (int) 0; //throw "too few electrons in String";
    //    xout << "String::destroy before="<<str()<<", orbital="<<orbital<<std::endl;
    int phase=(m_orbitals.size()%2) ? 1 : -1;
    for (std::vector<orbital_type>::const_iterator i = m_orbitals.begin(); i!=m_orbitals.end(); ++i) {
      if (*i==orbital)  {
        ms2-=m_spin;
        nelec--;
        symmetry^=orbitalSpace->orbital_symmetries[*i-1];
        m_orbitals.erase(i);
        //            xout << "String::destroy succeeds, after="<<str()<<", phase="<<phase<<std::endl;
        return phase;
      }
      phase=-phase;
    }
    //    xout << "String::destroy fails, after="<<str()<<std::endl;
    return (int)0; // exclusion principle
  }
  /*!
     \brief advance to the canonically next string
    \param sym the desired symmetry of the string. Negative signifies any will do.
     \return whether successful; false if you try to advance the canonically last string
    */
  bool next(int sym=-1);

  /*!
     * \brief Set the string to the vacuum
     */
  void nullify();    /*!
     \brief
    Set to the first string with n objects
    \param n the number of objects. If 0, use whatever we have presently.
    \param sym the desired symmetry of the string. Negative signifies any will do.
    \return false if it was not possible to make even one string, otherwise true
    */
  bool first(int n=0, int sym=-1);
  /*!
   * \brief Holds orbital labels. 8-bit unsigned is fine for up to 511 orbitals
   */
  typedef uint8_t orbital_type;
  /*!
   * \brief Holds hash keys
   */
  typedef size_t key_type;
  const std::vector<orbital_type> &orbitals() const;  /*!< The orbitals that make up the string */
  /*!
   * \brief Hash key that can be associated with this object
   * \return
   */
  const key_type& key() {return m_key;}
  template<class T>
  /*!
   * \brief Generate the hash key that can be associated with this object by summing a partial weight array.
   * \param partialWeightArray The partial weight array.
   */
  void keygen(const std::vector<std::vector<T> > & partialWeightArray)
  {
    m_key=0;
    for (size_t k=0; k<m_orbitals.size(); k++)
      m_key += partialWeightArray[k][m_orbitals[k]-1];
  }

  std::string str(int verbosity=0, unsigned int columns=UINT_MAX) const;
  /*!
     * \brief Calculate the spatial symmetry
     * \param nocheck If false, check whether the result is equal to the maintained symmetry variable
     * \return
     */
  unsigned int computed_symmetry(bool nocheck=false) const;
  static String exhausted; /*!< returned by next() when we're already on the last string */


  static key_type keyUnassigned; ///< conventional null value for key
  static size_t StringNotFound; ///< conventional null value for index

  /*!
     * \brief Find the location of this String in a given StringSet
     * \param set the StringSet that hopefully contains this String
     * \return the offset in set or StringNotFound if not in set
     */
  size_t index(const StringSet& set) const;

  /*!
   * \brief operator == test whether two String objects are identical
   * \param other
   * \return
   */
  bool operator==(const String& other) const;

  /*!
   * \brief serialise Produce a byte representation of the object's data
   * \return
   */
  std::vector<char> serialise() const;
private:
  key_type m_key; ///< \brief Hash key that can be associated with this object
  std::vector<orbital_type> m_orbitals; /*!< The orbitals that make up the string */
  int_least8_t m_spin; ///< \brief spin 1=alpha, -1=beta
};
}

using namespace gci;

#endif // GCISTRING_H
