#include "gciString.h"
#include "gciStringSet.h"

#include <iostream>
#include <sstream>
#include <string.h>

// using String = gci::String;
namespace molpro {
namespace gci {

size_t String::StringNotFound = (size_t)-1; ///< conventional null value for index
size_t String::keyUnassigned = (size_t)-1;  ///< conventional null value for key

String::String(const State *State, const int Spin) : m_spin(Spin) {
  if (State != nullptr) {
    orbitalSpace.reset(new OrbitalSpace(*State->orbitalSpace));
  }
  nullify();
}

String::String(const std::vector<char> bytestream, const State *State) {
  if (State == nullptr) {
    orbitalSpace = nullptr;
  } else {
    orbitalSpace.reset(new OrbitalSpace(*State->orbitalSpace));
    //    cout << "assigned orbitalSpace="<<orbitalSpace<<std::endl;
    //    cout <<orbitalSpace->str()<<std::endl;
  }
  nullify();
  size_t offset = 0;
  memcpy(&nelec, &bytestream[offset], sizeof(nelec));
  offset += sizeof(nelec);
  memcpy(&m_spin, &bytestream[offset], sizeof(m_spin));
  offset += sizeof(m_spin);
  m_orbitals.resize((size_t)nelec);
  memcpy(&m_orbitals[0], &bytestream[offset], sizeof(orbital_type) * m_orbitals.size());
  offset += sizeof(orbital_type) * m_orbitals.size();
  //  cout <<"constructed String from bytestream nelec="<<nelec<<", spin="<<spin<<", orbitals=";
  //  for (size_t i=0; i<nelec; i++) cout <<orbitals()[i]<<","; cout <<std::endl;
  //  cout <<str(1)<<std::endl;
  //  cout <<"back from str"<<std::endl;
}

std::vector<char> String::serialise() const {
  std::vector<char> bytestream(sizeof(nelec) + sizeof(m_spin) + nelec * sizeof(orbital_type));
  size_t offset = 0;
  memcpy(&bytestream[offset], &nelec, sizeof(nelec));
  offset += sizeof(nelec);
  memcpy(&bytestream[offset], &m_spin, sizeof(m_spin));
  offset += sizeof(m_spin);
  memcpy(&bytestream[offset], &m_orbitals[0], sizeof(orbital_type) * m_orbitals.size());
  offset += sizeof(orbital_type) * m_orbitals.size();
  //  cout <<"serialised String to bytestream nelec="<<nelec<<", spin="<<spin<<", orbitals=";
  //  for (size_t i=0; i<nelec; i++) cout <<orbitals()[i]<<","; cout <<std::endl;
  //  cout <<str()<<std::endl;
  return bytestream;
}

// int String::create(unsigned int orbital) {
////      cout << "String::create before="<<str()<<", orbital="<<orbital<<std::endl;
//  //        cout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
//  //        cout << "hamiltonian "<<(hamiltonian!=nullptr)<<std::endl;
//  if (orbitalSpace==nullptr)
//    throw std::logic_error("String::create missing orbitalSpace");
//  //        cout << "basisSize "<<hamiltonian->total()<<std::endl;
//  if (orbitalSpace==nullptr || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total()) throw
//  std::range_error("invalid orbital");
//  //    cout <<"make iterator "<<std::endl;
////  std::vector<unsigned int>::iterator ilast=orbitals_.begin();
//  //    cout <<"iterator OK"<<std::endl;

//  int phase=((m_orbitals.size()/2)*2 == m_orbitals.size()) ? 1 : -1;
//  //    cout <<"phase="<<phase<<std::endl;
//  //    cout <<"spin="<<spin<<std::endl;
//  for (std::vector<orbital_type>::iterator i = m_orbitals.begin(); i!=m_orbitals.end(); ++i) {
//    if (*i==orbital) return 0; // exclusion principle
//    if (*i > orbital){
//      ms2+=m_spin;
//      nelec++;
//      symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
////                  cout <<"create orbital="<<*i <<" with symmetry="<<orbitalSpace->orbital_symmetries[*i-1]<<",
/// giving total symmetry"<<symmetry<<std::endl;
//      m_orbitals.insert(i,orbital);
////                  cout << "String::create inserts, after="<<str()<<", phase="<<phase<<std::endl;
//      return phase;
//    }
//    phase=-phase;
//  }
//  ms2+=m_spin;
//  nelec++;
//  symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
//  m_orbitals.insert(m_orbitals.end(),orbital);
////      cout << "String::create final append, after="<<str()<<", phase="<<phase<<std::endl;
//  return phase;
//}

// int String::destroy(unsigned int orbital) {
//  if (orbitalSpace==nullptr || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total() ) throw
//  std::range_error("invalid orbital"); if (m_orbitals.size() <= 0) return (int) 0; //throw "too few electrons in
//  String";
//  //    cout << "String::destroy before="<<str()<<", orbital="<<orbital<<std::endl;
////  int phase=1;
//  int phase=((m_orbitals.size()/2)*2 == m_orbitals.size()) ? -1 : 1;
//  for (std::vector<orbital_type>::iterator i = m_orbitals.begin(); i!=m_orbitals.end(); ++i) {
//    if (*i==orbital)  {
//      ms2-=m_spin;
//      nelec--;
//      symmetry^=orbitalSpace->orbital_symmetries[*i-1];
//      m_orbitals.erase(i);
//      //            cout << "String::destroy succeeds, after="<<str()<<", phase="<<phase<<std::endl;
//      return phase;
//    }
//    phase=-phase;
//  }
//  //    cout << "String::destroy fails, after="<<str()<<std::endl;
//  return (int)0; // exclusion principle
//}

void String::nullify() {
  m_orbitals.clear();
  ms2 = 0;
  nelec = 0;
  symmetry = 0;
  m_key = keyUnassigned;
}

const std::vector<String::orbital_type> &String::orbitals() const { return m_orbitals; }

std::string String::str(int verbosity, unsigned int columns) const {
  std::string result;
  //      cout <<"String::str orbitals_[0]" <<orbitals_[0]<<std::endl;
  //  cout << "String::str length of orbitals_="<<orbitals_.size()<<std::endl;
  if (verbosity >= 0) {
    for (std::vector<orbital_type>::const_iterator i = m_orbitals.begin(); i != m_orbitals.end(); ++i) {
      if (i != m_orbitals.begin())
        result.append(",");
      std::stringstream ss;
      int ispin = (int)(*i) * (int)m_spin;
      ss << ispin;
      std::string rr;
      ss >> rr;
      result.append(rr);
    }
    //      cout << "end of loop occ orbital"<<std::endl;
    std::stringstream ss;
    ss << " [";
    if (m_key != keyUnassigned)
      ss << m_key << ".";
    ss << computed_symmetry() + 1 << "]"; // internally symmetries are implemented 0-7, but externally as 1-8
    std::string rr;
    ss >> rr;
    result.append(rr);
  }
  return result;
}

unsigned int String::computed_symmetry(bool nocheck) const {
  unsigned int s = 0;
  //  cout << "in computed_symmetry"<<std::endl;
  //  cout <<orbitalSpace->str()<<std::endl;
  for (int i = 0; i < (int)m_orbitals.size(); i++) {
    s ^= orbitalSpace->orbital_symmetries[m_orbitals[i] - 1];
    //        cout <<"orbital "<<orbitals_[i]<<",  symmetry="<<hamiltonian->orbital_symmetries[orbitals_[i]-1]<<" total
    //        symmetry now "<<s<<std::endl;
  }
  if (s != symmetry && !nocheck) {
    cout << "s=" << s << ", symmetry=" << symmetry << std::endl;
    throw std::runtime_error("String symmetry messed up");
  }
  //  cout << "computed_symmetry"<<s<<std::endl;
  return s;
}

bool String::next(int sym) {
  if (sym < 0) { // generate the next string of any symmetry
    orbital_type limit = orbitalSpace->total();
    std::vector<orbital_type>::iterator k;
    orbital_type floor = 0;
    for (std::vector<orbital_type>::reverse_iterator i = m_orbitals.rbegin(); i != m_orbitals.rend(); ++i) {
      floor = ++(*i);
      k = i.base();
      if (*i <= limit)
        break;
      limit--;
    }
    if (limit <= orbitalSpace->total() - m_orbitals.size())
      return false; // we ran out of boxes to put the objects into
    while (k != m_orbitals.rbegin().base())
      *(k++) = ++floor;
    symmetry = computed_symmetry(true);
    //    cout << "String::next returns with symmetry="<<symmetry<<std::endl;
    return true;
  } else { // call myself until we get the symmetry required
    bool notexhausted;
    while ((notexhausted = next()) && (symmetry = computed_symmetry()) != (unsigned int)sym)
      ;
    return notexhausted;
  }
}

bool String::first(int n, int sym) {
  //    cout << "String::first " << n <<" nelec="<<nelec<< std::endl;
  if (n <= 0)
    n = nelec;
  if (n <= 0)
    n = m_orbitals.size();
  nullify();
  //    cout << n <<std::endl;
  for (unsigned int i = 1; i <= (unsigned int)n; i++)
    create(i);
  nelec = n;
  //    cout << "String::first first go, sym, symmetry "<<sym<<symmetry<<std::endl;
  if (sym < 0 || computed_symmetry() == (unsigned int)sym)
    return true;
  return next(sym);
}

String String::exhausted;

size_t String::index(const StringSet &set) const {
  const auto i = set.addressMap.find(m_key);
  return (i == set.addressMap.end()) ? StringNotFound : i->second;
}

bool String::operator==(const String &other) const { return m_key == other.m_key; }
} // namespace gci
} // namespace molpro
