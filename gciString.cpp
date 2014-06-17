#include "gciString.h"
#include "gciStringSet.h"

#include <iostream>
#include <sstream>

size_t gci::String::StringNotFound=ULLONG_MAX; ///< conventional null value for index
size_t gci::String::keyUnassigned=ULLONG_MAX; ///< conventional null value for key

String::String(State* State, int Spin)
{
  if (State == NULL) {
    orbitalSpace=NULL;
  } else {
    orbitalSpace=State->orbitalSpace;
  }
  nullify();
  spin=Spin;
}

int String::create(unsigned int orbital) {
//      xout << "String::create before="<<str()<<", orbital="<<orbital<<std::endl;
  //        xout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
  //        xout << "hamiltonian "<<(hamiltonian!=NULL)<<std::endl;
  if (orbitalSpace==NULL)
    throw "String::create missing orbitalSpace";
  //        xout << "basisSize "<<hamiltonian->total()<<std::endl;
  if (orbitalSpace==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total()) throw "invalid orbital";
  //    xout <<"make iterator "<<std::endl;
//  std::vector<unsigned int>::iterator ilast=orbitals_.begin();
  //    xout <<"iterator OK"<<std::endl;

  int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
  //    xout <<"phase="<<phase<<std::endl;
  //    xout <<"spin="<<spin<<std::endl;
  for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
    if (*i==orbital) return 0; // exclusion principle
    if (*i > orbital){
      ms2+=spin;
      nelec++;
      symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
//                  xout <<"create orbital="<<*i <<" with symmetry="<<orbitalSpace->orbital_symmetries[*i-1]<<", giving total symmetry"<<symmetry<<std::endl;
      orbitals_.insert(i,orbital);
//                  xout << "String::create inserts, after="<<str()<<", phase="<<phase<<std::endl;
      return phase;
    }
    phase=-phase;
  }
  ms2+=spin;
  nelec++;
  symmetry^=orbitalSpace->orbital_symmetries[orbital-1];
  orbitals_.insert(orbitals_.end(),orbital);
//      xout << "String::create final append, after="<<str()<<", phase="<<phase<<std::endl;
  return phase;
}

int String::destroy(unsigned int orbital) {
  if (orbitalSpace==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) orbitalSpace->total() ) throw "invalid orbital";
  if (orbitals_.size() <= 0) return (int) 0; //throw "too few electrons in String";
  //    xout << "String::destroy before="<<str()<<", orbital="<<orbital<<std::endl;
//  int phase=1;
  int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? -1 : 1;
  for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
    if (*i==orbital)  {
      ms2-=spin;
      nelec--;
      symmetry^=orbitalSpace->orbital_symmetries[*i-1];
      orbitals_.erase(i);
      //            xout << "String::destroy succeeds, after="<<str()<<", phase="<<phase<<std::endl;
      return phase;
    }
    phase=-phase;
  }
  //    xout << "String::destroy fails, after="<<str()<<std::endl;
  return (int)0; // exclusion principle
}

void String::nullify()
{
  orbitals_.clear();
  ms2=0;
  nelec=0;
  symmetry=0;
  key=keyUnassigned;
}

std::vector<unsigned int> String::orbitals() const {
  return orbitals_;
}

std::string String::str(int verbosity) const {
  std::string result;
  //    xout <<"String::str orbitals_[0]" <<orbitals_[0]<<std::endl;
  if (verbosity >=0) {
    for (std::vector<unsigned int>::const_iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
      if (i!=orbitals_.begin()) result.append(",");
      std::stringstream ss;
      int ispin=(int)(*i)*(int)spin;
      ss << ispin;
      std::string rr;
      ss >> rr;
      result.append(rr);
    }
    std::stringstream ss;
    ss << " [";
    if (key != keyUnassigned) ss << key << ".";
    ss << computed_symmetry()+1 <<"]"; // internally symmetries are implemented 0-7, but externally as 1-8
    std::string rr;
    ss >> rr;
    result.append(rr)
        ;    }
  return result;
}

unsigned int String::computed_symmetry(bool nocheck) const
{
  unsigned int s=0;
  for (int i=0; i<(int)orbitals_.size(); i++) {
    s ^= orbitalSpace->orbital_symmetries[orbitals_[i]-1];
    //        xout <<"orbital "<<orbitals_[i]<<",  symmetry="<<hamiltonian->orbital_symmetries[orbitals_[i]-1]<<" total symmetry now "<<s<<std::endl;
  }
  if (s!=symmetry && !nocheck) {
    xout << "s="<<s<<", symmetry="<<symmetry<<std::endl;
    throw "String symmetry messed up";
  }
  return s;
}

bool String::next(int sym) {
  if (sym<0)
  { // generate the next string of any symmetry
    unsigned int limit=orbitalSpace->total();
    std::vector<unsigned int>::iterator k;
    unsigned int floor=0;
    for (std::vector<unsigned int>::reverse_iterator i = orbitals_.rbegin(); i!=orbitals_.rend(); ++i) {
      floor=++(*i);
      k=i.base();
      if (*i <= limit) break;
      limit--;
    }
    if (limit <= orbitalSpace->total()-orbitals_.size()) return false; // we ran out of boxes to put the objects into
    while (k!=orbitals_.rbegin().base())
      *(k++)=++floor;
    symmetry=computed_symmetry(true);
    //    xout << "String::next returns with symmetry="<<symmetry<<std::endl;
    return true;
  }
  else { // call myself until we get the symmetry required
    bool notexhausted;
    while ((notexhausted=next())&&(symmetry=computed_symmetry())!=(unsigned int)sym) ;
    return notexhausted;
  }
}

bool String::first(int n, int sym) {
  //    xout << "String::first " << n <<" nelec="<<nelec<< std::endl;
  if (n <=0 ) n=nelec;
  if (n <=0 ) n=orbitals_.size();
  nullify();
  //    xout << n <<std::endl;
  for (unsigned int i=1;i<=(unsigned int)n;i++)
    create(i);
  nelec=n;
  //    xout << "String::first first go, sym, symmetry "<<sym<<symmetry<<std::endl;
  if (sym<0 || computed_symmetry() == (unsigned int) sym) return true;
  return next(sym);
}

String String::exhausted;


size_t String::index(const StringSet& set) const
{
  std::map<size_t,size_t>::const_iterator i;
  i = set.addressMap.find(key);
  return (i == set.addressMap.end()) ? StringNotFound : i->second;
}

bool String::operator ==(const String& other) const
{
  return key == other.key;
}
