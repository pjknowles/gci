#include "gciStringSet.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>

using StringSet = gci::StringSet;

StringSet::StringSet() : memory::vector<String>() {
//      xout <<parallel_rank<<"StringSet default constructor"<<std::endl;
}

StringSet::StringSet(const String &prototype, bool all, int sym) : memory::vector<String>() {
//      xout <<parallel_rank<<"StringSet prototype constructor "<<all<<", sym="<<sym<<std::endl;
  // copy prototype
  proto = prototype;
//  xout << "prototype="<<proto<<std::endl;
  symmetry = sym;
  setupPartialWeightArray();
  if (all) complete(sym);
}

StringSet::StringSet(const StringSet &referenceSpace, int annihilations, int creations, int sym, bool parallel)
    : memory::vector<String>() {
  //std::cout<<"crashing"<<std::endl;std::cout.flush();MPI_Abort(MPI_COMM_COMPUTE,12345);
  addByOperators(referenceSpace, annihilations, creations, sym, parallel);
//  xout << "StringSet constructor from referenceSpace size()="<<size()<<", rank="<<parallel_rank<<", parallel="<<parallel<<std::endl;
}

StringSet::StringSet(const std::vector<StringSet> &referenceSpaces,
                     int annihilations,
                     int creations,
                     int sym,
                     bool parallel)
    : memory::vector<String>() {
//  xout << "StringSet constructor from referenceSpaces referenceSpaces:" << std::endl;
//  for (const auto &ss: referenceSpaces) for (const auto &s: ss) xout << "referenceSpace String " << s << std::endl;
  addByOperators(referenceSpaces, annihilations, creations, sym, parallel);
//  xout << "StringSet constructor from referenceSpaces size()=" << size() << ", rank=" << parallel_rank << ", parallel="
//       << parallel << std::endl;
}

void StringSet::addByOperators(const std::vector<StringSet> &referenceSpaces,
                               int annihilations,
                               int creations,
                               int sym,
                               bool parallel) {
  auto p = profiler->push("StringSet::addByOperators[]");
  //std::cout <<"referenceSpaces="<<&referenceSpaces<<std::endl;std::cout.flush();
//  std::cout << "addByOperators: referenceSpaces.size()=" << referenceSpaces.size() << std::endl;
//  std::cout.flush();
  if (parallel) {
    size_t ntask = 0;
    for (const auto &referenceSpace : referenceSpaces) {
      ntask += referenceSpace.size();
    }
    //  xout << "parallel_rank="<<parallel_rank<<", ntask="<<ntask<<std::endl;
    DivideTasks(ntask);
  }
//  xout << "about to call addByOperators " << annihilations << creations << sym << std::endl;
  for (const auto &referenceSpace : referenceSpaces) {
//    xout << "Stringset " << *referenceSpace << std::endl;
    addByOperators(referenceSpace, annihilations, creations, sym, parallel);
  }
//  xout << "size="<<size()<<std::endl;
#ifdef HAVE_MPI_H
  if (parallel) {
    EndTasks();
    auto pp = profiler->push("StringSet::addByOperators:distribute");
    std::vector<char> serialised;
    for (const auto &s : *this) {
      std::vector<char> serialised1 = s.serialise();
      for (std::vector<char>::const_iterator c = serialised1.begin(); c != serialised1.end(); c++)
        serialised.push_back(*c);
    }
    //  xout << parallel_rank<<"serialised "<<serialised.size()<<" bytes from "<<size()<<" String objects="<<std::endl;
    //  xout << parallel_rank<<"addressMap.size()"<< addressMap.size()<<std::endl;
    //  for (size_t i=0; i<parallel_size; i++){
    //  MPI_Barrier(MPI_COMM_COMPUTE);
    //  if (i!=parallel_rank)continue;
    //  xout << "Barrier crossed parallel_size="<<parallel_size<<", parallel_rank="<<parallel_rank<<std::endl;
    //  xout << "end of addByOperators creation phase, size()="<<size()<<std::endl;
    //  xout << "unmerged addressMap rank="<<parallel_rank<<", size="<<addressMap.size()<<std::endl;
    //  for (std::map<size_t,size_t> ::const_iterator a=addressMap.begin(); a!=addressMap.end(); a++)
    //    xout <<parallel_rank<<": "<<(*a).first<<","<<(*a).second<<std::endl;
    //  }
    // aggregate on master process
    //    String refString = this->at(0);
    size_t initial_symmetry;
    for (initial_symmetry = 0; initial_symmetry < 8; initial_symmetry++)
      if (!referenceSpaces[initial_symmetry].empty()) break;
    if (initial_symmetry > 7) {
      xout << "Something has gone wrong with discovering reference space in gci::StringSet::addByOperators"
           << std::endl;
      exit(1);
    }
    String refString = referenceSpaces[initial_symmetry][0];
    int bytestreamsize;
    //  std::cout << "parallel_rank = "<<parallel_rank<<std::endl;std::cout.flush();
    if (parallel_rank > 0) {
      //    std::cout << "slave "<<std::endl;    std::cout.flush();
      auto len = (int) size();
      MPI_Send(&len, (int) 1, MPI_INT, 0, 0, MPI_COMM_COMPUTE);
      //    xout << "slave sends len="<<len<<std::endl;
      if (len > 0) {
        bytestreamsize = (int) serialised.size() / size();
        MPI_Send(&bytestreamsize, (int) 1, MPI_INT, 0, 1, MPI_COMM_COMPUTE);
        MPI_Send(&serialised[0], len * bytestreamsize, MPI_BYTE, 0, 2, MPI_COMM_COMPUTE);
      }
      MPI_Bcast(&len, (int) 1, MPI_INT, 0, MPI_COMM_COMPUTE);
      MPI_Bcast(&bytestreamsize, (int) 1, MPI_INT, 0, MPI_COMM_COMPUTE);
      //        std::cout << "slave after receving broadcast len "<<len<<", bytestreamsize="<<bytestreamsize<<std::endl; std::cout.flush();
      serialised.resize(static_cast<unsigned long>(len * bytestreamsize));
      clear();
      addressMap.clear();
      //    xout << "slave ready to bcast"<<std::endl;
      MPI_Bcast(&serialised[0], len * bytestreamsize, MPI_BYTE, 0, MPI_COMM_COMPUTE);
      for (size_t k = 0; k < (size_t) len; k++) {
        std::vector<char> s(static_cast<unsigned long>(bytestreamsize));
        memcpy(&s[0], &serialised[k * bytestreamsize], static_cast<size_t>(bytestreamsize));
        //            std::cout << "slave construct string"<<std::endl;
        String ss(s, &refString);
        //            std::cout << "slave insert string "<<ss.str()<<std::endl; std::cout.flush();
        insert(ss);
        //        std::cout << "slave inserted string"<<std::endl; std::cout.flush();
      }
    } else {
      //    xout <<"master size()="<<size()<<", bytestreamsize="<<serialised.size()/(size() ? size() : 1) <<", serialised.size()"<<serialised.size()<<std::endl;
      for (int iproc = 1; iproc < parallel_size; iproc++) {
        int len;
        MPI_Status status;
        MPI_Recv(&len, (int) 1, MPI_INT, iproc, 0, MPI_COMM_COMPUTE, &status);
        //    xout << "master receives len="<<len<<std::endl;
        if (len > 0) {
          MPI_Recv(&bytestreamsize, (int) 1, MPI_INT, iproc, 1, MPI_COMM_COMPUTE, &status);
          //      xout <<"received len="<<len<<", bytestreamsize="<<bytestreamsize<<std::endl;
          serialised.resize((size_t) len * bytestreamsize);
          MPI_Recv(&serialised[0], len * bytestreamsize, MPI_BYTE, iproc, 2, MPI_COMM_COMPUTE, &status);
          for (size_t k = 0; k < (size_t) len; k++) {
            std::vector<char> s(static_cast<unsigned long>(bytestreamsize));
            memcpy(&s[0], &serialised[k * bytestreamsize], static_cast<size_t>(bytestreamsize));
            //        xout << "master construct string"<<std::endl;
            String ss(s, &refString);
            //        xout <<"master before insert"<<std::endl;
            insert(ss);
            //        xout <<"master after insert"<<std::endl;
          }
          serialised.clear();
          for (const auto &s : *this) {
            std::vector<char> serialised1 = s.serialise();
            for (std::vector<char>::const_iterator c = serialised1.begin(); c != serialised1.end(); c++)
              serialised.push_back(*c);
          }
        }
      }
      auto len = (int) size();
      //        xout <<"master after serialising global list"<<std::endl;
      //    xout << "master ready to bcast len="<<len<<", bytestreamsize="<<bytestreamsize<<std::endl;
      bytestreamsize = serialised.size() / std::max((int) size(), 1);
      MPI_Bcast(&len, (int) 1, MPI_INT, 0, MPI_COMM_COMPUTE);
      MPI_Bcast(&bytestreamsize, (int) 1, MPI_INT, 0, MPI_COMM_COMPUTE);
      //xout <<"master after broadcasting len"<<std::endl;
      //    xout << "master ready to bcast"<<std::endl;
      MPI_Bcast(&serialised[0], len * bytestreamsize, MPI_BYTE, 0, MPI_COMM_COMPUTE);
      //        xout <<"master after broadcasting global list"<<std::endl;
    }
    //    xout << "Reached end of forked code rank="<<parallel_rank<<std::endl;
    //  for (size_t i=0; i<parallel_size; i++){
    //  MPI_Barrier(MPI_COMM_COMPUTE);
    //  if (i!=parallel_rank)continue;
    //  xout << "Barrier crossed parallel_size="<<parallel_size<<", parallel_rank="<<parallel_rank<<std::endl;
    //  xout << "end of addByOperators size()="<<size()<<std::endl;
    //  xout << "merged addressMap rank="<<parallel_rank<<", size="<<addressMap.size()<<std::endl;
    //  for (std::map<size_t,size_t> ::const_iterator a=addressMap.begin(); a!=addressMap.end(); a++)
    //    xout <<parallel_rank<<": "<<(*a).first<<","<<(*a).second<<std::endl;
    //  }
  }
#endif
}

void StringSet::addByOperators(const StringSet &referenceSpace,
                               int annihilations,
                               int creations,
                               int sym,
                               bool parallel) {
  size_t count = 0;
  size_t countall = 0;
  bool first = empty();
  symmetry = sym;
//  xout << "in StringSet::addByOperators, referenceSpace=" << referenceSpace.str(5) << std::endl;
//  xout << "referenceSpace.symmetry=" << referenceSpace.symmetry << ", nelec=" << referenceSpace.proto.nelec
//       << ", annihilations" << annihilations << ", creations" << creations << std::endl;
  if ((int) referenceSpace.proto.nelec + creations - annihilations < 0
      || (int) referenceSpace.proto.nelec + creations - annihilations > (int) referenceSpace.proto.orbitals().size())
    return; // null space because not enough electrons or holes left
//  auto p = profiler->push("addByOperators");
//  xout << "still here" << std::endl;
  int symexc =
      (referenceSpace.symmetry >= 0 && sym >= 0) ? referenceSpace.symmetry ^ sym : -1; // use symmetry if we can
  for (auto from : referenceSpace) {
//    xout << "next reference" << std::endl;
    countall++;
    if (parallel)
      if (!NextTask()) continue;
//    if (parallel && countall%parallel_size != parallel_rank) continue;
    count++;
//    xout << "from=" << from.str(5) << std::endl;
    if (annihilations + creations == 1) {
      for (int i = 0; i < (int) from.orbitalSpace->orbital_symmetries.size(); i++) {
        if (from.orbitalSpace->orbital_symmetries[i] == (unsigned int) symexc || symexc == -1) {
          String tt = from;
//                profiler->start("destroy or create");
          int phase = (annihilations > 0) ? tt.destroy(i + 1) : tt.create(i + 1);
//                profiler->stop("destroy or create");
          if (phase) {
            if (first) {
              proto = tt;
              setupPartialWeightArray();
              insert(tt);
              proto = tt;
              first = false;
            } else {
//                profiler->push("insert");
              insert(tt);
            }
          }
        }
      }
    } else if (annihilations + creations == 2) {
//     xout << "annihilations="<<annihilations<<", creations="<<creations<<", first="<<first<<std::endl;
      for (int j = 0; j < (int) from.orbitalSpace->orbital_symmetries.size(); j++) {
        String a = from;
        int phasea = (annihilations > 0) ? a.destroy(j + 1) : a.create(j + 1);
        if (phasea) {
          for (int i = 0; i < (int) from.orbitalSpace->orbital_symmetries.size(); i++) {
            if (from.orbitalSpace->orbital_symmetries[i]
                == ((unsigned int) symexc ^ from.orbitalSpace->orbital_symmetries[j]) || symexc == -1) {
              String tt = a;
              int phase = (annihilations > 1) ? phasea * tt.destroy(i + 1) : phasea * tt.create(i + 1);
              if (phase) {
                if (first) {
                  proto = tt;
                  setupPartialWeightArray();
                  insert(tt);
                  proto = tt;
                  first = false;
                } else
                  insert(tt);
              }
            }
          }
        }
      }
    }
  }
//  xout << "addByOperators parallel="<<parallel<<", parallel_rank="<<parallel_rank<<", count="<<count<<", countall="<<countall<<std::endl;
}

void StringSet::setupPartialWeightArray() { // set up partial weight array for addressing binomial distributions
//  auto p = profiler->push("StringSet::setupPartialWeightArray");
  int nitem = proto.nelec;
  int nbox = proto.orbitalSpace->total();
  PartialWeightArray = std::vector<std::vector<int> >(static_cast<unsigned long>(nitem),
                                                      std::vector<int>(static_cast<unsigned long>(nbox)));
  for (int k = 0; k < nitem; k++) {
    for (int l = 0; l < nbox; l++)
      PartialWeightArray[k][l] = 0;
    for (int l = k; l < nbox - nitem + k; l++)
      PartialWeightArray[k][l + 1] = binomial_coefficient(static_cast<unsigned long>(nbox - l - 1),
                                                          static_cast<unsigned long>(nitem - k - 1))
          + PartialWeightArray[k][l];
  }
  for (int k = 0; k < nitem - 1; k++) {
    for (int l = k; l < nbox - nitem + k + 1; l++)
      PartialWeightArray[k][l] = PartialWeightArray[k][l] - PartialWeightArray[k + 1][l + 1];
  }
  if (nitem > 0)
    for (int l = nitem - 1; l < nbox; l++)
      PartialWeightArray[nitem - 1][l] = l - nitem + 1;

  //    xout << "PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<PartialWeightArray[k][l]; xout <<std::endl; }
}

long StringSet::binomial_coefficient(unsigned long n, unsigned long k) {
  unsigned long i;
  long b;
  if (0 == k || n == k) {
    return 1;
  }
  if (k > n) {
    return 0;
  }
  if (k > (n - k)) {
    k = n - k;
  }
  if (1 == k) {
    return n;
  }
  b = 1;
  for (i = 1; i <= k; ++i) {
    b *= (n - (k - i));
    if (b < 0) return -1; /* Overflow */
    b /= i;
  }
  return b;
}

void StringSet::complete(int sym) {
  auto p = profiler->push("StringSet::complete");
//      xout <<"StringSet::complete prototype"<<proto.str(1)<<std::endl;
  String string(&proto);
  this->erase(this->begin(), this->end());
  if (string.first(proto.nelec, sym)) {
//        xout <<"StringSet::complete symmetry="<<sym<<" first String: "<<string<<std::endl;
    do {
//              xout << "in StringSet::complete about to push_back " << string.str(1) <<std::endl;
      this->insert(string);
    } while (string.next(sym));
  }
//      xout << "in StringSet::complete final list: " <<std::endl ;
//      for (iterator s=this->begin(); s!=this->end(); s++) xout << s->str()<<std::endl;

}

//void StringSet::insert(String& s)
//{
////    std::cout <<parallel_rank<< "StringSet::insert "<<s.str()<<std::endl;std::cout.flush();
////    xout <<parallel_rank<< "addressMap has "<<addressMap.size()<<" entries; size()="<<size()<<std::endl;
//  s.keygen(PartialWeightArray);
////  xout << "s.key="<<s.key<<", s.orbitals().size()="<<s.orbitals().size()<<std::endl;
//  if (addressMap.count(s.key())) {
////        std::cout <<parallel_rank<<" "<<size()<<" "<<addressMap.count(s.key)<< "StringSet::insert found existing"<<std::endl;std::cout.flush();
//    if (addressMap[s.key()] >= size()) throw std::logic_error("something wrong in StringSet reset");
//    at(addressMap[s.key()]) = s;
//  } else {
////        if (s.orbitals().size()==0) std::cout <<parallel_rank<< "StringSet::insert found new"<<std::endl;std::cout.flush();
//    addressMap[s.key()]=size();
//    memory::vector<String>::push_back(s);
////      if (s.orbitals().size()==0) std::cout <<parallel_rank<<"StringSet::push_back " <<s <<" size()=" <<size()<<std::endl;std::cout.flush();
//  }
//  //std::cout <<parallel_rank<< "StringSet::insert finished "<<std::endl;std::cout.flush();
//}

std::string StringSet::str(int verbosity, unsigned int columns) const {
  std::ostringstream s;
  if (verbosity >= -1) {
    s << "StringSet size=" << size();
  }
  if (verbosity > 0) {
    for (const auto &i : *this)
      s << std::endl << i.str();
  }
  return s.str();
}

std::vector<gci::ExcitationSet> StringSet::allExcitations(StringSet &to, int annihilations, int creations) {
  std::vector<ExcitationSet> set;
  for (auto ff : *this) {
    set.emplace_back(ff, to, annihilations, creations);
  }
  return set;
}

std::vector<double> StringSet::occupationNumbers() {
  std::vector<double> result;
  if (!this->empty()) {
    String firstString = this->at(0);
    result.resize(this->size() * firstString.orbitalSpace->total(), (double) 0);
    int stringoffset = 0;
    for (auto &s : *this) {
      auto orbitals = s.orbitals();
      for (unsigned char &orbital : orbitals) {
        //                xout << "StringSet::occupationNumbers stringoffset="<<stringoffset<<" *i="<<*i<<std::endl;
        result[stringoffset + (orbital - 1) * this->size()] = (double) 1;
      }
      stringoffset++;
    }
  }
  return result;
}

std::vector<size_t> StringSet::index(const StringSet &set) const {
  std::vector<size_t> result;
  for (const auto &s : *this)
    result.push_back(s.index(set));
  return result;
}
