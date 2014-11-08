#include "gciStringSet.h"
#include <iostream>
#include <sstream>

StringSet::StringSet() : std::vector<String>()
{
  //    xout <<"StringSet default constructor"<<std::endl;
}

StringSet::StringSet(String prototype, bool all, int sym) : std::vector<String>()
{
  //    xout <<"StringSet prototype constructor "<<all<<std::endl;
  // copy prototype
  proto = prototype;
  symmetry = sym;
  setupPartialWeightArray();
  if (all) complete(sym);
}

void StringSet::makekey(String &s)
{
  s.key=0;
  for (int k=0; k<(int)s.orbitals_.size(); k++)
    s.key+= PartialWeightArray[k][s.orbitals_[k]-1];
}

StringSet::StringSet(const StringSet &referenceSpace, int annihilations, int creations, int sym)
{
  addByOperators(referenceSpace, annihilations, creations, sym);
}

StringSet::StringSet(const std::vector<StringSet>& referenceSpaces, int annihilations, int creations, int sym)
{
  addByOperators(referenceSpaces, annihilations, creations, sym);
}

#include <string.h>
void StringSet::addByOperators(const std::vector<StringSet> &referenceSpaces, int annihilations, int creations, int sym)
{
  profiler.start("StringSet::addByOperators[]");
//  xout <<"referenceSpaces.size()="<<referenceSpaces.size()<<std::endl;
  size_t ntask=0;
  for (std::vector<StringSet>::const_iterator s=referenceSpaces.begin(); s!=referenceSpaces.end(); s++)
    ntask+=(*s).size();
//  xout << "parallel_rank="<<parallel_rank<<", ntask="<<ntask<<std::endl;
  DivideTasks(ntask);
  for (std::vector<StringSet>::const_iterator referenceSpace=referenceSpaces.begin(); referenceSpace != referenceSpaces.end(); referenceSpace++)
    addByOperators(*referenceSpace, annihilations, creations, sym, true);
  EndTasks();
  profiler.start("StringSet::addByOperators:distribute");
  std::vector<char> serialised;
  for (StringSet::const_iterator s=begin(); s!=end(); s++) {
    std::vector<char> serialised1=s->serialise();
    for (std::vector<char>::const_iterator c=serialised1.begin(); c!=serialised1.end();c++)
      serialised.push_back(*c);
  }
//  xout << "serialised "<<serialised.size()<<" bytes from "<<size()<<" String objects="<<std::endl;
//  xout << "addressMap.size()"<< addressMap.size()<<std::endl;
#ifdef GCI_MPI
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
  String refString = referenceSpaces[0][0];
  int bytestreamsize;
  if (parallel_rank>0) {
    int len=(int)size();
    MPI_Send(&len,(int) 1,MPI_INT,0,0,MPI_COMM_COMPUTE);
    bytestreamsize=(int)serialised.size()/size();
    MPI_Send(&bytestreamsize,(int) 1,MPI_INT,0,1,MPI_COMM_COMPUTE);
    MPI_Send(&serialised[0],len*bytestreamsize,MPI_BYTE,0,2,MPI_COMM_COMPUTE);
    MPI_Bcast(&len,(int) 1,MPI_INT,0,MPI_COMM_COMPUTE);
    MPI_Barrier(MPI_COMM_COMPUTE);
//    xout << "slave after receving broadcast len "<<len<<std::endl;
    serialised.resize(len*bytestreamsize);
    clear();
    MPI_Bcast(&serialised[0],len*bytestreamsize,MPI_BYTE,0,MPI_COMM_COMPUTE);
    MPI_Barrier(MPI_COMM_COMPUTE);
    for (size_t k=0; k<len; k++) {
        std::vector<char> s(bytestreamsize); memcpy(&s[0],&serialised[k*bytestreamsize],bytestreamsize);
        xout << "slave construct string"<<std::endl;
        String ss(s,&refString);
        xout << "slave insert string "<<ss.str()<<std::endl;
        insert(ss);
        xout << "slave inserted string"<<std::endl;
    }
  } else {
      xout <<"master size()="<<size()<<", bytestreamsize="<<serialised.size()/size()<<", serialised.size()"<<serialised.size()<<std::endl;
    for (int iproc=1; iproc < parallel_size;iproc++) {
      int len;
      MPI_Status status;
      MPI_Recv(&len,(int) 1,MPI_INT,iproc,0,MPI_COMM_COMPUTE,&status);
      MPI_Recv(&bytestreamsize,(int) 1,MPI_INT,iproc,1,MPI_COMM_COMPUTE,&status);
      xout <<"received len="<<len<<", bytestreamsize="<<bytestreamsize<<std::endl;
      serialised.resize((size_t)len*bytestreamsize);
      MPI_Recv(&serialised[0],(int) len*bytestreamsize,MPI_BYTE,iproc,2,MPI_COMM_COMPUTE,&status);
      for (size_t k=0; k<len; k++) {
        std::vector<char> s(bytestreamsize); memcpy(&s[0],&serialised[k*bytestreamsize],bytestreamsize);
        xout << "master construct string"<<std::endl;
        String ss(s,&refString);
        xout <<"master before insert"<<std::endl;
        insert(ss);
        xout <<"master after insert"<<std::endl;
      }
      serialised.clear();
      for (StringSet::const_iterator s=begin(); s!=end(); s++) {
        std::vector<char> serialised1=s->serialise();
        for (std::vector<char>::const_iterator c=serialised1.begin(); c!=serialised1.end();c++)
          serialised.push_back(*c);
      }
    }
      int len=(int)size();
        xout <<"master after serialising global list"<<std::endl;
      MPI_Bcast(&len,(int) 1,MPI_INT,0,MPI_COMM_COMPUTE);
    MPI_Barrier(MPI_COMM_COMPUTE);
        xout <<"master after broadcasting len"<<std::endl;
      MPI_Bcast(&serialised[0],len*bytestreamsize,MPI_BYTE,0,MPI_COMM_COMPUTE);
    MPI_Barrier(MPI_COMM_COMPUTE);
        xout <<"master after broadcasting global list"<<std::endl;
  }
    MPI_Barrier(MPI_COMM_COMPUTE);
    xout << "Reached end of forked code rank="<<parallel_rank<<std::endl;
    MPI_Barrier(MPI_COMM_COMPUTE);
#ifdef GCI_MPI
//  for (size_t i=0; i<parallel_size; i++){
//  MPI_Barrier(MPI_COMM_COMPUTE);
//  if (i!=parallel_rank)continue;
//  xout << "Barrier crossed parallel_size="<<parallel_size<<", parallel_rank="<<parallel_rank<<std::endl;
//  xout << "end of addByOperators size()="<<size()<<std::endl;
//  xout << "merged addressMap rank="<<parallel_rank<<", size="<<addressMap.size()<<std::endl;
//  for (std::map<size_t,size_t> ::const_iterator a=addressMap.begin(); a!=addressMap.end(); a++)
//    xout <<parallel_rank<<": "<<(*a).first<<","<<(*a).second<<std::endl;
//  }
#endif
#endif
  profiler.stop("StringSet::addByOperators:distribute");
  profiler.stop("StringSet::addByOperators[]");
}

void StringSet::addByOperators(const StringSet &referenceSpace, int annihilations, int creations, int sym, bool parallel)
{
  profiler.start("addByOperators");
  size_t count=0;
  size_t countall=0;
  bool first=size()==0;
  symmetry = sym;
  //    xout << "in StringSet creator, referenceSpace="<<referenceSpace.str(5)<<std::endl;
  if ((int) referenceSpace.proto.nelec + creations - annihilations < 0
      || (int) referenceSpace.proto.nelec + creations - annihilations > (int) referenceSpace.proto.orbitals_.size())
    return; // null space because not enough electrons or holes left
  int symexc = (referenceSpace.symmetry>=0 && sym >=0) ? referenceSpace.symmetry ^ sym : -1 ; // use symmetry if we can
  for (StringSet::const_iterator s = referenceSpace.begin(); s != referenceSpace.end(); s++) {
    countall++;
    if (! NextTask()) continue;
    if (parallel && ! NextTask()) continue;
//    if (parallel && countall%parallel_size != parallel_rank) continue;
    count++;
    String from = *s;
    //        xout << "from="<<from.str(5)<<std::endl;
    if (annihilations + creations ==1) {
      for (int i=0; i<(int)from.orbitalSpace->orbital_symmetries.size(); i++) {
        if (from.orbitalSpace->orbital_symmetries[i]==(unsigned int)symexc || symexc==-1) {
          String tt = from;
          int phase = (annihilations > 0) ? tt.destroy(i+1) : tt.create(i+1);
          if (phase) {
            if (first) {
              proto = tt; setupPartialWeightArray(); insert(tt); proto = tt; first = false;
            } else
              insert(tt);
          }
        }
      }
    }
    else if (annihilations+creations==2) {
      for (int j=0; j<(int)from.orbitalSpace->orbital_symmetries.size(); j++) {
        String a = from;
        int phasea = (annihilations > 0) ?  a.destroy(j+1) : a.create(j+1);
        if (phasea) {
          for (int i=0; i<(int)from.orbitalSpace->orbital_symmetries.size(); i++) {
            if (from.orbitalSpace->orbital_symmetries[i]==((unsigned int)symexc^from.orbitalSpace->orbital_symmetries[j]) || symexc==-1) {
              String tt = a;
              int phase = (annihilations > 1) ? phasea*tt.destroy(i+1) : phasea*tt.create(i+1);
              if (phase) {
                if (first) {
                  proto = tt; setupPartialWeightArray(); insert(tt); proto = tt; first = false;
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
  profiler.stop("addByOperators",count);
}

void StringSet::setupPartialWeightArray()
{ // set up partial weight array for addressing binomial distributions
  int nitem = proto.nelec;
  int nbox = proto.orbitalSpace->total();
  PartialWeightArray = std::vector< std::vector<int> > (nitem, std::vector<int>(nbox));
  for (int k=0;k<nitem;k++) {
    for (int l=0;l<nbox;l++)
      PartialWeightArray[k][l]=0;
    for (int l=k;l<nbox-nitem+k;l++)
      PartialWeightArray[k][l+1] = binomial_coefficient(nbox-l-1,nitem-k-1)+PartialWeightArray[k][l];
  }
  for (int k=0;k<nitem-1;k++) {
    for (int l=k;l<nbox-nitem+k+1;l++)
      PartialWeightArray[k][l] = PartialWeightArray[k][l] - PartialWeightArray[k+1][l+1];
  }
  if (nitem > 0)
    for (int l=nitem-1;l<nbox;l++)
      PartialWeightArray[nitem-1][l] = l-nitem+1;

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

void StringSet::complete(int sym)
{
  profiler.start("StringSet::complete");
  //    xout <<"StringSet::complete prototype"<<proto.printable(1)<<std::endl;
  String string(&proto);
  this->erase(this->begin(),this->end());
  if (string.first(proto.nelec,sym)) {
    //    xout <<"StringSet::complete symmetry="<<sym<<" first String: "<<string.printable()<<std::endl;
    do {
      //        xout << "in StringSet::complete about to push_back " << string.printable(1) <<std::endl;
      this->insert(string);
    } while (string.next(sym));
  }
  //    xout << "in StringSet::complete final list: " <<std::endl ;
  //    for (iterator s=this->begin(); s!=this->end(); s++) xout << s->printable()<<std::endl;

  profiler.stop("StringSet::complete");
}

void StringSet::insert(String& s)
{
//  xout << "StringSet::insert "<<s.str()<<std::endl;
  s.key=0;
  for (int k=0; k<(int)s.orbitals_.size(); k++)
    s.key+= PartialWeightArray[k][s.orbitals_[k]-1];
  if (addressMap.count(s.key)) {
//    xout << "StringSet::insert found existing"<<std::endl;
    if (addressMap[s.key] >= size()) throw "something wrong in StringSet reset";
    at(addressMap[s.key]) = s;
  } else {
//    xout << "StringSet::insert found new"<<std::endl;
    addressMap[s.key]=size();
//            xout <<"StringSet::push_back " <<s <<" size()=" <<size()<<std::endl;
    std::vector<String>::push_back(s);
  }
//  xout << "StringSet::insert finished "<<std::endl;
}

std::string StringSet::str(int verbosity) const
{
  std::ostringstream s;
  if (verbosity >= -1) {
    s << "StringSet size=" << size();
  }
  if (verbosity >0 )
  {
    for (StringSet::const_iterator i=begin(); i!=end(); i++)
      s << std::endl << i->str();
  }
  return s.str();
}

std::vector<ExcitationSet> StringSet::allExcitations(StringSet &to, int annihilations, int creations)
{
  std::vector<ExcitationSet> set;
  for (iterator f=begin(); f!=end(); f++) {
    String ff = *f;
    set.push_back(ExcitationSet(ff,to,annihilations,creations));
  }
  return set;
}

std::vector<double> StringSet::occupationNumbers()
{
  std::vector<double> result;
  if (this->size()) {
    String firstString=this->at(0);
    result.resize(this->size()* firstString.orbitalSpace->total(), (double) 0);
    int stringoffset=0;
    for (StringSet::iterator s=this->begin(); s!=this->end(); s++)
    {
      std::vector<unsigned int> orbitals = s->orbitals();
      for (std::vector<unsigned int>::iterator i=orbitals.begin(); i !=orbitals.end(); i++) {
        //                xout << "StringSet::occupationNumbers stringoffset="<<stringoffset<<" *i="<<*i<<std::endl;
        result[stringoffset+(*i-1)*this->size()]=(double)1;
      }
      stringoffset++;
    }
  }
  return result;
}

std::vector<size_t> StringSet::index(const StringSet &set) const
{
  std::vector<size_t> result;
  for (const_iterator s=begin(); s!=end(); s++)
    result.push_back(s->index(set));
  return result;
}
