#ifdef MOLPRO
#include "common/molpro_config.h"
#else
#define _GNU_SOURCE
#endif
#include <algorithm>
#include <sstream>
#include <iostream>
#include <deque>
#include <queue>
#include <string>
#include <string.h>
#include <iomanip>
#include "Profiler.h"
#include "memory.h"

Profiler::Profiler()
{
}

Profiler::Profiler(std::string name, const int level)
{
  reset(name);
  active(level);
}

void Profiler::reset(const std::string name)
{
  Name=name;
  stopall();
  results.clear();
  active();
  level=0;
  start("TOP");
}

void Profiler::active(const int level, const int stopPrint)
{
    activeLevel=level;
    stopPrint_=stopPrint;
}

#include <assert.h>
void Profiler::start(const std::string name)
{
  level++;
  if (level>activeLevel) return;
  assert(level==resourcesStack.size()+1);
  struct resources now=getResources();now.name=name;now.calls=1;
  if (! resourcesStack.empty())
    totalise(now,0,0);
#ifdef MEMORY_H
  // memory accounting:
  // the statistic for a segment is the maximum used, ie memory_used(1) minus the actual start memory, ie memory_used(0)
  // memory_reset_maximum_stack()  is used to reset the memory manager's notion of high water
  // memoryStack[01] are used to store the values from memory_used in start()
  // in start(),
  //  (a) remember current memory in memoryStack[01]
  // in totalise(),
  //  (a) set the my maximum memory to max(existing,maximum memory used - memoryStack0[me]);
  // in stop(),
  //  (a) reset maximum stack to memoryStack1[me]
  // in accumulate(), add max(children) to parent
  //if (! memoryStack0.empty()) {
    //resourcesStack.back().stack = std::max((int64_t)memory_used(1)-memoryStack0.back(),resourcesStack.back().stack);
  //}
  memoryStack0.push_back(memory_used(0));
  memoryStack1.push_back(memory_used(1));
#endif
  resourcesStack.push_back(now);
  startResources.push_back(now);
//    std::cout <<std::endl<< "start now.wall="<<now.wall<<std::endl;
}

void Profiler::totalise(const struct resources now, const long operations, const int calls)
{
  resources diff=now;
  diff-=resourcesStack.back();
  diff.name=resourcesStack.back().name;
  diff.operations=operations;
#ifdef MEMORY_H
  diff.stack=memory_used(1)-memoryStack0.back();
#endif
  std::string key;
  for(std::vector<resources>::const_reverse_iterator r=resourcesStack.rbegin(); r!= resourcesStack.rend(); r++) key=r->name+":"+key;
  key.erase(key.end()-1,key.end());
  diff.name=key;
  results[key] += diff;
  results[key].calls += calls;
}

void Profiler::stop(const std::string name, long operations)
{
  level--;
  if (level > 0 && level>=activeLevel) return;
  assert(level==resourcesStack.size()-1);
  assert(name=="" || name == resourcesStack.back().name);
  struct resources now=getResources();now.operations=operations;
  totalise(now,operations,1);

  if (stopPrint_>-1) {
    struct resources diff=now;
//    std::cout <<std::endl<< "now.wall="<<now.wall<<std::endl;
//    std::cout << "then.wall="<<startResources.back().wall<<std::endl;
    diff-=startResources.back();
    diff.name="";
    for(std::vector<resources>::const_reverse_iterator r=resourcesStack.rbegin(); r!= resourcesStack.rend(); r++) diff.name=r->name+":"+diff.name;
    diff.name.erase(diff.name.end()-1,diff.name.end());
    std::cout <<"Profiler \""<<Name<<"\" stop: "<< diff.str(0,stopPrint_,false,3,"All") <<std::endl;
  }

#ifdef MEMORY_H
  memoryStack0.pop_back();
  memoryStack1.pop_back();
  memory_reset_maximum_stack(memoryStack1.back());
#endif
  resourcesStack.pop_back();
  startResources.pop_back();
  if (! resourcesStack.empty()) {now.name=resourcesStack.back().name; resourcesStack.back()=now;}
}

void Profiler::declare(const std::string name)
{
  if (results.count(name)==0) {
  struct resources tt; tt.cpu=0; tt.wall=0; tt.name=name; tt.operations=0; tt.calls=0;
  results[name] = tt;
  }
}

void Profiler::stopall()
{
  while (! resourcesStack.empty()) stop();
}

#include <cmath>
#ifdef GCI_PARALLEL
#include "mpi.h"
#endif
#ifdef MOLPRO
#include "mpp/CxMpp.h"
#include "cic/ItfMpp.h"
itf::FMppInt interface(itf::FMppInt::MPP_NeedSharedFs|itf::FMppInt::MPP_GlobalDeclaration);
//extern "C" {
//int64_t get_iprocs_cxx_();
//}
#endif
Profiler::resultMap Profiler::totals() const
{
  Profiler thiscopy=*this; // take a copy so that we can run stopall yet be const, and so that we can sum globally
  thiscopy.stopall();
  while(thiscopy.results.erase(""));
  for (resultMap::iterator s=thiscopy.results.begin(); s!=thiscopy.results.end(); ++s) {
#ifdef GCI_PARALLEL
      int type=1, len=1;
      double val=s->second.wall;
      MPI_Allreduce(&val,&((*s).second.wall),len,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      val=s->second.cpu;
      MPI_Allreduce(&val,&((*s).second.cpu),len,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      int calls=s->second.calls;
      MPI_Allreduce(&calls,&((*s).second.calls),len,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      long operations=s->second.operations;
      MPI_Allreduce(&operations,&((*s).second.operations),len,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
      int64_t stack = s->second.stack;
      MPI_Allreduce(&stack,&((*s).second.stack),len,MPI_LONG_LONG_INT,MPI_MAX,MPI_COMM_WORLD);
#else
#ifdef MOLPRO
    // only '+' works in Molpro runtime
    //    interface.GlobalSum(&((*s).second.wall),(std::size_t)1,(uint)0,(const char*) "max");
    interface.GlobalSum(&((*s).second.cpu),(std::size_t)1);
    // Molpro interface presently only does doubles
    double value=(double)(*s).second.calls;
    interface.GlobalSum(&value,(std::size_t)1);
    (*s).second.calls=(int)value;
    value=(double)(*s).second.operations;
    interface.GlobalSum(&value,(std::size_t)1);
    (*s).second.operations=(long)value;
    // global max of (*s).second.stack
    value=(double)(*s).second.stack;
    interface.GlobalSum(&value,(std::size_t)1, (uint)0, "max");
    (*s).second.stack=(int64_t)value;
#endif
#endif
  }
  thiscopy.accumulate(thiscopy.results);
  return thiscopy.results;
}

std::string Profiler::resources::str(const int width, const int verbosity, const bool cumulative, const int precision, const std::string defaultName) const
{
  std::stringstream ss;
  std::vector<std::string> prefixes;
  prefixes.push_back(""); prefixes.push_back("k"); prefixes.push_back("M"); prefixes.push_back("G"); prefixes.push_back("T"); prefixes.push_back("P"); prefixes.push_back("E"); prefixes.push_back("Z"); prefixes.push_back("Y");
  const struct resources * r=this;
  std::string name=r->name;
  if (cumulative) r=(r->cumulative);
  if (name.find(":") == std::string::npos)
    if (defaultName!="")
      name=defaultName;
    else
      name = cumulative ? "All" : "(other)";
  else
    name.replace(0,name.find(":")+1,"");
  size_t wid = width > 0 ?  width : name.size();
  ss.precision(precision);
  ss <<std::right <<std::setw(wid) << name <<":";
  if (r->calls > 0) ss<<" calls="<<r->calls<<",";
  ss<<" cpu="<<std::fixed<<r->cpu<<","<<" wall="<<r->wall;
  double ops=r->operations;
  double wall=r->wall;
  if (ops>(double)0 && wall>(double)0) {
    ops /= wall;
    int shifter = ops > 1 ? (int)(log10(ops)/3) : 0 ; shifter = shifter >= (int) prefixes.size() ? (int) prefixes.size()-1 : shifter;  ops *= pow((double)10, -shifter*3);
    ss<<", "<<ops<<" "<<prefixes[shifter]<<"op/s";
  }
  size_t stack=r->stack;
  if (stack > (size_t)0) {
    ss<<", stack="<<stack;
  }
  return ss.str();
}

std::string Profiler::str(const int verbosity, const bool cumulative, const int precision) const
{
  if (verbosity<0) return "";
  resultMap localResults=totals();
  typedef std::pair<std::string,Profiler::resources> data_t;
  std::priority_queue<data_t, std::deque<data_t>, compareResources<data_t>  > q(localResults.begin(),localResults.end());
  std::stringstream ss;
  size_t maxWidth=0;
  for (resultMap::const_iterator s=localResults.begin(); s!=localResults.end(); ++s)
    if ((*s).first.size() > maxWidth) maxWidth=(*s).first.size();
  maxWidth-=localResults.begin()->first.size()+1; // assumes the first node is the top level
  if (maxWidth < 7) maxWidth=7;
  ss << "Profiler \""<<Name<<"\""; if(cumulative) ss<<" (cumulative)"; if (activeLevel < INT_MAX) ss <<" to depth "<<activeLevel; ss <<std::endl;
  for ( ; ! q.empty(); q.pop())
      ss << q.top().second.str(maxWidth,verbosity,cumulative,precision) <<std::endl;
  return ss.str();
}

void Profiler::accumulate(resultMap& results)
{
  for (resultMap::iterator r=results.begin(); r!=results.end(); ++r) r->second.name=r->first;
  for (resultMap::iterator parent=results.begin(); parent!=results.end(); ++parent) {
      parent->second.cumulative = new Profiler::resources;
      *parent->second.cumulative-=*parent->second.cumulative;
      // nb 'child' includes the parent itself
    for (resultMap::iterator child=results.begin(); child!=results.end(); ++child) {
      if (parent->first.size() <= child->first.size() && parent->first == child->first.substr(0,parent->first.size())) {
        *parent->second.cumulative += child->second;
	if (parent->first != child->first)
	  parent->second.cumulative->stack=std::max(parent->second.cumulative->stack,parent->second.stack+child->second.stack);
        parent->second.cumulative->calls = parent->second.calls;
      }
    }
  }
}

std::ostream& operator<<(std::ostream& os, Profiler & obj)
{
  return os << obj.str();
}

#include <time.h>
#include <sys/time.h>
static int init=1;
static double wallbase;
struct Profiler::resources Profiler::getResources()
{
  struct Profiler::resources result;
  result.operations=0;
  result.calls=0;
  result.cpu=(double)clock()/CLOCKS_PER_SEC;
  struct timeval time;
  result.wall=(double)0;
  if (!gettimeofday(&time,NULL)) {
    result.wall = (double)time.tv_sec + (double)time.tv_usec * .000001;
    if (init) wallbase=result.wall;
    init=0;
    result.wall-=wallbase;
  }
#ifdef MEMORY_H
  result.stack=(size_t) memory_used(1);
#else
  result.stack=0;
  result.cumulative=NULL;
#endif
  return result;
}

/*!
 * \brief Profiler::resources::operator += add another object to this one
 * \param w2 object to add
 * \return a copy of the object
 */
struct Profiler::resources& Profiler::resources::operator+=( const struct Profiler::resources &w2)
{
  cpu += w2.cpu;
  wall += w2.wall;
  operations += w2.operations;
  if (w2.stack > stack) stack = w2.stack; // choose maximum stack of the two objects
  return *this;
}

struct Profiler::resources Profiler::resources::operator+(const struct Profiler::resources &w2)
{
  struct Profiler::resources result=*this;
  result += w2;
  return result;
}

struct Profiler::resources& Profiler::resources::operator-=( const struct Profiler::resources &w2)
{
  cpu -= w2.cpu;
  wall -= w2.wall;
  operations -= w2.operations;
  stack -= w2.stack;
  return *this;
}

struct Profiler::resources Profiler::resources::operator-(const struct Profiler::resources &w2)
{
  struct Profiler::resources result=*this;
  result -= w2;
  return result;
}

// C binding
extern "C" {
#include <stdlib.h>
#include <string.h>
void* profilerNew(char* name) { return new Profiler(name); }
void profilerReset(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->reset(std::string(name)); }
void profilerActive(void* profiler, int level, int stopPrint) { Profiler* obj=(Profiler*)profiler; obj->active(level,stopPrint); }
void profilerStart(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->start(std::string(name)); }
void profilerDeclare(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->declare(std::string(name)); }
void profilerStop(void* profiler, char* name, long operations) { Profiler* obj=(Profiler*)profiler; obj->stop(std::string(name),operations); }
char* profilerStr(void* profiler, int verbosity, int cumulative, int precision) { Profiler* obj=(Profiler*)profiler; std::string res = obj->str(verbosity,bool(cumulative), precision); char* result = (char*)malloc(res.size()+1); strcpy(result, res.c_str()); return result; }
void profilerStrSubroutine(void*profiler, char* result, int maxResult, int verbosity, int cumulative, int precision) { strncpy(result, profilerStr(profiler, verbosity, cumulative, precision),maxResult-1);}
}
