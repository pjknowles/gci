#include <algorithm>
#include <sstream>
#include <iostream>
#include <deque>
#include <queue>
#include <string>
#include <iomanip>
#include "Profiler.h"

Profiler::Profiler()
{
}
Profiler::Profiler(std::string name)
{
  reset(name);
}

void Profiler::reset(const std::string name)
{
  Name=name;
  stopall();
  results.clear();
  start("* Other");
}

void Profiler::start(const std::string name)
{
  struct times now=getTimes();
  if (! stack.empty())
    stack.top()+=now;
  struct times minusNow; minusNow.cpu=-now.cpu; minusNow.wall=-now.wall; minusNow.name=name; minusNow.operations=0;
  stack.push(minusNow);
}

#include <assert.h>
void Profiler::stop(const std::string name, long operations)
{
//  if (operations>0) std::cout << "Profiler::stop "<<stack.top().name<<":"<<name<<" operations="<<operations<<std::endl;
  assert(name=="" || name == stack.top().name);
  struct times now=getTimes();now.operations=operations;
//  std::cout << "stack.top().operations="<<stack.top().operations<<std::endl;
  stack.top()+=now;
//  if (operations>0) std::cout << "stack.top().operations="<<stack.top().operations<<std::endl;
  results[stack.top().name] += stack.top();
  results[stack.top().name].calls++;
  stack.pop();
//  std::cout <<"now="<<now.operations<<std::endl;
//  if (stack.size()>=1) std::cout<<"stop before subtracting from top of stack " << stack.top().name << " " <<stack.top().operations <<std::endl;
  if (! stack.empty()) stack.top()-=now;
}

void Profiler::declare(const std::string name)
{
  if (results.count(name)==0) {
  struct times tt; tt.cpu=0; tt.wall=0; tt.name=name; tt.operations=0; tt.calls=0;
  results[name] = tt;
  }
}

void Profiler::stopall()
{
  while (! stack.empty()) stop();
}

#include <cmath>
#ifdef GCI_PARALLEL
#define HAVE_PPIDD
#endif
#ifdef HAVE_PPIDD
extern "C" {
#include "ppidd_c.h"
}
#endif
#ifdef MOLPRO
#include "mpp/CxMpp.h"
#include "cic/ItfMpp.h"
itf::FMppInt interface;
//extern "C" {
//int64_t get_iprocs_cxx_();
//}
#endif
std::string Profiler::str(const int verbosity, const int precision)
{
  if (verbosity<0) return "";
  stopall();
  resultMap results=this->results; // local copy that we can sum globally
  while(results.erase(""));
  for (resultMap::iterator s=results.begin(); s!=results.end(); ++s) {
#ifdef GCI_PARALLEL
    int64_t type=1, len=1;
    char* opm=strdup("max");
    PPIDD_Gsum(&type,&((*s).second.wall),&len,opm);
    char* op=strdup("+");
    PPIDD_Gsum(&type,&((*s).second.cpu),&len,op);
    int64_t value=(int64_t)(*s).second.calls; type=0;
    PPIDD_Gsum(&type,&value,&len,op);
    (*s).second.calls=(int)value;
    value=(int64_t)(*s).second.operations; type=0;
    PPIDD_Gsum(&type,&value,&len,op);
    (*s).second.operations=(long)value;
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
#endif
#endif
  }
  typedef std::pair<std::string,Profiler::times> data_t;
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#else
  std::priority_queue<data_t, std::deque<data_t>, compareTimes<data_t>  > q(results.begin(),results.end());
#endif
  std::stringstream ss;
  size_t maxWidth=0;
  long maxOperations=0;
  Profiler::times totalTimes;totalTimes.operations=0;
  for (resultMap::const_iterator s=results.begin(); s!=results.end(); ++s) {
      if ((*s).second.operations > maxOperations) maxOperations=(*s).second.operations;
      if ((*s).first.size() > maxWidth) maxWidth=(*s).first.size();
      totalTimes += (*s).second;
  }
  totalTimes.calls=1;
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#else
  q.push(data_t("* TOTAL",totalTimes));
#endif
  ss << "Profiler "<<Name<<std::endl;
  std::vector<std::string> prefixes;
  prefixes.push_back(""); prefixes.push_back("k"); prefixes.push_back("M"); prefixes.push_back("G");
  prefixes.push_back("T"); prefixes.push_back("P"); prefixes.push_back("E"); prefixes.push_back("Z"); prefixes.push_back("Y");
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#else
  while (! q.empty()) {
    ss.precision(precision);
    ss <<std::right <<std::setw(maxWidth) << q.top().first <<": calls="<<q.top().second.calls<<", cpu="<<std::fixed<<q.top().second.cpu<<", wall="<<q.top().second.wall;
    double ops=q.top().second.operations;
    double wall=q.top().second.wall;
    if (ops>(double)0 && wall>(double)0) {
      ops /= wall;
      int shifter = ops > 1 ? (int)(log10(ops)/3) : 0 ; shifter = shifter >= (int) prefixes.size() ? (int) prefixes.size()-1 : shifter;  ops *= pow((double)10, -shifter*3);
      ss<<", "<<ops<<" "<<prefixes[shifter]<<"op/s";
    }
      ss <<std::endl;
    q.pop();
  }
#endif
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, Profiler & obj)
{
  return os << obj.str();
}

#include <time.h>
#include <sys/time.h>
struct Profiler::times Profiler::getTimes()
{
  struct Profiler::times result;
  result.operations=0;
  result.cpu=(double)clock()/CLOCKS_PER_SEC;
  struct timeval time;
  result.wall=(double)0;
  if (!gettimeofday(&time,NULL))
    result.wall = (double)time.tv_sec + (double)time.tv_usec * .000001;
  return result;
}

struct Profiler::times& Profiler::times::operator+=( const struct Profiler::times &w2)
{
  cpu += w2.cpu;
  wall += w2.wall;
  operations += w2.operations;
  return *this;
}

struct Profiler::times Profiler::times::operator+(const struct Profiler::times &w2)
{
  struct Profiler::times result=*this;
  result += w2;
  return result;
}

struct Profiler::times& Profiler::times::operator-=( const struct Profiler::times &w2)
{
  cpu -= w2.cpu;
  wall -= w2.wall;
  operations -= w2.operations;
  return *this;
}

struct Profiler::times Profiler::times::operator-(const struct Profiler::times &w2)
{
  struct Profiler::times result=*this;
  result -= w2;
  return result;
}

// C binding
extern "C" {
#include "ProfilerC.h"
#include <stdlib.h>
#include <string.h>
void* profilerNew(char* name) { return new Profiler(name); }
void profilerReset(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->reset(std::string(name)); }
void profilerStart(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->start(std::string(name)); }
void profilerDeclare(void* profiler, char* name) { Profiler* obj=(Profiler*)profiler; obj->declare(std::string(name)); }
void profilerStop(void* profiler, char* name, long operations) { Profiler* obj=(Profiler*)profiler; obj->stop(std::string(name),operations); }
char* profilerStr(void* profiler) { Profiler* obj=(Profiler*)profiler; char* result = (char*)malloc(obj->str().size()+1); strcpy(result, obj->str().c_str()); return result; }
  void profilerStrSubroutine(void*profiler, char* result, int maxResult) { strncpy(result, profilerStr(profiler),maxResult-1);}
}
