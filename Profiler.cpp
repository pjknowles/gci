#include <algorithm>
#include <sstream>
#include <iostream>
#include <deque>
#include <queue>
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
//  std::cout << "Profiler::start "<<name<<std::endl;
  struct times now=getTimes();
  if (! stack.empty())
    stack.top()+=now;
//  if (stack.size()==1) std::cout<<"adjusted top of stack " << stack.top().name << " " <<stack.top().wall <<std::endl;
  struct times minusNow; minusNow.cpu=-now.cpu; minusNow.wall=-now.wall; minusNow.name=name;
  stack.push(minusNow);
//  if (stack.size()==1) std::cout<<"starting stack " << stack.top().name << " " <<stack.top().wall <<std::endl;
}

#include <assert.h>
void Profiler::stop(const std::string name, size_t operations)
{
//  std::cout << "Profiler::stop "<<stack.top().name<<":"<<name<<std::endl;
  assert(name=="" || name == stack.top().name);
  struct times now=getTimes();now.operations=operations;
  stack.top()+=now;
  results[stack.top().name] += stack.top();
  results[stack.top().name].calls++;
//  if (stack.size()==1) {
//    std::cout<<"stop added to top of stack " << stack.top().name << " " <<stack.top().wall <<std::endl;
//    std::cout<<"results now "<<results[stack.top().name].wall <<std::endl;
//  }
  stack.pop();
  if (! stack.empty()) stack.top()-=now;
//  if (stack.size()==1) std::cout<<"stop subtracted from top of stack " << stack.top().name << " " <<stack.top().wall <<std::endl;
}

void Profiler::stopall()
{
  while (! stack.empty()) stop();
}

std::string Profiler::str(const int verbosity, const int precision)
{
  if (verbosity<0) return "";
  stopall();
  typedef std::pair<std::string,Profiler::times> data_t;
  std::priority_queue<data_t, std::deque<data_t>, compareTimes<data_t>  > q(results.begin(),results.end());
  std::stringstream ss;
  size_t maxWidth=0;
  size_t maxOperations;
  Profiler::times totalTimes;
  for (resultMap::const_iterator s=results.begin(); s!=results.end(); ++s) {
      if ((*s).second.operations > maxOperations) maxOperations=(*s).second.operations;
      if ((*s).first.size() > maxWidth) maxWidth=(*s).first.size();
      totalTimes += (*s).second;
  }
  totalTimes.calls=1;
  q.push(data_t("* TOTAL",totalTimes));
  ss << "Profiler "<<Name<<std::endl;
  while (! q.empty()) {
    ss.precision(precision);
    ss <<std::right <<std::setw(maxWidth) << q.top().first <<": calls="<<q.top().second.calls<<", cpu="<<std::fixed<<q.top().second.cpu<<", wall="<<q.top().second.wall;
    if (q.top().second.operations>0)
      ss<<", operations="<<q.top().second.operations;
      ss <<std::endl;
    q.pop();
  }
  return ss.str();
}

std::ostream& operator<<(std::ostream& os, Profiler & obj)
{
  return os << obj.str();
}

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
struct Profiler::times& Profiler::times::operator-=( const struct Profiler::times &w2)
{
  cpu -= w2.cpu;
  wall -= w2.wall;
  operations -= w2.operations;
  return *this;
}

struct Profiler::times operator+(const struct Profiler::times &w1, const struct Profiler::times &w2)
{
  struct Profiler::times result=w1;
  result += w2;
  return result;
}
struct Profiler::times operator-(const struct Profiler::times &w1, const struct Profiler::times &w2)
{
  struct Profiler::times result=w1;
  result -= w2;
  return result;
}
